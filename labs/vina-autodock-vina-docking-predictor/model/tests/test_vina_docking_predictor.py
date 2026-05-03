from __future__ import annotations

import io
import json
from pathlib import Path
import subprocess
import time

import pytest
from biosim.signals import (AcceptedSignalProfile, ArraySignal, BioSignal, EventSignal, RecordSignal, ScalarSignal, SignalSpec)
import yaml


def _set_required_inputs(
    module,
    BioSignal,
    *,
    receptor_pdbqt_path: str | None = None,
    ligand_pdbqt_path: str | None = None,
    run_options: dict | None = None,
):
    signals = {}
    if receptor_pdbqt_path is not None:
        signals["receptor_pdbqt_path"] = _make_signal(source="test", name="receptor_pdbqt_path", value=receptor_pdbqt_path, emitted_at=0.0, spec=None)
    if ligand_pdbqt_path is not None:
        signals["ligand_pdbqt_path"] = _make_signal(source="test", name="ligand_pdbqt_path", value=ligand_pdbqt_path, emitted_at=0.0, spec=None)
    if run_options is not None:
        signals["run_options"] = _make_signal(source="test", name="run_options", value=run_options, emitted_at=0.0, spec=None)
    module.set_inputs(signals)


def _fake_vina_output() -> str:
    return """MODEL 1
REMARK VINA RESULT:      -13.23      0.000      0.000
ROOT
ATOM      1  N   UNL     1      16.600  51.810  14.798  1.00  0.00    -0.322 N
ATOM      2  C   UNL     1      15.629  51.747  15.784  1.00  0.00     0.255 C
ENDROOT
ENDMDL
MODEL 2
REMARK VINA RESULT:      -11.29      0.986      1.681
ROOT
ATOM      1  N   UNL     1      16.700  51.900  14.900  1.00  0.00    -0.322 N
ATOM      2  C   UNL     1      15.700  51.800  15.900  1.00  0.00     0.255 C
ENDROOT
ENDMDL
"""


def test_instantiation(biosim, tmp_path):
    from src.vina_docking_predictor import VinaDockingPredictor

    module = VinaDockingPredictor(work_dir=str(tmp_path))
    assert module.integration_step > 0
    assert module.runtime_mode == "managed"
    assert set(module.inputs()) == {"receptor_pdbqt_path", "ligand_pdbqt_path", "run_options"}
    assert set(module.outputs()) == {"pose_summary", "docking_summary", "structure_artifacts", "run_metadata"}
    assert module.vina_version == "1.2.7"


def test_missing_inputs_surface_error_metadata(biosim, tmp_path):
    from src.vina_docking_predictor import VinaDockingPredictor

    module = VinaDockingPredictor(work_dir=str(tmp_path))
    module.advance_window(0.0, 0.1)

    outputs = module.get_outputs()
    assert _signal_value(outputs["run_metadata"])["status"] == "error"
    assert "receptor_pdbqt_path" in _signal_value(outputs["run_metadata"])["error"]
    assert module.visualize() is None


def test_model_relative_path_resolution_uses_checked_in_asset(biosim):
    from src.vina_docking_predictor import VinaDockingPredictor

    module = VinaDockingPredictor()
    resolved = module._resolve_required_input_path(
        "data/1iep/1iep_receptor.pdbqt",
        input_name="receptor_pdbqt_path",
    )
    assert Path(resolved).is_file()
    assert Path(resolved).name == "1iep_receptor.pdbqt"


def test_run_options_validation_requires_box_fields(biosim, tmp_path):
    from src.vina_docking_predictor import VinaDockingPredictor
    from biosim.signals import BioSignal

    module = VinaDockingPredictor(work_dir=str(tmp_path))
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={"exhaustiveness": 8},
    )
    module.advance_window(0.0, 0.1)

    metadata = _signal_value(module.get_outputs()["run_metadata"])
    assert metadata["status"] == "error"
    assert "run_options.box_center" in metadata["error"]


def test_run_options_validation_rejects_unknown_keys(biosim, tmp_path):
    from src.vina_docking_predictor import VinaDockingPredictor
    from biosim.signals import BioSignal

    module = VinaDockingPredictor(work_dir=str(tmp_path))
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={
            "box_center": [15.19, 53.903, 16.917],
            "box_size": [20.0, 20.0, 20.0],
            "not_supported": 1,
        },
    )
    module.advance_window(0.0, 0.1)

    metadata = _signal_value(module.get_outputs()["run_metadata"])
    assert metadata["status"] == "error"
    assert "unsupported run_options key" in metadata["error"]


def test_platform_tag_resolution_covers_supported_pairs(biosim, monkeypatch):
    from src.vina_docking_predictor import VinaDockingPredictor

    module = VinaDockingPredictor()
    monkeypatch.setattr("src.vina_docking_predictor.platform.system", lambda: "Darwin")
    monkeypatch.setattr("src.vina_docking_predictor.platform.machine", lambda: "arm64")
    assert module._platform_tag() == "mac_aarch64"

    monkeypatch.setattr("src.vina_docking_predictor.platform.system", lambda: "Linux")
    monkeypatch.setattr("src.vina_docking_predictor.platform.machine", lambda: "x86_64")
    assert module._platform_tag() == "linux_x86_64"


def test_build_command_supports_scoring_and_seed(biosim):
    from src.vina_docking_predictor import VinaDockingPredictor

    module = VinaDockingPredictor()
    command = module._build_command(
        vina_executable="/tmp/vina",
        receptor_path="/tmp/receptor.pdbqt",
        ligand_path="/tmp/ligand.pdbqt",
        config_path=Path("/tmp/config.txt"),
        output_path=Path("/tmp/output.pdbqt"),
        options={
            "box_center": [1.0, 2.0, 3.0],
            "box_size": [4.0, 5.0, 6.0],
            "exhaustiveness": 32,
            "n_poses": 5,
            "energy_range": 4.5,
            "cpu": 2,
            "seed": 7,
            "scoring": "vinardo",
        },
    )
    assert "--scoring" in command
    assert "vinardo" in command
    assert "--seed" in command
    assert "7" in command
    assert "--num_modes" in command
    assert "5" in command


def test_managed_runtime_bootstraps_and_parses_outputs(biosim, tmp_path, monkeypatch):
    from src.vina_docking_predictor import VinaDockingPredictor
    from biosim.signals import BioSignal

    downloads: list[tuple[str, Path]] = []

    def fake_download(url: str, target: Path) -> None:
        downloads.append((url, target))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        target.chmod(0o755)

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if command[0].endswith("vina_split"):
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="fallback")
        raise AssertionError(f"Unexpected command: {command}")

    def fake_run_command(*, command, cwd, timeout, env, phase, start_message, heartbeat_message, completion_message):  # noqa: ARG001
        output_path = Path(command[command.index("--out") + 1])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(_fake_vina_output(), encoding="utf-8")
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="Scoring function : vina\nPerforming docking (random seed: 12345) ...\n",
            stderr="",
        )

    monkeypatch.setattr(VinaDockingPredictor, "_download_binary", staticmethod(fake_download))
    monkeypatch.setattr(VinaDockingPredictor, "_run_command_with_live_logs", staticmethod(fake_run_command))
    monkeypatch.setattr(subprocess, "run", fake_run)

    module = VinaDockingPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "runtime"),
        cache_dir=str(tmp_path / "cache"),
    )
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={
            "box_center": [15.19, 53.903, 16.917],
            "box_size": [20.0, 20.0, 20.0],
            "exhaustiveness": 32,
            "n_poses": 2,
            "cpu": 0,
            "scoring": "vina",
        },
    )
    module.advance_window(0.0, 0.5)

    outputs = module.get_outputs()
    metadata = _signal_value(outputs["run_metadata"])
    pose_summary = _signal_value(outputs["pose_summary"])
    docking_summary = _signal_value(outputs["docking_summary"])
    artifacts = _signal_value(outputs["structure_artifacts"])

    assert metadata["status"] == "completed"
    assert metadata["runtime_bootstrapped"] is True
    assert metadata["seed"] == 12345
    assert docking_summary["top_pose_affinity_kcal_mol"] == -13.23
    assert len(pose_summary) == 2
    assert pose_summary[0]["rank"] == 1
    assert Path(artifacts["docking_output_file"]).is_file()
    assert Path(artifacts["top_pose_file"]).is_file()
    assert Path(artifacts["top_complex_file"]).is_file()
    assert Path(artifacts["stdout_tail_file"]).is_file()
    assert Path(artifacts["stderr_tail_file"]).is_file()
    assert len(downloads) == 2

    visuals = module.visualize()
    assert visuals is not None
    assert visuals[0]["render"] == "structure3d"
    assert visuals[0]["data"]["format"] == "pdb"
    assert visuals[1]["render"] == "table"
    assert visuals[1]["data"]["columns"] == ["Rank", "Affinity", "RMSD l.b.", "RMSD u.b.", "Pose File"]


def test_managed_runtime_reuses_cached_binary(biosim, tmp_path, monkeypatch):
    from src.vina_docking_predictor import VinaDockingPredictor

    downloads: list[tuple[str, Path]] = []

    def fake_download(url: str, target: Path) -> None:
        downloads.append((url, target))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        target.chmod(0o755)

    monkeypatch.setattr(VinaDockingPredictor, "_download_binary", staticmethod(fake_download))

    module = VinaDockingPredictor(
        runtime_dir=str(tmp_path / "runtime"),
        cache_dir=str(tmp_path / "cache"),
    )
    metadata: dict[str, object] = {}
    module._ensure_managed_runtime(metadata)
    assert len(downloads) == 2

    metadata = {}
    module._ensure_managed_runtime(metadata)
    assert len(downloads) == 2


def test_run_command_with_live_logs_emits_progress_and_streams_lines(biosim, tmp_path, monkeypatch, capsys):
    from src.vina_docking_predictor import VinaDockingPredictor

    class FakePopen:
        def __init__(self, *args, **kwargs):  # noqa: D401, ARG002
            self.stdout = io.StringIO("stdout line 1\nstdout line 2\n")
            self.stderr = io.StringIO("stderr line 1\n")
            self.returncode = None
            self._polls = 0

        def poll(self):
            self._polls += 1
            if self._polls >= 4:
                self.returncode = 0
                return 0
            return None

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(subprocess, "Popen", FakePopen)

    module = VinaDockingPredictor(
        work_dir=str(tmp_path),
        progress_heartbeat_s=0.01,
    )
    completed = module._run_command_with_live_logs(
        command=["vina", "--help"],
        cwd=tmp_path,
        timeout=10.0,
        env=None,
        phase="docking",
        start_message="Starting AutoDock Vina docking",
        heartbeat_message="AutoDock Vina docking is still running",
        completion_message="AutoDock Vina docking finished",
    )

    captured = capsys.readouterr()
    assert completed.returncode == 0
    assert "stdout line 1" in captured.out
    assert "stderr line 1" in captured.err
    progress_events = [
        json.loads(line.removeprefix("BSIM_PROGRESS:"))
        for line in captured.out.splitlines()
        if line.startswith("BSIM_PROGRESS:")
    ]
    phases = {event["phase"] for event in progress_events}
    assert "docking" in phases
    assert any("still running" in event["message"].lower() for event in progress_events)


def test_subprocess_failure_surfaces_metadata(biosim, tmp_path, monkeypatch):
    from src.vina_docking_predictor import VinaDockingPredictor
    from biosim.signals import BioSignal

    def fake_download(url: str, target: Path) -> None:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        target.chmod(0o755)

    def fake_run_command(*, command, cwd, timeout, env, phase, start_message, heartbeat_message, completion_message):  # noqa: ARG001
        return subprocess.CompletedProcess(command, 3, stdout="", stderr="boom")

    monkeypatch.setattr(VinaDockingPredictor, "_download_binary", staticmethod(fake_download))
    monkeypatch.setattr(VinaDockingPredictor, "_run_command_with_live_logs", staticmethod(fake_run_command))

    module = VinaDockingPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "runtime"),
    )
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={
            "box_center": [15.19, 53.903, 16.917],
            "box_size": [20.0, 20.0, 20.0],
        },
    )
    module.advance_window(0.0, 0.5)

    metadata = _signal_value(module.get_outputs()["run_metadata"])
    assert metadata["status"] == "error"
    assert metadata["returncode"] == 3
    assert "non-zero" in metadata["error"]


def test_repeat_advance_does_not_rerun_until_reset(biosim, tmp_path, monkeypatch):
    from src.vina_docking_predictor import VinaDockingPredictor
    from biosim.signals import BioSignal

    calls = {"dock": 0}

    def fake_download(url: str, target: Path) -> None:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        target.chmod(0o755)

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if command[0].endswith("vina_split"):
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="fallback")
        raise AssertionError(f"Unexpected command: {command}")

    def fake_run_command(*, command, cwd, timeout, env, phase, start_message, heartbeat_message, completion_message):  # noqa: ARG001
        calls["dock"] += 1
        output_path = Path(command[command.index("--out") + 1])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(_fake_vina_output(), encoding="utf-8")
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="Performing docking (random seed: 12345) ...\n",
            stderr="",
        )

    monkeypatch.setattr(VinaDockingPredictor, "_download_binary", staticmethod(fake_download))
    monkeypatch.setattr(VinaDockingPredictor, "_run_command_with_live_logs", staticmethod(fake_run_command))
    monkeypatch.setattr(subprocess, "run", fake_run)

    module = VinaDockingPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "runtime"),
    )
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={
            "box_center": [15.19, 53.903, 16.917],
            "box_size": [20.0, 20.0, 20.0],
        },
    )
    module.advance_window(0.0, 0.2)
    module.advance_window(0.0, 0.3)
    assert calls["dock"] == 1

    module.reset()
    _set_required_inputs(
        module,
        BioSignal,
        receptor_pdbqt_path="data/1iep/1iep_receptor.pdbqt",
        ligand_pdbqt_path="data/1iep/1iep_ligand.pdbqt",
        run_options={
            "box_center": [15.19, 53.903, 16.917],
            "box_size": [20.0, 20.0, 20.0],
        },
    )
    module.advance_window(0.0, 0.4)
    assert calls["dock"] == 2


def test_example_files_parse_and_reference_real_interface(biosim):
    repo_root = Path(__file__).resolve().parents[4]
    minimal = yaml.safe_load((repo_root / "examples" / "vina-minimal" / "config.yaml").read_text(encoding="utf-8"))
    wiring = yaml.safe_load((repo_root / "examples" / "vina-wiring" / "lab.yaml").read_text(encoding="utf-8"))

    assert minimal["model"]["path"] == "../../labs/vina-autodock-vina-docking-predictor/model"
    assert minimal["model"]["inputs"]["receptor_pdbqt_path"] == "data/1iep/1iep_receptor.pdbqt"
    assert minimal["model"]["inputs"]["ligand_pdbqt_path"] == "data/1iep/1iep_ligand.pdbqt"
    assert minimal["model"]["inputs"]["run_options"]["exhaustiveness"] == 32
    assert wiring["models"][0]["path"] == "../../labs/vina-autodock-vina-docking-predictor/model"
    assert wiring["models"][0]["parameters"]["default_receptor_pdbqt_path"] == "data/1iep/1iep_receptor.pdbqt"
    assert wiring["models"][0]["parameters"]["default_run_options"]["scoring"] == "vina"


def _schema_type(value):
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int) and not isinstance(value, bool):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    return "json"


def _signal_value(signal):
    value = signal.value
    if isinstance(value, dict) and set(value.keys()) == {"payload"}:
        return value["payload"]
    return value


def _generic_input_spec(description=None):
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


def _make_signal(*, source, name, value, emitted_at, spec=None):
    if spec is None:
        if isinstance(value, dict):
            spec = SignalSpec.record(schema={str(key): _schema_type(item) for key, item in value.items()})
        elif isinstance(value, (list, tuple)):
            spec = SignalSpec.record(schema={"payload": "json"})
        else:
            spec = SignalSpec.scalar(dtype=_schema_type(value))

    if spec.signal_type == "scalar":
        return ScalarSignal(source=source, name=name, value=value, emitted_at=emitted_at, spec=spec)
    if spec.signal_type == "array":
        return ArraySignal(source=source, name=name, value=value, emitted_at=emitted_at, spec=spec)
    if spec.signal_type == "event":
        event_value = value
        if spec.schema is not None and not (isinstance(value, dict) and set(value.keys()) == set(spec.schema.keys())):
            event_value = {"payload": value}
        return EventSignal(source=source, name=name, value=event_value, emitted_at=emitted_at, spec=spec)

    record_value = value
    if not isinstance(value, dict) or set(value.keys()) != set((spec.schema or {}).keys()):
        record_value = {"payload": value}
    return RecordSignal(source=source, name=name, value=record_value, emitted_at=emitted_at, spec=spec)
