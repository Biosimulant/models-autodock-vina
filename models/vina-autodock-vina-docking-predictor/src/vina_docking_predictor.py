# SPDX-FileCopyrightText: 2026-present Biosimulant Team
#
# SPDX-License-Identifier: MIT
"""AutoDock Vina single-complex docking wrapper."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.request
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

from biosim import BioModule
from biosim.signals import BioSignal


_ALLOWED_RUN_OPTIONS = {
    "box_center",
    "box_size",
    "exhaustiveness",
    "n_poses",
    "energy_range",
    "cpu",
    "seed",
    "scoring",
}
_PATH_LIKE_SUFFIXES = {".pdbqt"}
_VINA_RESULT_RE = re.compile(
    r"^REMARK VINA RESULT:\s*(?P<affinity>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<rmsd_lb>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<rmsd_ub>[-+]?\d+(?:\.\d+)?)\s*$"
)
_SEED_RE = re.compile(r"random seed:\s*(?P<seed>[-+]?\d+)", re.IGNORECASE)
_ELEMENT_OVERRIDES = {
    "A": "C",
    "C": "C",
    "N": "N",
    "NA": "N",
    "OA": "O",
    "O": "O",
    "S": "S",
    "SA": "S",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "SI": "Si",
    "HD": "H",
    "H": "H",
    "MG": "Mg",
    "ZN": "Zn",
    "CA": "Ca",
    "FE": "Fe",
}


@dataclass(frozen=True)
class _RuntimeBinaries:
    vina: str
    vina_split: str
    platform_tag: str


def _coerce_string(value: Any, *preferred_keys: str) -> Optional[str]:
    if isinstance(value, str):
        text = value.strip()
        return text or None
    if isinstance(value, Mapping):
        for key in preferred_keys:
            candidate = value.get(key)
            if isinstance(candidate, str):
                text = candidate.strip()
                if text:
                    return text
    return None


def _coerce_run_options(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        return {}
    out: dict[str, Any] = {}
    for key, item in value.items():
        if isinstance(key, str):
            out[key] = item
    return out


def _tail_text(value: str, *, max_chars: int = 12_000) -> str:
    return value[-max_chars:] if len(value) > max_chars else value


class VinaDockingPredictor(BioModule):
    """Run classic AutoDock Vina for a single receptor and ligand PDBQT pair."""

    def __init__(
        self,
        default_receptor_pdbqt_path: Optional[str] = None,
        default_ligand_pdbqt_path: Optional[str] = None,
        default_run_options: Optional[Mapping[str, Any]] = None,
        runtime_mode: str = "managed",
        runtime_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        work_dir: Optional[str] = None,
        vina_version: str = "1.2.7",
        vina_release_base_url: str = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download",
        vina_executable: Optional[str] = None,
        vina_split_executable: Optional[str] = None,
        command_timeout_s: float = 3_600.0,
        download_timeout_s: float = 1_200.0,
        progress_heartbeat_s: float = 15.0,
        min_dt: float = 0.01,
    ) -> None:
        self.min_dt = min_dt
        self.runtime_mode = runtime_mode
        self.vina_version = vina_version
        self.vina_release_base_url = vina_release_base_url.rstrip("/")
        self.vina_executable = vina_executable
        self.vina_split_executable = vina_split_executable
        self.command_timeout_s = float(command_timeout_s)
        self.download_timeout_s = float(download_timeout_s)
        self.progress_heartbeat_s = max(0.0, float(progress_heartbeat_s))
        self.work_dir = Path(work_dir).expanduser().resolve() if work_dir else None

        self.model_root = Path(__file__).resolve().parents[1]
        repo_root = self.model_root.parents[1]
        self.runtime_dir = (
            Path(runtime_dir).expanduser().resolve()
            if runtime_dir
            else (repo_root / ".runtime" / "vina").resolve()
        )
        self.cache_dir = (
            Path(cache_dir).expanduser().resolve()
            if cache_dir
            else (self.runtime_dir / "downloads").resolve()
        )

        self._receptor_pdbqt_path: Optional[str] = _coerce_string(
            default_receptor_pdbqt_path, "path"
        )
        self._ligand_pdbqt_path: Optional[str] = _coerce_string(
            default_ligand_pdbqt_path, "path"
        )
        self._run_options: dict[str, Any] = _coerce_run_options(default_run_options)
        self._outputs: dict[str, BioSignal] = {}
        self._cached_payloads: dict[str, Any] = {}
        self._last_signature: Optional[str] = None

    def inputs(self) -> set[str]:
        return {"receptor_pdbqt_path", "ligand_pdbqt_path", "run_options"}

    def outputs(self) -> set[str]:
        return {"pose_summary", "docking_summary", "structure_artifacts", "run_metadata"}

    def reset(self) -> None:
        self._outputs = {}
        self._cached_payloads = {}
        self._last_signature = None

    def set_inputs(self, signals: dict[str, BioSignal]) -> None:
        changed = False

        receptor_signal = signals.get("receptor_pdbqt_path")
        if receptor_signal is not None:
            receptor_pdbqt_path = _coerce_string(receptor_signal.value, "path")
            if receptor_pdbqt_path != self._receptor_pdbqt_path:
                self._receptor_pdbqt_path = receptor_pdbqt_path
                changed = True

        ligand_signal = signals.get("ligand_pdbqt_path")
        if ligand_signal is not None:
            ligand_pdbqt_path = _coerce_string(ligand_signal.value, "path")
            if ligand_pdbqt_path != self._ligand_pdbqt_path:
                self._ligand_pdbqt_path = ligand_pdbqt_path
                changed = True

        run_signal = signals.get("run_options")
        if run_signal is not None:
            run_options = _coerce_run_options(run_signal.value)
            if run_options != self._run_options:
                self._run_options = run_options
                changed = True

        if changed:
            self._last_signature = None

    def advance_to(self, t: float) -> None:
        metadata: dict[str, Any] = {
            "status": "running",
            "runtime_mode": self.runtime_mode,
            "runtime_dir": str(self.runtime_dir),
            "cache_dir": str(self.cache_dir),
            "vina_version": self.vina_version,
            "runtime_bootstrapped": False,
            "stdout": "",
            "stderr": "",
        }

        self._emit_progress("inputs", "Validating AutoDock Vina inputs")
        try:
            resolved_options = self._resolved_options()
            receptor_path = self._resolve_required_input_path(
                self._receptor_pdbqt_path,
                input_name="receptor_pdbqt_path",
            )
            ligand_path = self._resolve_required_input_path(
                self._ligand_pdbqt_path,
                input_name="ligand_pdbqt_path",
            )
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = str(exc)
            self._emit_progress("error", metadata["error"])
            self._set_error_payload(str(exc), metadata=metadata)
            self._emit_outputs(t)
            return

        signature = json.dumps(
            {
                "receptor_pdbqt_path": receptor_path,
                "ligand_pdbqt_path": ligand_path,
                "run_options": resolved_options,
            },
            sort_keys=True,
        )
        if signature == self._last_signature and self._cached_payloads:
            self._emit_progress("cache", "Reusing cached AutoDock Vina outputs for unchanged inputs")
            self._emit_outputs(t)
            return

        run_root = self._create_run_root()
        output_dir = run_root / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        ligand_stem = Path(ligand_path).stem
        docking_output_path = output_dir / f"{ligand_stem}_vina_out.pdbqt"
        config_path = output_dir / "vina.config.txt"
        stdout_tail_file = output_dir / "stdout_tail.log"
        stderr_tail_file = output_dir / "stderr_tail.log"

        metadata["run_root"] = str(run_root)
        metadata["output_dir"] = str(output_dir)
        metadata["receptor_pdbqt_path"] = receptor_path
        metadata["ligand_pdbqt_path"] = ligand_path
        metadata["docking_output_path"] = str(docking_output_path)
        metadata["config_file"] = str(config_path)

        try:
            self._emit_progress("runtime", "Preparing AutoDock Vina runtime")
            runtime = self._prepare_runtime(metadata)
            metadata["resolved_vina_executable"] = runtime.vina
            metadata["resolved_vina_split_executable"] = runtime.vina_split
            metadata["platform_tag"] = runtime.platform_tag

            self._emit_progress("config", "Writing AutoDock Vina config")
            self._write_config(config_path, resolved_options)
            command = self._build_command(
                vina_executable=runtime.vina,
                receptor_path=receptor_path,
                ligand_path=ligand_path,
                config_path=config_path,
                output_path=docking_output_path,
                options=resolved_options,
            )
            metadata["command"] = command

            completed = self._run_command_with_live_logs(
                command=command,
                cwd=output_dir,
                timeout=self.command_timeout_s,
                env=self._command_env(),
                phase="docking",
                start_message="Starting AutoDock Vina docking",
                heartbeat_message="AutoDock Vina docking is still running",
                completion_message="AutoDock Vina docking finished",
            )
            metadata["returncode"] = completed.returncode
            metadata["stdout"] = _tail_text(completed.stdout or "")
            metadata["stderr"] = _tail_text(completed.stderr or "")
            stdout_tail_file.write_text(metadata["stdout"], encoding="utf-8")
            stderr_tail_file.write_text(metadata["stderr"], encoding="utf-8")
            parsed_seed = self._parse_seed(completed.stdout or "")
            if parsed_seed is not None:
                metadata["seed"] = parsed_seed
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = f"failed to execute AutoDock Vina: {exc}"
            self._emit_progress("error", metadata["error"])
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        if completed.returncode != 0:
            metadata["status"] = "error"
            metadata["error"] = "AutoDock Vina returned a non-zero exit code"
            self._emit_progress("error", metadata["error"])
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        try:
            self._emit_progress("postprocess", "Collecting ranked docking poses")
            pose_records = self._parse_pose_records(docking_output_path)
            pose_files = self._split_pose_files(
                vina_split_executable=runtime.vina_split,
                output_path=docking_output_path,
                output_dir=output_dir,
            )
            pose_records = self._attach_pose_files(pose_records, pose_files)

            pose_summary_file = output_dir / "pose_summary.json"
            docking_summary = self._build_docking_summary(
                pose_records=pose_records,
                options=resolved_options,
                seed=metadata.get("seed"),
            )
            docking_summary_file = output_dir / "docking_summary.json"
            self._write_json(pose_summary_file, {"poses": pose_records})
            self._write_json(docking_summary_file, docking_summary)

            top_pose_path = Path(str(pose_records[0]["pose_file"])).resolve()
            top_complex_file = output_dir / "top_rank_complex.pdb"
            self._emit_progress("postprocess", "Building merged receptor-ligand complex")
            self._build_top_rank_complex(
                receptor_pdbqt_path=Path(receptor_path),
                ligand_pose_path=top_pose_path,
                output_path=top_complex_file,
            )

            self._emit_progress("outputs", "Publishing AutoDock Vina artifacts and visuals")
            artifacts = self._build_structure_artifacts(
                output_dir=output_dir,
                pose_records=pose_records,
                docking_output_path=docking_output_path,
                config_path=config_path,
                top_complex_file=top_complex_file,
                pose_summary_file=pose_summary_file,
                docking_summary_file=docking_summary_file,
                stdout_tail_file=stdout_tail_file,
                stderr_tail_file=stderr_tail_file,
            )
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = f"expected AutoDock Vina outputs were not found: {exc}"
            self._emit_progress("error", metadata["error"])
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        metadata["status"] = "completed"
        metadata["seed"] = metadata.get("seed", resolved_options.get("seed"))
        self._cached_payloads = {
            "pose_summary": pose_records,
            "docking_summary": docking_summary,
            "structure_artifacts": artifacts,
            "run_metadata": metadata,
        }
        self._last_signature = signature
        self._emit_progress("completed", "AutoDock Vina outputs are ready")
        self._emit_outputs(t)

    def get_outputs(self) -> dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[list[dict[str, Any]]]:
        run_metadata = self._cached_payloads.get("run_metadata")
        artifacts = self._cached_payloads.get("structure_artifacts")
        docking_summary = self._cached_payloads.get("docking_summary")
        poses = self._cached_payloads.get("pose_summary")

        if not isinstance(run_metadata, Mapping) or run_metadata.get("status") != "completed":
            return None
        if not isinstance(artifacts, Mapping) or not isinstance(docking_summary, Mapping):
            return None
        if not isinstance(poses, list):
            return None

        top_complex_raw = artifacts.get("top_complex_file")
        if not isinstance(top_complex_raw, str) or not top_complex_raw:
            return None

        top_complex_path = Path(top_complex_raw).expanduser().resolve()
        rows: list[list[str]] = []
        for row in poses:
            if not isinstance(row, Mapping):
                continue
            rows.append(
                [
                    str(row.get("rank", "")),
                    "" if row.get("affinity_kcal_mol") is None else str(row.get("affinity_kcal_mol")),
                    "" if row.get("rmsd_lb") is None else str(row.get("rmsd_lb")),
                    "" if row.get("rmsd_ub") is None else str(row.get("rmsd_ub")),
                    Path(str(row.get("pose_file", ""))).name,
                ]
            )

        return [
            {
                "render": "structure3d",
                "description": "Top-ranked AutoDock Vina complex for the latest docking run.",
                "data": {
                    "title": "Top-Ranked Docked Complex",
                    "source": {
                        "kind": "artifact",
                        "artifact_id": self._structure_artifact_id(top_complex_path),
                        "path": str(top_complex_path),
                    },
                    "format": "pdb",
                    "annotations": [
                        {
                            "label": "Top Pose Affinity (kcal/mol)",
                            "value": docking_summary.get("top_pose_affinity_kcal_mol"),
                        },
                        {
                            "label": "Scoring",
                            "value": docking_summary.get("scoring"),
                        },
                        {
                            "label": "Pose Count",
                            "value": docking_summary.get("pose_count"),
                        },
                    ],
                    "initial_view": {"reset_camera": True},
                },
            },
            {
                "render": "table",
                "description": "Ranked pose summary from the latest AutoDock Vina run.",
                "data": {
                    "title": "AutoDock Vina Pose Summary",
                    "columns": ["Rank", "Affinity", "RMSD l.b.", "RMSD u.b.", "Pose File"],
                    "rows": rows,
                },
            },
        ]

    def _resolved_options(self) -> dict[str, Any]:
        resolved: dict[str, Any] = {
            "exhaustiveness": 8,
            "n_poses": 9,
            "energy_range": 3.0,
            "cpu": 0,
            "scoring": "vina",
        }

        for key in self._run_options:
            if key not in _ALLOWED_RUN_OPTIONS:
                raise ValueError(f"unsupported run_options key: {key}")

        resolved["box_center"] = self._coerce_triplet(
            self._run_options.get("box_center"),
            option_name="box_center",
        )
        resolved["box_size"] = self._coerce_triplet(
            self._run_options.get("box_size"),
            option_name="box_size",
        )

        for key in ("exhaustiveness", "n_poses"):
            if key in self._run_options:
                value = self._run_options.get(key)
                if not isinstance(value, int) or value <= 0:
                    raise ValueError(f"run_options.{key} must be a positive integer")
                resolved[key] = value

        if "energy_range" in self._run_options:
            energy_range = self._run_options.get("energy_range")
            if not isinstance(energy_range, (int, float)) or float(energy_range) <= 0:
                raise ValueError("run_options.energy_range must be a positive number")
            resolved["energy_range"] = float(energy_range)

        if "cpu" in self._run_options:
            cpu = self._run_options.get("cpu")
            if not isinstance(cpu, int) or cpu < 0:
                raise ValueError("run_options.cpu must be a non-negative integer")
            resolved["cpu"] = cpu

        if "seed" in self._run_options:
            seed = self._run_options.get("seed")
            if not isinstance(seed, int):
                raise ValueError("run_options.seed must be an integer")
            resolved["seed"] = seed

        if "scoring" in self._run_options:
            scoring = _coerce_string(self._run_options.get("scoring"))
            if scoring not in {"vina", "vinardo"}:
                raise ValueError("run_options.scoring must be 'vina' or 'vinardo'")
            resolved["scoring"] = scoring

        return resolved

    def _coerce_triplet(self, value: Any, *, option_name: str) -> list[float]:
        if not isinstance(value, (list, tuple)) or len(value) != 3:
            raise ValueError(f"run_options.{option_name} must be a list of three numbers")
        out: list[float] = []
        for item in value:
            if not isinstance(item, (int, float)):
                raise ValueError(f"run_options.{option_name} must contain only numbers")
            out.append(float(item))
        return out

    def _create_run_root(self) -> Path:
        base_dir = self.work_dir
        if base_dir is not None:
            base_dir.mkdir(parents=True, exist_ok=True)
        root = Path(
            tempfile.mkdtemp(
                prefix="vina-run-",
                dir=str(base_dir) if base_dir is not None else None,
            )
        )
        return root.resolve()

    def _prepare_runtime(self, metadata: dict[str, Any]) -> _RuntimeBinaries:
        mode = self.runtime_mode.strip().lower()
        if mode == "managed":
            return self._ensure_managed_runtime(metadata)
        if mode == "external":
            binaries = self._resolve_external_binaries()
            self._emit_progress("runtime", f"Using external AutoDock Vina binaries at {binaries.vina}")
            return binaries
        raise ValueError(f"unsupported runtime_mode: {self.runtime_mode}")

    def _ensure_managed_runtime(self, metadata: dict[str, Any]) -> _RuntimeBinaries:
        self.runtime_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        platform_tag = self._platform_tag()
        bin_dir = (self.runtime_dir / "bin" / self.vina_version / platform_tag).resolve()
        bin_dir.mkdir(parents=True, exist_ok=True)
        vina_target = bin_dir / "vina"
        vina_split_target = bin_dir / "vina_split"

        downloaded = False
        for executable_name, target in (("vina", vina_target), ("vina_split", vina_split_target)):
            if target.is_file() and os.access(target, os.X_OK):
                continue
            downloaded = True
            asset_name = f"{executable_name}_{self.vina_version}_{platform_tag}"
            url = f"{self.vina_release_base_url}/v{self.vina_version}/{asset_name}"
            self._emit_progress("runtime", f"Downloading {asset_name}")
            self._download_binary(url, target)

        if downloaded:
            metadata["runtime_bootstrapped"] = True
            self._emit_progress("runtime", "Managed AutoDock Vina runtime is ready")
        else:
            self._emit_progress("runtime", "Reusing cached managed AutoDock Vina runtime")

        return _RuntimeBinaries(
            vina=str(vina_target),
            vina_split=str(vina_split_target),
            platform_tag=platform_tag,
        )

    def _resolve_external_binaries(self) -> _RuntimeBinaries:
        vina = self._resolve_external_executable(self.vina_executable, "vina")
        vina_split = self._resolve_external_executable(self.vina_split_executable, "vina_split")
        return _RuntimeBinaries(vina=vina, vina_split=vina_split, platform_tag="external")

    def _resolve_external_executable(self, configured: Optional[str], fallback_name: str) -> str:
        if isinstance(configured, str) and configured.strip():
            candidate = Path(configured).expanduser()
            if candidate.is_file():
                return str(candidate.resolve())
        resolved = shutil.which(fallback_name)
        if resolved:
            return resolved
        raise FileNotFoundError(
            f"could not resolve external executable `{fallback_name}`; provide an explicit path"
        )

    def _download_binary(self, url: str, target: Path) -> None:
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp_target = (self.cache_dir / f"{target.name}.{hashlib.sha256(url.encode('utf-8')).hexdigest()[:12]}.part").resolve()
        if tmp_target.exists():
            tmp_target.unlink()
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "biosimulant-vina-wrapper/0.1"},
            )
            with urllib.request.urlopen(request, timeout=self.download_timeout_s) as response:
                with tmp_target.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
            tmp_target.chmod(0o755)
            shutil.move(str(tmp_target), str(target))
            target.chmod(0o755)
        finally:
            if tmp_target.exists():
                tmp_target.unlink()

    def _platform_tag(self) -> str:
        system = platform.system().lower()
        machine = platform.machine().lower()
        if machine in {"x86_64", "amd64"}:
            arch = "x86_64"
        elif machine in {"arm64", "aarch64"}:
            arch = "aarch64"
        else:
            raise RuntimeError(f"unsupported architecture for AutoDock Vina managed runtime: {machine}")
        if system == "linux":
            return f"linux_{arch}"
        if system == "darwin":
            return f"mac_{arch}"
        raise RuntimeError(f"unsupported platform for AutoDock Vina managed runtime: {system}")

    def _write_config(self, path: Path, options: dict[str, Any]) -> None:
        center_x, center_y, center_z = options["box_center"]
        size_x, size_y, size_z = options["box_size"]
        contents = (
            f"center_x = {center_x:.3f}\n"
            f"center_y = {center_y:.3f}\n"
            f"center_z = {center_z:.3f}\n"
            f"size_x = {size_x:.3f}\n"
            f"size_y = {size_y:.3f}\n"
            f"size_z = {size_z:.3f}\n"
        )
        path.write_text(contents, encoding="utf-8")

    def _build_command(
        self,
        *,
        vina_executable: str,
        receptor_path: str,
        ligand_path: str,
        config_path: Path,
        output_path: Path,
        options: dict[str, Any],
    ) -> list[str]:
        command = [
            vina_executable,
            "--receptor",
            receptor_path,
            "--ligand",
            ligand_path,
            "--config",
            str(config_path.resolve()),
            "--scoring",
            str(options["scoring"]),
            "--exhaustiveness",
            str(options["exhaustiveness"]),
            "--num_modes",
            str(options["n_poses"]),
            "--energy_range",
            str(options["energy_range"]),
            "--cpu",
            str(options["cpu"]),
            "--out",
            str(output_path.resolve()),
        ]
        if "seed" in options:
            command.extend(["--seed", str(options["seed"])])
        return command

    def _command_env(self) -> dict[str, str]:
        env = dict(os.environ)
        env.setdefault("PYTHONUNBUFFERED", "1")
        env.setdefault("XDG_CACHE_HOME", str((self.cache_dir / "xdg").resolve()))
        return env

    def _run_command_with_live_logs(
        self,
        *,
        command: list[str],
        cwd: Path,
        timeout: float,
        env: Optional[dict[str, str]],
        phase: str,
        start_message: str,
        heartbeat_message: str,
        completion_message: str,
    ) -> subprocess.CompletedProcess[str]:
        self._emit_progress(phase, start_message)

        started_at = time.monotonic()
        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=env,
        )
        stdout_lines: list[str] = []
        stderr_lines: list[str] = []

        def _pump(stream, buffer: list[str], writer) -> None:
            if stream is None:
                return
            try:
                for line in iter(stream.readline, ""):
                    buffer.append(line)
                    writer.write(line)
                    writer.flush()
            finally:
                stream.close()

        stdout_thread = threading.Thread(
            target=_pump,
            args=(process.stdout, stdout_lines, sys.stdout),
            daemon=True,
        )
        stderr_thread = threading.Thread(
            target=_pump,
            args=(process.stderr, stderr_lines, sys.stderr),
            daemon=True,
        )
        stdout_thread.start()
        stderr_thread.start()

        tick = 0
        try:
            while True:
                returncode = process.poll()
                if returncode is not None:
                    break
                elapsed = time.monotonic() - started_at
                if elapsed > timeout:
                    process.kill()
                    raise TimeoutError(f"timed out after {int(timeout)}s: {' '.join(command)}")
                if self.progress_heartbeat_s > 0:
                    current_tick = int(elapsed / self.progress_heartbeat_s)
                    if current_tick > tick:
                        tick = current_tick
                        self._emit_progress(
                            phase,
                            f"{heartbeat_message} ({int(elapsed)}s elapsed)",
                            tick=tick,
                            duration=float(int(elapsed)),
                        )
                time.sleep(0.1)
        finally:
            stdout_thread.join(timeout=5.0)
            stderr_thread.join(timeout=5.0)

        stdout = "".join(stdout_lines)
        stderr = "".join(stderr_lines)
        completed = subprocess.CompletedProcess(
            command,
            process.returncode or 0,
            stdout=stdout,
            stderr=stderr,
        )
        if completed.returncode == 0:
            self._emit_progress(phase, completion_message)
        else:
            self._emit_progress(
                phase,
                f"{completion_message} (exit code {completed.returncode})",
            )
        return completed

    def _parse_seed(self, stdout: str) -> Optional[int]:
        match = _SEED_RE.search(stdout)
        if match is None:
            return None
        try:
            return int(match.group("seed"))
        except ValueError:
            return None

    def _parse_pose_records(self, output_path: Path) -> list[dict[str, Any]]:
        pose_records: list[dict[str, Any]] = []
        current_rank: Optional[int] = None
        current_result: tuple[float, float, float] | None = None

        for raw in output_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if raw.startswith("MODEL"):
                parts = raw.split()
                try:
                    current_rank = int(parts[1]) if len(parts) > 1 else len(pose_records) + 1
                except ValueError:
                    current_rank = len(pose_records) + 1
                current_result = None
                continue
            match = _VINA_RESULT_RE.match(raw.strip())
            if match is not None:
                current_result = (
                    float(match.group("affinity")),
                    float(match.group("rmsd_lb")),
                    float(match.group("rmsd_ub")),
                )
                continue
            if raw.startswith("ENDMDL"):
                if current_result is None:
                    raise RuntimeError(f"missing Vina result annotation in {output_path}")
                pose_records.append(
                    {
                        "rank": current_rank or len(pose_records) + 1,
                        "affinity_kcal_mol": current_result[0],
                        "rmsd_lb": current_result[1],
                        "rmsd_ub": current_result[2],
                    }
                )
                current_rank = None
                current_result = None

        if not pose_records:
            raise FileNotFoundError(f"no AutoDock Vina pose records found under {output_path}")
        pose_records.sort(key=lambda item: int(item["rank"]))
        return pose_records

    def _split_pose_files(
        self,
        *,
        vina_split_executable: str,
        output_path: Path,
        output_dir: Path,
    ) -> list[Path]:
        split_dir = output_dir / "split_poses"
        split_dir.mkdir(parents=True, exist_ok=True)
        before = {path.name for path in split_dir.glob("*.pdbqt")}
        completed = subprocess.run(
            [vina_split_executable, "--input", str(output_path.resolve())],
            cwd=str(split_dir),
            capture_output=True,
            text=True,
            timeout=min(self.command_timeout_s, 600.0),
            check=False,
        )
        after = sorted(path for path in split_dir.glob("*.pdbqt") if path.name not in before)
        if completed.returncode == 0 and after:
            return self._sort_pose_files(after)
        return self._manual_split_pose_files(output_path, split_dir)

    def _manual_split_pose_files(self, output_path: Path, split_dir: Path) -> list[Path]:
        pose_files: list[Path] = []
        current_rank: Optional[int] = None
        current_lines: list[str] = []
        for raw in output_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if raw.startswith("MODEL"):
                if current_lines and current_rank is not None:
                    pose_files.append(self._write_pose_file(split_dir, current_rank, current_lines))
                parts = raw.split()
                try:
                    current_rank = int(parts[1]) if len(parts) > 1 else len(pose_files) + 1
                except ValueError:
                    current_rank = len(pose_files) + 1
                current_lines = [raw]
                continue
            if current_rank is None:
                continue
            current_lines.append(raw)
            if raw.startswith("ENDMDL"):
                pose_files.append(self._write_pose_file(split_dir, current_rank, current_lines))
                current_rank = None
                current_lines = []

        if current_lines and current_rank is not None:
            pose_files.append(self._write_pose_file(split_dir, current_rank, current_lines))
        if not pose_files:
            single_pose = split_dir / "rank_1.pdbqt"
            single_pose.write_text(output_path.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
            pose_files.append(single_pose)
        return self._sort_pose_files(pose_files)

    def _write_pose_file(self, split_dir: Path, rank: int, lines: list[str]) -> Path:
        pose_path = split_dir / f"rank_{rank}.pdbqt"
        pose_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return pose_path.resolve()

    def _sort_pose_files(self, pose_files: list[Path]) -> list[Path]:
        def _key(path: Path) -> tuple[int, str]:
            match = re.search(r"(\d+)", path.stem)
            return (int(match.group(1)) if match else 10_000, path.name)

        return sorted((path.resolve() for path in pose_files), key=_key)

    def _attach_pose_files(
        self,
        pose_records: list[dict[str, Any]],
        pose_files: list[Path],
    ) -> list[dict[str, Any]]:
        if not pose_files:
            raise RuntimeError("no split AutoDock Vina pose files were produced")
        out: list[dict[str, Any]] = []
        for index, record in enumerate(pose_records):
            if index >= len(pose_files):
                raise RuntimeError("split AutoDock Vina pose count did not match parsed score count")
            next_record = dict(record)
            next_record["pose_file"] = str(pose_files[index])
            out.append(next_record)
        return out

    def _build_docking_summary(
        self,
        *,
        pose_records: list[dict[str, Any]],
        options: dict[str, Any],
        seed: Any,
    ) -> dict[str, Any]:
        top = pose_records[0]
        return {
            "top_pose_rank": top["rank"],
            "top_pose_affinity_kcal_mol": top.get("affinity_kcal_mol"),
            "pose_count": len(pose_records),
            "scoring": options["scoring"],
            "exhaustiveness": options["exhaustiveness"],
            "n_poses": options["n_poses"],
            "energy_range": options["energy_range"],
            "cpu": options["cpu"],
            "seed": seed,
            "box_center": list(options["box_center"]),
            "box_size": list(options["box_size"]),
        }

    def _build_structure_artifacts(
        self,
        *,
        output_dir: Path,
        pose_records: list[dict[str, Any]],
        docking_output_path: Path,
        config_path: Path,
        top_complex_file: Path,
        pose_summary_file: Path,
        docking_summary_file: Path,
        stdout_tail_file: Path,
        stderr_tail_file: Path,
    ) -> dict[str, Any]:
        artifacts: dict[str, Any] = {
            "output_dir": str(output_dir.resolve()),
            "docking_output_file": str(docking_output_path.resolve()),
            "config_file": str(config_path.resolve()),
            "top_complex_file": str(top_complex_file.resolve()),
            "pose_summary_file": str(pose_summary_file.resolve()),
            "docking_summary_file": str(docking_summary_file.resolve()),
            "stdout_tail_file": str(stdout_tail_file.resolve()),
            "stderr_tail_file": str(stderr_tail_file.resolve()),
        }
        for pose in pose_records:
            rank = int(pose["rank"])
            artifacts[f"rank_{rank}_file"] = str(Path(str(pose["pose_file"])).resolve())
            if rank == 1:
                artifacts["top_pose_file"] = artifacts[f"rank_{rank}_file"]
        return artifacts

    def _build_top_rank_complex(
        self,
        *,
        receptor_pdbqt_path: Path,
        ligand_pose_path: Path,
        output_path: Path,
    ) -> None:
        receptor_lines = self._pdb_lines_from_pdbqt(
            receptor_pdbqt_path,
            serial_start=1,
            ligand=False,
        )
        ligand_lines = self._pdb_lines_from_pdbqt(
            ligand_pose_path,
            serial_start=len(receptor_lines) + 1,
            ligand=True,
        )
        output_path.write_text(
            "\n".join(receptor_lines + ["TER"] + ligand_lines + ["END"]) + "\n",
            encoding="utf-8",
        )

    def _pdb_lines_from_pdbqt(
        self,
        path: Path,
        *,
        serial_start: int,
        ligand: bool,
    ) -> list[str]:
        lines: list[str] = []
        serial = serial_start
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not raw.startswith(("ATOM", "HETATM")):
                continue
            atom_name = (raw[12:16].strip() if len(raw) >= 16 else "") or "C"
            record = "HETATM" if ligand else (raw[0:6].strip() or "ATOM")
            res_name = "LIG" if ligand else ((raw[17:20].strip() if len(raw) >= 20 else "") or "REC")
            chain = "Z" if ligand else ((raw[21].strip() if len(raw) >= 22 else "") or "A")
            res_seq = 1 if ligand else self._safe_int(raw[22:26], default=1)
            x = self._safe_float(raw[30:38], default=0.0)
            y = self._safe_float(raw[38:46], default=0.0)
            z = self._safe_float(raw[46:54], default=0.0)
            occupancy = self._safe_float(raw[54:60], default=1.0)
            temp_factor = self._safe_float(raw[60:66], default=0.0)
            autodock_type = raw[77:].strip().split()[-1] if len(raw) > 77 and raw[77:].strip() else ""
            element = self._infer_element(atom_name, autodock_type)
            lines.append(
                f"{record:<6}{serial:5d} {atom_name[:4]:>4} {res_name[:3]:>3} {chain[:1]}{res_seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}          {element:>2}"
            )
            serial += 1
        if not lines:
            raise RuntimeError(f"no atoms were extracted from {path}")
        return lines

    def _infer_element(self, atom_name: str, autodock_type: str) -> str:
        if autodock_type:
            normalized_type = autodock_type.strip().upper()
            if normalized_type in _ELEMENT_OVERRIDES:
                return _ELEMENT_OVERRIDES[normalized_type]
        trimmed = "".join(ch for ch in atom_name if ch.isalpha())
        if not trimmed:
            return "C"
        if len(trimmed) >= 2 and trimmed[:2].upper() in {"CL", "BR", "SI", "ZN", "MG", "FE", "CA"}:
            return trimmed[:2].title()
        return trimmed[0].upper()

    def _safe_float(self, raw: str, *, default: float) -> float:
        try:
            return float(raw.strip())
        except Exception:  # noqa: BLE001
            return default

    def _safe_int(self, raw: str, *, default: int) -> int:
        try:
            return int(raw.strip())
        except Exception:  # noqa: BLE001
            return default

    def _resolve_required_input_path(self, raw: Optional[str], *, input_name: str) -> str:
        if raw is None:
            raise ValueError(f"{input_name} input is required")
        for candidate in self._path_candidates(raw):
            if candidate.exists() and candidate.is_file():
                return str(candidate.resolve())
        raise FileNotFoundError(f"{input_name} file does not exist: {raw}")

    def _path_candidates(self, raw: str) -> list[Path]:
        path = Path(raw).expanduser()
        if path.is_absolute():
            return [path]
        return [(self.model_root / path).resolve(), (Path.cwd() / path).resolve()]

    def _looks_like_path(self, value: str) -> bool:
        if value.startswith(".") or "/" in value or "\\" in value:
            return True
        return Path(value).suffix.lower() in _PATH_LIKE_SUFFIXES

    def _write_json(self, path: Path, payload: Any) -> None:
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")

    def _emit_progress(self, phase: str, message: str, **payload: Any) -> None:
        event: dict[str, Any] = {"phase": phase, "message": message}
        for key, value in payload.items():
            if value is not None:
                event[key] = value
        print(f"BSIM_PROGRESS:{json.dumps(event, sort_keys=True)}", flush=True)

    def _set_error_payload(self, error_message: str, *, metadata: Optional[dict[str, Any]] = None) -> None:
        next_metadata = dict(metadata or {})
        next_metadata.setdefault("status", "error")
        next_metadata.setdefault("error", error_message)
        self._cached_payloads = {
            "pose_summary": [],
            "docking_summary": {},
            "structure_artifacts": {},
            "run_metadata": next_metadata,
        }

    def _emit_outputs(self, t: float) -> None:
        self._outputs = {}
        for name in self.outputs():
            self._outputs[name] = BioSignal(
                source="vina",
                name=name,
                value=self._cached_payloads.get(name, {}),
                time=t,
            )

    def _structure_artifact_id(self, path: Path) -> str:
        digest = hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:16]
        return f"structure-{digest}"
