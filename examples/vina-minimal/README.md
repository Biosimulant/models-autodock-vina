# vina-minimal

Direct local runner config for the AutoDock Vina `1iep` receptor / ligand `PDBQT` assets.

## What it runs

A single docking job through the `vina-autodock-vina-docking-predictor` wrapper, against the checked-in `1iep` `PDBQT` pair under `labs/vina-autodock-vina-docking-predictor/model/data/1iep/`. The full input set (search box, exhaustiveness, scoring) is in [config.yaml](config.yaml).

## Requirements

- Python 3.10+ with `pyyaml` and the `biosim` package importable on the path.
- Internet access on the first run: the wrapper downloads the pinned official AutoDock Vina `1.2.7` binaries into `.runtime/vina/bin/...` (about 1 to 2 minutes). Subsequent runs are offline.
- CPU only. No GPU required.

## Run via CLI

This example is shaped for the per-repo runner [examples/run_example.py](../run_example.py).

```bash
cd /path/to/models-autodock-vina
python3 examples/run_example.py vina-minimal
```

Useful flags accepted by `run_example.py`:

- `--config <path>`: use an alternate config file (skips the positional example name).
- `--work-dir <path>`: override the module work directory (defaults to `examples/vina-minimal/runs/`).
- `--runtime-dir <path>`: override the managed runtime directory (where Vina binaries are cached).
- `--output-json <path>`: also write the final BioSignal payloads to a JSON file.

## Run via Desktop

This example is a `config.yaml` for the direct runner; it is not a desktop lab on its own. To run the same docking job in Biosimulant Desktop, use the sibling [vina-wiring](../vina-wiring/) example, which expresses the same inputs as a portable schema 2.0 `lab.yaml`. Import it once with `biosimulant labs import examples/vina-wiring`, then run it from the desktop app. See [vina-wiring/README.md](../vina-wiring/README.md) for the full flow.

## Expected outputs

The runner prints the four BioSignal outputs as JSON to stdout:

- `pose_summary`: ranked Vina poses with affinities and RMSDs.
- `docking_summary`: aggregate stats for the run.
- `structure_artifacts`: paths to per-pose `PDB` files plus the merged `top_rank_complex.pdb`.
- `run_metadata`: status, timing, and runtime info (used to compute the exit code).

Per-run artifacts (Vina logs, generated `PDBQT` outputs, merged complex) are written under `examples/vina-minimal/runs/`. The script exits `0` when `run_metadata.status == "completed"`.
