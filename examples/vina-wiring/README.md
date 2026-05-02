# vina-wiring

Portable single-model `lab.yaml` wiring example for AutoDock Vina.

## What it runs

The same `1iep` docking job as [vina-minimal](../vina-minimal/), expressed as a schema 2.0 [lab.yaml](lab.yaml) so it can be opened in Biosimulant Desktop or run through the generic `biosim` CLI. The lab declares one model alias (`vina`) pointing at `labs/vina-autodock-vina-docking-predictor/model`, with the receptor / ligand paths and run options provided as `default_*` parameters. `wiring: []` (no inter-module wiring; this is a single-model lab).

## Requirements

- Python 3.10+ with `pyyaml` and the `biosim` package importable on the path.
- Internet access on the first run: the wrapper downloads the pinned official AutoDock Vina `1.2.7` binaries (about 1 to 2 minutes). Subsequent runs are offline.
- CPU only. No GPU required.

## Run via CLI

The generic `biosim` runner accepts both `config.yaml` and schema 2.0 `lab.yaml`:

```bash
cd /path/to/models-autodock-vina
python -m biosim examples/vina-wiring/lab.yaml --duration 0.01
```

Add `--simui --open` to launch the SimUI web dashboard in a browser instead of running headless:

```bash
python -m biosim examples/vina-wiring/lab.yaml --simui --open
```

Other flags from `python -m biosim`: `--communication-step`, `--port`, `--host`.

## Run via Desktop

1. Build a `.bsilab` package from this directory using the `biosim` packager:

   ```bash
   cd /path/to/models-autodock-vina
   python -m biosim pack build examples/vina-wiring --out vina-wiring.bsilab
   ```

2. Open `vina-wiring.bsilab` in Biosimulant Desktop. The `.bsilab` extension is registered by the desktop app's file association, so you can:
   - double-click the file in Finder / Explorer,
   - drag it onto the running app window, or
   - use the in-app "Open Lab" entry.

3. Run from the lab view. Remote runs honor the manifest-declared init overrides that force `runtime_mode: managed`, and the wrapper streams live `stdout` / `stderr` plus `BSIM_PROGRESS:` milestones from the sandbox into the desktop run log.

## Expected outputs

Same four BioSignal outputs as vina-minimal: `pose_summary`, `docking_summary`, `structure_artifacts`, `run_metadata`. The desktop renders the merged `top_rank_complex.pdb` through its built-in `structure3d` viewer.
