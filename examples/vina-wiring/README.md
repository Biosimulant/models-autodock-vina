# vina-wiring

Portable single-model `lab.yaml` wiring example for AutoDock Vina.

## What it runs

The same `1iep` docking job as [vina-minimal](../vina-minimal/), expressed as a schema 2.0 [lab.yaml](lab.yaml) so it can be imported into Biosimulant Desktop and run through the canonical `biosimulant` CLI. The lab declares one model alias (`vina`) pointing at `labs/vina-autodock-vina-docking-predictor/model`, with the receptor / ligand paths and run options provided as `default_*` parameters. `wiring: []` (no inter-module wiring; this is a single-model lab).

## Requirements

- The `biosimulant` CLI on your `PATH`. Install it from the desktop app (Settings > CLI Tools > Install CLI), or from the standalone installer at `https://biosimulant.com/install.sh`. Verify with `biosimulant doctor`.
- Internet access on the first run: the wrapper downloads the pinned official AutoDock Vina `1.2.7` binaries (about 1 to 2 minutes). Subsequent runs are offline.
- CPU only. No GPU required.

## Run via CLI

Import the lab into the local Biosimulant data directory, then start a run:

```bash
cd /path/to/models-autodock-vina
biosimulant labs import examples/vina-wiring
biosimulant labs list --json
```

`labs list` returns the imported lab's id. Use it from the desktop GUI, or feed it into raw mode:

```bash
biosimulant raw run_lab --input '{"lab_id":"<id-from-labs-list>"}' --json
```

The CLI is headless by default. Add `--report-file run-report.html --open` (global flags supported on `biosimulant run lab` and `runs open`) to also write a static HTML report and launch your browser when the run finishes.

## Run via Desktop

1. Import the lab once:

   ```bash
   biosimulant labs import examples/vina-wiring
   ```

2. Launch Biosimulant Desktop (or, if the app is already running, refresh the labs list). The imported lab appears under your local labs.
3. Open the lab and hit **Run**. The wrapper streams live `stdout` / `stderr` plus `BSIM_PROGRESS:` milestones from the sandbox into the desktop run log, and remote runs honor the manifest-declared init overrides that force `runtime_mode: managed`.

## Expected outputs

Same four BioSignal outputs as vina-minimal: `pose_summary`, `docking_summary`, `structure_artifacts`, `run_metadata`. The desktop renders the merged `top_rank_complex.pdb` through its built-in `structure3d` viewer.
