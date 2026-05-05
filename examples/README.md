# models-autodock-vina examples

## Included Examples

- [`vina-minimal`](vina-minimal/README.md): direct local runner config for the checked-in AutoDock Vina `1iep` receptor / ligand `PDBQT` assets. Run via the per-repo CLI.
- [`vina-wiring`](vina-wiring/README.md): portable single-model `lab.yaml` wiring example. Import with the `biosimulant` CLI and run from Biosimulant Desktop.

## Quick Start

### CLI (biosimulant)

Import the wiring lab into the local data directory, then list it:

```bash
biosimulant labs import examples/vina-wiring
biosimulant labs list --json
```

### Desktop

After importing with `biosimulant labs import`, open Biosimulant Desktop. The imported lab shows up under your local labs; open it and click **Run**.

See the per-example READMEs above for input details, expected outputs, and runtime requirements.
