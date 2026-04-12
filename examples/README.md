# models-autodock-vina examples

## Included Examples

- `vina-minimal`: direct local runner config for the checked-in AutoDock Vina
  `1iep` receptor/ligand `PDBQT` assets
- `vina-wiring`: a portable single-model `space.yaml` wiring example
- `build_vina_bsispace.py`: exports a portable desktop lab and optionally stages
  it through Hub for a real remote CPU validation run

## Quick Start

Run the local example directly:

```bash
python3 examples/run_example.py vina-minimal
```

Build the portable desktop package:

```bash
python3 examples/build_vina_bsispace.py
```
