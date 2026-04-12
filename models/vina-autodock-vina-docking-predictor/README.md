# AutoDock Vina Docking Predictor

`vina-autodock-vina-docking-predictor` is a native `biosim.BioModule` wrapper
around the official AutoDock Vina `1.2.7` CLI for **single-complex** docking.

## Public Interface

Inputs:
- `receptor_pdbqt_path`: path to a prepared receptor `.pdbqt`
- `ligand_pdbqt_path`: path to a prepared ligand `.pdbqt`
- `run_options`: map with required `box_center` and `box_size`, plus optional
  `exhaustiveness`, `n_poses`, `energy_range`, `cpu`, `seed`, and `scoring`

Outputs:
- `pose_summary`
- `docking_summary`
- `structure_artifacts`
- `run_metadata`

## Runtime Model

Local default:
- `runtime_mode: managed`
- downloads the pinned official AutoDock Vina `1.2.7` release binaries into
  `.runtime/vina/bin/...`

Remote default:
- manifest-declared remote init overrides force `runtime_mode: managed`
- the wrapper downloads the same pinned binaries under the mounted remote
  runtime/cache roots and executes CPU-backed docking there
- the wrapper streams live child stdout and stderr while also emitting
  `BSIM_PROGRESS:` milestones for the desktop/web feedback loop

The first real run needs:
- internet access to download the official release binaries
- prepared receptor and ligand `PDBQT` files

## Example Assets

The checked-in assets under `data/1iep/` use the official upstream AutoDock
Vina basic-docking `1iep` example:

- `1iep_receptor.pdbqt`
- `1iep_ligand.pdbqt`
- `1iep_receptor.box.txt`
- `1iep_receptor.box.pdb`
- `1iep_receptorH.pdb`
- `1iep_ligand.sdf`
