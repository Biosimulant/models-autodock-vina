# AutoDock Vina: inspect a single-complex docking result

Dock one prepared ligand into one prepared rigid receptor using AutoDock Vina 1.2.7 on CPU. The Lab returns ranked PDBQT poses, a merged receptor–ligand PDB for the 3D viewer, resolved search settings, random seed and input-file SHA-256 checksums. It does not prepare raw structures or run virtual-screening batches.

## Start with the bundled example

The packaged `1iep` inputs use a search-box center of `[15.19, 53.903, 16.917]` Å and size `[20, 20, 20]` Å. Defaults are exhaustiveness 32, up to 9 poses, energy range 3 kcal/mol, two CPU threads and seed 2026. The fixed seed is a reproducibility setting, not a tuned biological parameter. The default example is documented in the [official Vina tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html).

The existing screenshots in `assets/` show an earlier release; they are not evidence for version 1.1.0. Current results must come from an actual run.

## Supply a different complex

- `receptor_pdbqt_path`: a prepared receptor PDBQT file accessible in the execution environment.
- `ligand_pdbqt_path`: a prepared ligand PDBQT file accessible in the execution environment.
- `run_options`: an object with partial overrides of the Lab defaults. Specify a suitable `box_center` and `box_size` for the new receptor; the bundled box is specific to the example.

Supported options are `box_center`, `box_size`, `exhaustiveness`, `n_poses`, `energy_range`, `cpu`, `seed` and `scoring`. Scoring may be `vina` or `vinardo`; this wrapper does not implement AutoDock4 map input. Box sizes must be positive and finite, box centers finite, pose count and exhaustiveness positive integers, energy range positive and finite, CPU count a nonnegative integer and seed a signed 32-bit integer. A CPU count of zero delegates thread selection to Vina. Invalid settings stop the job with an explicit error.

For example, `{"exhaustiveness": 16, "seed": 17}` preserves the preset box and other defaults. Each new options object merges with the original defaults, not a previous run's overrides. File paths provided as inputs replace their defaults; an explicitly blank path fails instead of silently using the example.

## Interpret the output

Inspect the 3D pose and download the PDBQT files alongside the score table. Lower scores rank poses within the selected scoring function and run. They are not measured binding affinities or proof that a compound binds. RMSD lower/upper bounds compare poses with this run's best pose; they do not measure agreement with a crystal structure. The legacy JSON field `affinity_kcal_mol` is retained for compatibility and contains the docking score.

Accuracy depends on the target, structure preparation, search region and sampling. A plausible-looking pose alone does not establish biological activity. Validate the protocol using independent target-specific evidence before making selection decisions. The [official FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html) describes these limitations and the need to hold all inputs and settings fixed for seed-based reproducibility.

## Outputs and provenance

`pose_summary` contains pose ranks, scores, RMSD bounds and file paths. `docking_summary` contains the resolved settings and interpretation limits. `structure_artifacts` lists the original docking output, individual poses, merged PDB, configuration, summary files and logs. `run_metadata` records the actual input and executable checksums, seed, options, runtime location and `status: completed` or `status: error`.

The Python modules use Biosimulant 0.0.34 and Python 3.12. Both are finite, once-before-run modules; the duration and communication interval do not change the number of docking jobs. Managed mode downloads the official versioned Vina binaries. External mode requires preinstalled `vina` and `vina_split`. No GPU is required.

The code, documentation and source-file notices supplied with the repository remain applicable. `source-provenance.json` pins the bundled input bytes without claiming independent experimental validation. See `MTS.md` for the execution and validation contract.
