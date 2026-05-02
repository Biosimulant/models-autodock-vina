# models-autodock-vina

> Storage-only repo: each former root model now lives in `labs/<slug>/model/` and is wrapped by
> `labs/<slug>/lab.yaml`. This repo has no repo-level import catalog and no composed labs at the root.

Curated collection of **AutoDock Vina-family docking models** for the **biosim**
platform.

This repository currently ships one native Python `biosim.BioModule` wrapper:
`vina-autodock-vina-docking-predictor`, a single-complex AutoDock Vina docking
module for prepared receptor and ligand `PDBQT` inputs.

## What's Inside

### Wrapper Sublabs

| Sublab | Description |
|---|---|
| `vina-autodock-vina-docking-predictor` | Native AutoDock Vina `1.2.7` wrapper for single receptor/ligand docking workflows. |

## Scope

This repository is for:
- native AutoDock Vina-family wrappers that implement the `biosim.BioModule` contract
- single-complex docking runs that emit ranked pose summaries plus file-backed structural artifacts
- portable examples that can be exported to `.bsilab` and validated against remote CPU execution

This repository is not for:
- receptor or ligand preparation from `PDB`, `SDF`, or `SMILES`
- AutoDock-GPU workflows
- virtual screening or batch scheduling orchestration

## Remote Execution

The Vina model uses the existing generic remote execution path:

- the wrapper bootstraps a managed native runtime by downloading the pinned official AutoDock Vina `1.2.7` release binaries for the current OS and architecture
- remote runs force `runtime_mode: managed` through manifest-declared remote init overrides
- the wrapper creates its managed runtime under the mounted remote cache/work roots and runs classic CPU-backed Vina there
- the wrapper streams live CLI stdout and stderr, emits `BSIM_PROGRESS:` milestones, and derives a merged `top_rank_complex.pdb` artifact for the existing `structure3d` renderer

The release-grade validation target is Linux CPU execution on Modal.

## Examples

See [examples/README.md](examples/README.md) for the example inventory, including
the remote `.bsilab` builder used for desktop end-to-end validation.
