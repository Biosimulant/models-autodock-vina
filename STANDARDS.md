# Model & Space Standards

Acceptance criteria and conventions for every model and example in this
repository.

New contributions must satisfy all items marked `REQUIRED`. Items marked
`RECOMMENDED` are strongly encouraged and will be enforced over time.

## Model Standards

### Directory Layout

`REQUIRED` for every model:

```text
models/<family>-<subject>-<role>/
├── model.yaml
├── README.md
├── src/
│   └── <module>.py
└── tests/
    ├── conftest.py
    └── test_<module>.py
```

### Manifest (`model.yaml`)

`REQUIRED` fields:

```yaml
schema_version: "2.0"
title: "<Family>: <ComponentName>"
description: "<One-line description>"
standard: other
tags: [boltz, ...]
authors: ["Biosimulant Team"]
biosim:
  entrypoint: "src.<module>:<ClassName>"
```

Rules:
- `standard` must be `other`
- `biosim.entrypoint` must resolve to an importable callable
- any declared package dependency must use exact `==` pinning
- external runtime requirements that are not Python wheels must be documented
  in the model README

### Source Code

`REQUIRED` for every Python file in `src/`:
- SPDX header
- module docstring
- `from __future__ import annotations`
- type annotations on public methods
- a `BioModule` subclass implementing the BioSim contract

### Boltz-Specific Runtime Rules

`REQUIRED`:
- large structural outputs must stay file-backed rather than being embedded
  directly into manifests or checked into the repo
- models must document managed-runtime behavior plus any optional external
  executable and hardware/runtime assumptions
- subprocess boundaries must be explicit in code and tested
- example files must stay aligned with the current constructor parameters,
  input ports, and output ports

`RECOMMENDED`:
- treat Boltz-family wrappers as run-once modules using cached outputs on later
  ticks
- emit compact structured summaries plus absolute artifact paths
- record command, output directory, and stderr/stdout in run metadata
- emit `structure3d` visuals when a model produces `mmCIF` or `PDB` artifacts
  so compatible BioSim clients can inspect the structure directly

### Tests

`REQUIRED` for each model:
- default instantiation
- input and output port checks
- output keys match `outputs()`
- reset behavior
- mocked subprocess success path
- mocked subprocess failure path
- missing expected output artifact handling
- repeat `advance_to(t)` does not rerun with unchanged inputs

## Example Standards

Examples live under `examples/` and must:
- use the actual public module interface
- avoid ONNX terminology unless the model actually uses ONNX
- prefer real, runnable small examples over placeholder-only examples
- document what must exist locally before the example can run

Example folders should include:
- a runnable YAML config or `space.yaml`
- brief documentation of expected outputs

## Acceptance Checklist

Before merging:
- `python scripts/validate_manifests.py`
- `python scripts/check_entrypoints.py`
- `bash scripts/check_public_boundary.sh`
- `pytest`
- example files reviewed against the current module interface
