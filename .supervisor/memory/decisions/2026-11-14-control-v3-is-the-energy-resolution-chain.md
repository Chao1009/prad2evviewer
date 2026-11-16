---
domain: science
summary: Control v3 (the benchmark defaults) is the current energy-resolution chain; new packages are compared against it.
affects: [analysis/pyscripts/energy_resolution_benchmark.py]
verified: 2026-11-14
source: user 2026-11-14; analysis/pyscripts/notes/plan_energy_resolution_improvement_20261112.md
supersedes:
---
The benchmark defaults define control v3, the reference chain for energy-resolution work. Out-of-sample CB sigma/E is 3.160, 1.746, 1.470, 1.443 and 1.987 % at 0.73, 2.24, 3.49 (old target) and 3.49 (Ta X17) GeV.
Why: it includes the adopted shape-factor refit, per-energy I(phi) correction and ring-scoped gain constants, each validated out of sample.
Rejected: earlier controls (v1, v2) and the legacy chain, kept only as comparison arms through the documented `--set` reproduction string.
