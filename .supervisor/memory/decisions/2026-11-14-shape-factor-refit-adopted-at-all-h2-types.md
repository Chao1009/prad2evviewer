---
domain: science
summary: The shape-factor refit after gain changes is adopted at all four H2 types; the 2.24 GeV tail_lo_f guard is waived.
affects: [analysis/pyscripts/energy_resolution_benchmark.py]
verified: 2026-11-14
source: user 2026-11-14; analysis/pyscripts/notes/plan_energy_resolution_improvement_20261112.md
supersedes:
---
The refit is on for every H2 type, and the 2.24 GeV tail_lo_f guard does not block it.
Why: the user accepted the out-of-sample gain over the guard's tail penalty.
Rejected: adopting only at the types that pass the guard.
