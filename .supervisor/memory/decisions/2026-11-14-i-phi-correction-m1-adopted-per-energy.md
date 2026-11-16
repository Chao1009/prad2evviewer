---
domain: science
summary: The I(phi) correction method M1 is adopted, fitted separately per beam energy.
affects: [analysis/pyscripts/energy_resolution_benchmark.py]
verified: 2026-11-14
source: user 2026-11-14; analysis/pyscripts/notes/plan_energy_resolution_improvement_20261112.md
supersedes:
---
M1 of the I(phi) port to mf_pulse.py is adopted with per-energy constants. The ring-scoped gain constants are adopted with the gap clause of their gate waived, and the outermost-layer cut is not applied.
Why: the gains held out of sample; the user judged the waiver and the cut not worth their cost.
Rejected: a single constant for all energies; the outermost-layer cut.
