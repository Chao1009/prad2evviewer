---
domain: science
summary: Snake-scan analysis always uses the LMS ratio (HV changed), no tagger, beam path taken from the data; the farm pass is approved.
affects: [analysis/pyscripts/notes/plan_p14_snake_scan_20261114.md]
verified: 2026-11-14
source: user 2026-11-14; analysis/pyscripts/notes/plan_p14_snake_scan_20261114.md
supersedes:
---
Snake-scan gains are normalised with the LMS ratio (user: "always use LMS ratio"), because HyCal HV changed after the scan. The tagger was not read out, and the beam path comes from the data. The user approved the farm pass ("Prepare the farm pass now").
Why: HV changes invalidate absolute gains; no tagger gives no external energy.
Rejected: absolute-gain normalisation; a tagger energy reference.
