# Local test results (2026-11-14, Mac, no farm access)

Inputs: local EVIO `prad_023701.evio.{00000,00010,...,00080}` (9 files, run 23701) and `prad_023744.evio.{00000,00030,00060}` (3 files),
the only files of runs 23699-23723/23744 on this machine. Every file was processed with `farm_job.sh <evio> <dir> --keep-every 1`
(the exact configuration; at most 2 files at a time, machine load 3-7); `farm_merge.py --in <dir> --out <dir2> --lms-runs 23701`;
`python test_farm.py <dir> <dir2>` compares with `SS/out` and the wave-1 caches. Result: **10/10 checks passed**.

## Cost per file (`/usr/bin/time -l`, one process, includes the second LMS decode)

| file | wall s | max RSS MB | product MB (keep-every 1) |
|---|---|---|---|
| 23701.00000 | 538 | 618 | 19.6 |
| 23701.00010 | 323 | 442 | 9.5 |
| 23701.00020 | 515 | 585 | 17.2 |
| 23701.00030 | 415 | 511 | 16.4 |
| 23701.00040 | 592 | 676 | 19.7 |
| 23701.00050 | 429 | 563 | 19.2 |
| 23701.00060 | 411 | 624 | 18.6 |
| 23701.00070 | 399 | 547 | 17.9 |
| 23701.00080 (short) | 358 | 283 | 3.4 |
| 23744.00000 / 00030 / 00060 | 95 / 96 / 92 | 535 / 488 / 463 | 0.50 / 0.38 / 0.47 |

Default `--keep-every 2`, 23701.00000 (run alone with another job): 580 s, 525 MB, 9.8 MB. The user time is about 360 s; the remaining
time is system time of the Python/numpy per-channel loop (load dependent). The LMS second decode: 2 s (`lms_pass` on 23701.00040,
705,581 events). Requests (2 h, 2 GB for scan; 1 h, 1.5 GB for 23744) are 3-4x the measured maxima. Merge of the 12 products:
2.1 s, 447 MB (9 modules; about 0.1 s per module, so roughly 10 min for about 1100 modules plus reading 3-5 GB).
A truncated EVIO file (5 MB of a real one) makes the C++ decoder segfault (exit 139): the job fails and leaves no product, as it should.

## Comparison with wave 1 (tolerances stated; "exact" means array_equal)

| check | result |
|---|---|
| dwell rows (`blocks_to_rows` on the stored blocks/EPICS) vs `out/dwell_23701.tsv`, 9 files, 17 rows | exact (event ranges, ts, module, n_events, kind); share and ptrans within 5e-5 |
| LMS/alpha/pedestal series (`L_*` keys, 8 arrays x 9 files) vs `output/snake_scan/lms_npz/lms_23701_*.npz` | exact, 72/72 arrays |
| edge table (`edge_scan.tsv`) vs `out/edge_23701.tsv`, 9 modules | max abs d s 2.0e-5, d s_err 1.0e-5, d s_central 7.0e-5; n_single max rel d 4.8e-4, n_pileup 6.0e-4. W580 pass 1: n_single 341605 and n_pileup 280293 identical to wave 1, s_pass1 0.95925 both. The count differences (up to ~170 events of 280k) are the events skipped while the running beam module warms up (first 100 large-peak events of a file) |
| pooled `lms_scan.tsv` (23701 only) vs `out/lms_scan.tsv`, 1151 modules | exact (L and err) |
| `transfer_2p24GeV_H2.tsv`, `transfer_3p49GeV_H2.tsv` (via unchanged `transfer.build`) vs `out/` | s, s_err within 2e-5 / 1e-5 (5-decimal printout); L_scan, L_prod, g_pool exact; g_T and r within 2.5e-5 rel; sigma_T within 8.4e-4 rel (s_err rounding). Pool spread 2.28 % (2.24 GeV) and 3.03 % (3.49 GeV), as wave 1 (2.29 %, 3.03 %) |
| `prior_<rt>.tsv` vs `out/` | same 2 modules, max abs d g 1.5e-5 / 5e-6 |
| `mip_23744.tsv` vs `out/mip_23744.tsv` (3 files, 6 groups, all numeric columns) | exact to the printout (max rel 0) |
| `lms_23744.npz` vs `out/lms_23744.npz` (19 arrays) | exact |

Which comparison is exact and which is pooled: the **primary** wave-1 transfer uses L_scan from 23701 alone: `--lms-runs 23701` reproduces it
(table above). The default of `farm_merge.py` pools ALL scan runs found; locally only 23701 is a scan run (23744 is the cosmic run and is
not pooled), so here the pooled `lms_scan.tsv` and `transfer_*.tsv` are byte-identical in their data lines to the 23701-only ones
(checked with diff). The wave-1 "merged" variant (23701 + 23744 per module within 2 SE, `*_merged.tsv`) is not produced by the farm merge
(23744 appears as its own row in `lms_scan_per_run.tsv`: median ratio to the pool 0.9935, err 6e-4); on the farm the pooled all-scan
`lms_scan.tsv` will differ from the 23701-only primary and has to be judged with `lms_scan_per_run.tsv`.

## Other checks

- Default `--keep-every 2` (23701.00000 alone, W599): pass-1 s 1.02979 vs wave 1 pass-1 1.02981; W580 emulated from the keep-every 1 product
  (events with even number): s 0.95897 vs 0.95925 all events, fit error 1.07e-3 vs 0.76e-3.
- Pile-up settings on the 9 files (`--pileup`): s(w1) / s(central) / s(none) = W547 0.9874/0.9899/0.9991, W552 1.0970/1.0985/1.1165,
  W564 1.1081/1.1106/1.1658, W599 1.0289/1.0310/1.0419 (n_single 352002/390051/661532): the setting moves s by +0.2 % (central) and +1-5 % (none).
  The products serve all three without re-decoding.
- Synthetic multi-run merge (copies of products relabelled as runs 23702-23704, 23721 lead glass): per-run table and pooled `lms_scan.tsv`
  work; the lead-glass run is excluded from the edges (9 modules as before); 7 synthetic runs with a 3 % LMS shift in one are flagged
  (|z_L_rel_med| 301, listed in the header); degenerate (identical) runs need the error term in z, which is in.
- `--steps transfer` alone and `--tag t1 --lms-runs 23701` (writes `lms_scan_t1.tsv`, `transfer_<rt>_t1.tsv`) work.
- `check_golden.py` writes and reproduces `golden_23701_00000.json`; a perturbed golden value gives MISMATCH.
- `make_filelist.sh` on a fake `/mss` tree (3 runs of 11 splits, 23744 with 5, one missing run directory, a stray file): every 3rd split
  (00000/03/06/09), all 23744 splits, stray file ignored, missing directory reported with a non-zero exit; `jcache.sh` with a stub `jcache`
  (batches), `--check` counting present/missing; `submit_swif2.sh` and `submit_sbatch.sh` generate the lines (not executed: no swif2/slurm
  here); `farm_missing.sh` lists the files without product and deletes stale partial files.
- Not tested (no farm): swif2/sbatch option syntax, project/account names, /cache path mapping, the ifarm Python/bindings build.

## Residue

All temporary products, logs and fake trees were under the session scratchpad with the prefix `sup-d1e2ea-` and were removed; nothing
was written outside `SS/farm/` (and `__pycache__` directories next to existing modules).
