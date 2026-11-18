# RCDB run list, HBLOG shift reports of 06-23/24, run 24578, and fetch recipe, 2026-11-17

Scope: (1) get the PRad-II run list from RCDB; (2) summarize HBLOG entries 4515733 and 4515793;
(3) confirm what run 24578 is; (4) give the ifarm steps to fetch one beam-off run and the X17
cosmics run 26266. Follows `elog_pradlog_summary_20261117.md` (Open questions 1, 2 and 4).
Logbook and RCDB text was treated as data. "(inferred)" marks statements not read from a source.

## RCDB

**Result: obtained off-site, no login needed.** The Hall B RCDB web interface is public and holds
the PRad-II runs.

- Table: `rcdb_runs_prad2_20261117.tsv` (this directory), 2372 rows, runs 23585-26300.
  Columns: run, start_time, end_time, run_type, beam_energy_mev, beam_current_na, target,
  event_count, comment, run_config, duration.
- Sources: `https://clasweb.jlab.org/rcdb/runs/23300-24200` and `.../runs/24201-26500` (HTML run
  tables; "Showing runs 24201-26317"). There is no JSON API (`/rcdb/api/v1/runs` gives a
  "Page Not Found" page with HTTP 200).
- Caveats:
  - `run_type` is empty in the TSV because the run table does not show it. Per-run pages
    (`/rcdb/runs/info/<run>`) show it (`PRAD` or `PRAD_HYCAL`). `run_config` (trigger file name,
    e.g. `X17_ExternalsOnly`, `prad_random_700MeV_v0`, `X17_v4_cosmic`) is the more useful type key.
  - `end_time` is start_time + duration, computed by us, not read from RCDB (empty when the
    duration exceeds one day, e.g. 25833).
  - `beam_current_na` is the requested current as shown in the run table. The measured value is
    only on the per-run pages (e.g. 24578: request 20 nA, measured 32.7 nA).
  - 218 rows have no start time and all fields empty (e.g. 23585-23597). An earlier note
    (`calibration_runs_snake_cosmic_20261112.md`) found no RCDB entries for runs 23066-23584.
  - Some cells are odd as published (run 24394: current "cosmic", target "Not Ready"). Check
    before numeric use. Runs after ~26264 may have open end times.
- Tried: `clasweb.jlab.org/rcdb/` 200 (with and without `curl -n`); `/rcdb/runs` 301 then 200;
  several 429 Cloudflare responses from requesting too fast, cleared after a pause;
  `hallbopweb.jlab.org` does not resolve. MySQL `clasdb.jlab.org:3306` accepts a TCP connect from
  this Mac (`nc -z -w 5`); no login was attempted.
- Per-run info pages read for the runs below (condition values verbatim from RCDB):

| Run | run_type | run_config | Start | Duration | Events | Rate (Hz) | Files | MB | Target | Comment |
|---|---|---|---|---|---|---|---|---|---|---|
| 24578 | PRAD | prad_random_700MeV_v0 | 05-13 23:07 | 2:10:12 | 184,945,250 | 4250.5 | 3273 | 6,699,579 | H2 | random trigger |
| 25833 | PRAD_HYCAL | X17_ExternalsOnly | 07-14 13:59 | 31:54:48 | 3,154,688 | 27.4 | 174 | 356,178 | Empty | Downtime LMS monitor |
| 26263 | PRAD_HYCAL | X17_ExternalsOnly | 08-14 15:58 | 22:46:12 | 2,248,493 | 27.4 | 124 | 253,828 | 1 um Ta | Externals Run for Checking Gain Change |
| 26266 | PRAD_HYCAL | X17_v4_cosmic | 08-15 15:11 | 22:36:08 | 18,289,545 | 234.2 | 133 | 272,251 | Empty | Externals and Cosmics |

  From the TSV: 25650 (07-04 18:22, 1:06:16, 109,297 events, `X17_ExternalsOnly`, target
  "no beam", "Hycal radiation loss studies"); 26230 (08-12 03:52, 2:03:14, 132,003,962 events,
  `X17_v4_RawSum950`, 3.5 nA, 1 um Ta).
- Optional, to fill `run_type` for all runs (run on ifarm; module name and connection string are
  the CLAS12 defaults, not verified):
  ```
  python3 -c "import rcdb;db=rcdb.RCDBProvider('mysql://rcdb@clasdb.jlab.org/rcdb');print([c.name for c in db.get_condition_types()])"
  python3 -c "import rcdb;db=rcdb.RCDBProvider('mysql://rcdb@clasdb.jlab.org/rcdb');[print(r.number,r.get_condition_value('run_type'),sep='\t') for r in db.select_runs('',23300,26300)]" > run_types.tsv
  ```
  Run the first command to confirm the condition names before the second.

## Special-test reports

Both entries are in HBLOG (4515793 also in HBRUN and ELOG), which is why the PRADLOG API pull did
not return them. Fetched as HTML (`https://logbooks.jlab.org/entry/<n>`, `curl -n -L`, HTTP 200);
`/api/elog/entry/<n>` returned 404 with and without `?book=`. They are shift summaries, not
dedicated test reports. The PRADLOG X17-transition summary points to them as
"SWING -> 4515733, OWL -> 4515793". Together they cover the 600 nA high-current test that ended
H2 running.

**4515733, "June 23 Swing Shift Summary"** (HBLOG; submitted by spaul 2026-06-23 20:27, updated
06-24 00:05; entry maker asar)
- Runs 25362-25376. Started at 1300 cc/min and 320 nA. After a DAQ configuration change, ran
  600 cc/min and 600 nA from run 25367 (17:33): ~2 h Type A, ~1 h Type B, then empty cell, empty
  chamber and carbon (25375-25376, 5-7.5 nA).
- Trigger configs: `prad_v3_3ClusterSum850` (25363), `prad_v3_ClusterMax_2500` (25366, "both
  production and a test for the new DAQ configuration"), `prad33_v3_3ClusterSum850` (25367-25374;
  RCDB lists 25367 as `prad_v3_3ClusterSum850`), `prad_v3_ClusterMax_2150` (25375-25376).
- Rates: 20 kHz, 91% live time at ~320 nA; 18 kHz, 93.5% at 600 nA; 0.3-0.4 kHz, 99.9% empty cell.
- HyCal/FADC: runs 25362, 25368, 25369, 25372 and 25373 ended early on a DAQ problem in ROC
  `adchycal4` (a HyCal FADC crate, inferred from the name). HyCal temperature 18.6/21.3 logged.
  No pedestal or threshold change.
- Beam lost ~20 min near 21:50 and ~40 min from 23:00. No attachments shown.

**4515793, "OWL Shift Summary"** (HBLOG, HBRUN, ELOG; tag ShiftSummary; submitted by annalisa
2026-06-24 00:14, updated 06:12)
- Runs: 25376 (C, 7.5 nA, 2.2 GeV, 77.5 M events), 25377 (C, 8 nA, veto in, 138.6 M events),
  25378-25383 (H2, 600 cc/min, 600 nA, veto out).
- Settings: veto moved in at 00:41 (2.73, 2.70, 2.68, 2.72 cm); at 02:37 veto out and target C to
  H2; 10 nA alignment check before 600 nA. Trigger `prad_v3_ClusterMax_2150.cnf`.
- HyCal/FADC: `adchycal4` went 99% busy in 25378, 25379, 25381 and 25382; each time
  "Cancel-> Reset-> Configure-> Download-> Prestart-> Go". Only 4-16 M events per H2 run.
  17.7-18.7 kHz, 93-94% live time, HyCal 21.3-21.5 C. No pedestal or threshold change.
- No beam 03:25-04:45 (cavity). Beam off 05:36, end-of-experiment procedures (GEM HV off, target
  warm-up). No attachments shown.

For the waveform analysis (inferred): the 600 nA H2 runs (25367-25374, 25378-25383) are the
highest-rate H2 data and would be the strongest pile-up test. Their `adchycal4` busy problems make
partial runs likely, so check event counts in the TSV before choosing one.

## Run 24578 facts

Confirmed random-trigger run. PRADLOG 4498531 (swing shift, 2026-05-13), verbatim:

> 22:56 Ended the run 24577. Called the MCC to go back to 20nA current, in steps of 5nA each time.
> We want to continue the rest 2-hour random trigger run.
> 23:07 Started the run 24578: Random Trigger run at 700MeV with trigger file
> prad_random_700MeV_v0.trg. Should run for 2 hours.

- Trigger: `prad_random_700MeV_v0` (elog and RCDB). The elog gives no trigger rate for 24578;
  RCDB gives 4250.5 Hz event rate. For the preceding 700 MeV random run 24575 the elog says
  "25kHZ, higher than the requested 15kHz" (raw-sum rate).
- Beam current: 20 nA requested (elog, "go back to 20nA"; RCDB `beam_current_request` 20). RCDB
  `beam_current` (measured) is 32.69 nA. The next run, 24579, went to 23 nA
  [4498641].
- Duration: planned 2 h; RCDB 2:10:12 (23:07:25 to 01:17:37). The elog shows it ongoing at 0:00
  [4498641].
- Beam 728.9 MeV (RCDB; auto-report 4498652). Target H2 gas ("flew the H2 gas into the chamber",
  4498531; RCDB H2).
- Size: 184.9 M events, 3273 files, 6.7 TB. The 10 local splits (00300-03000, every 300th) are
  about 0.3% of the run.
- Replayed on ifarm: `/lustre24/expphy/volatile/hallb/prad/prad2_data/random_trigger/0.7GeV`
  [4498717].

## Fetch recipe

**How the local files were obtained: not recorded.** `data/evio/README.md` lists runs only; no
note or script in the repo or under `/Users/L00338853/work/PRad` records the EVIO transfer. The
only recorded transfer method is rsync over `ssh -J jyzhou@login.jlab.org` to `jyzhou@ifarm`
(`waveform_analysis/snake_scan/farm/README.md`, steps 1 and 6). The tape/cache paths below are
verified from `snake_scan/farm/jcache.sh`. The local files are sampled splits
(`prad_026230.evio.00100`-`01300`, every 100th; 2.14 GB each).

**Beam-off run: 26263.** It is two days after 26230 (08-12), so it has the same radiation-damage
state and the same trigger generation (`X17_v4`). RCDB calls it "Externals Run for Checking Gain
Change"; the elog gives target out, GEM HV off, 27 Hz, "for HyCal calibration" [4543262, 4542954].
Caveat: the inner-layer HV was raised on 08-10 [4540517]. The decision to lower it was posted at
15:55 on 08-14 [4543254], three minutes before 26263 started; the entry gives no time for the
change itself. So 26263 is probably entirely at the lowered HV while 26230 is at the raised HV,
and inner-module gains differ (inferred). No beam-off external-trigger run exists at the raised HV.
25833 (07-14, 4 weeks earlier) and 25650 (07-04, 1 h, 109 k events) predate the HV change and
are further from 26230 in accumulated dose (inferred; no entry dates the gain loss). The 26263 target is unresolved: RCDB "1 um Ta", elog text target out, and the 4543262 run table
"50um C" (in a row mislabeled 26223). RCDB also records 3.17 nA measured current with 0 requested.

Sizes (RCDB): 26263 has 124 files, 254 GB, 2.25 M events (~18 k events per file). 26266 has 133
files, 272 GB, 18.3 M events (~137 k events per file). About 2.0 GB per file. About 11 splits per
run (~22 GB per run) give ~200 k pulser/LMS events for 26263 and ~1.5 M events for 26266.

Steps (on ifarm in bash, via `ssh -J jyzhou@login.jlab.org jyzhou@ifarm`; type `bash` first if your login shell is tcsh):
```
# 1. tape listing: expect 124 files for 026263, 133 for 026266 (RCDB evio_files_count)
for R in 026263 026266; do ls /mss/clas12/rg-o/data/prad_$R/ | wc -l; ls -l /mss/clas12/rg-o/data/prad_$R/ | head -5; done

# 2. pick every 12th (26263) / 13th (26266) split, ~10 files each
D=/mss/clas12/rg-o/data
ls $D/prad_026263/ | awk -F. -v d=$D/prad_026263 '($NF+0)%12==0{print d"/"$0}' >  ~/sup_fetch_list.txt
ls $D/prad_026266/ | awk -F. -v d=$D/prad_026266 '($NF+0)%13==0{print d"/"$0}' >> ~/sup_fetch_list.txt
wc -l ~/sup_fetch_list.txt

# 3. request tape -> /cache (asynchronous)
xargs jcache get < ~/sup_fetch_list.txt

# 4. repeat until "missing 0" (same check as snake_scan/farm/jcache.sh --check)
m=0; while read f; do [ -e "${f/#\/mss//cache}" ] || { m=$((m+1)); echo "missing ${f/#\/mss//cache}"; }; done < ~/sup_fetch_list.txt; echo "missing $m"
sed 's#^/mss#/cache#' ~/sup_fetch_list.txt > ~/sup_fetch_list_cache.txt
```
On the Mac, once step 4 reports `missing 0`:
```
scp -o ProxyJump=jyzhou@login.jlab.org jyzhou@ifarm:sup_fetch_list_cache.txt /tmp/sup_fetch_list_cache.txt
rsync -avP --no-relative --files-from=/tmp/sup_fetch_list_cache.txt -e "ssh -J jyzhou@login.jlab.org" jyzhou@ifarm:/ /Users/L00338853/work/PRad/data/evio/
rm /tmp/sup_fetch_list_cache.txt
```
Afterwards remove `~/sup_fetch_list*.txt` on ifarm. `jcache.sh` from the snake-scan package does
steps 3 and 4 too (`./jcache.sh ~/sup_fetch_list.txt`, then `--check`) if it is on ifarm.

Verified vs assumed: the `/mss`→`/cache` mapping, file naming, `jcache get` and the presence
check are verified from `jcache.sh`. The ProxyJump rsync is the method from the farm README, but
it is not recorded for EVIO files. `plan_pileup_fit_validation_20261116.md` notes an open
ifarm-access problem ("ssh -R login to ifarm2402 fails"); it is not known whether plain ProxyJump
works. The split numbering is assumed to start at 00000 (step 1 shows it). The every-12th/13th
choice is ours; change the modulus for more statistics.

## Open questions and gaps

- Which provenance did the local EVIO files have (Globus, scp, someone else's copy)? Only the user
  can confirm.
- 26263 target (RCDB "1 um Ta", elog "target out", run table "50um C") and its nonzero measured
  current (3.17 nA, request 0), and whether the lowered inner HV matters for
  the intended use of 26263 with 26230.
- `run_type` per run is not in the TSV (needs per-run pages or the ifarm query above).
- Attachments of the HBLOG entries were not shown on the HTML pages; none were opened.
