---
domain: code
summary: Where the exotic-waveform study's scripts, outputs, write-up chapter, opt-in hook and caches live, and what each cache regenerates.
affects: [analysis/pyscripts/waveform_analysis/pileup_analysis/exotic_events/, analysis/pyscripts/output/pileup_rates/, analysis/pyscripts/waveform_analysis/writeup/figs/exotic_events/]
verified: 2026-11-18
source: analysis/pyscripts/notes/plan_pileup_fit_validation_20261116.md#exotic-events-impact-packages-approved-2026-11-18
supersedes:
---
Scripts exotic_p1…p13 and outputs sit in analysis/pyscripts/waveform_analysis/pileup_analysis/exotic_events/; exotic_summary.py rebuilds exotic_summary_20261118.pdf (253 pp, uncommitted). Chapter sec:exotic is in writeup/waveform_analysis_prad2_x17.tex, figures in writeup/figs/exotic_events/. EXOTIC_HOOK in energy_resolution_benchmark.py and moller_calib.py is opt-in; the default path is unchanged.
Caches in analysis/pyscripts/output/pileup_rates/, all regenerable: exotic_sample/ (sample scan); exotic_full/ (19 Ta npz, H2 candidate pickles without amplitudes); exotic_p7_cache.npz, exotic_p8_cache.npz (cluster join, energy error); exotic_p9/ (calibration arms); exotic_p10_cache.pkl (yields); exotic_p12_chan_cache.npz (slice rates); exotic_p13/ and exotic_p13_cache_v2.npz (drop policy); exotic_p13_cache.npz is superseded.
Hot modules: W335, W613, W457; high-gain W398.
