# prad2py: timing algorithm and log-normal CFD fit

Added 2026-11-06 (python/bind_dec.cpp only; C++ analyzer unchanged, CFD fraction stays 0.5).

    source analysis/pyscripts/setup_env.sh      # PYTHONPATH=build/python
    from prad2py import dec
    wa = dec.WaveAnalyzer(core.daq_cfg.wave_cfg)         # or p.wave_ana from _common.setup_pipeline
    ped_mean, ped_rms, peaks, fits = wa.analyze_fit(samples_uint16)

- `peaks[i]` is the usual `dec.Peak`; new field `pk.time_algo`: 0 peaking sub-sample,
  1 linear CFD, 2 log-normal fit CFD. `pk.time` is the time in ns.
- `fits[i]` is a `dec.LogNormalFitResult` aligned with `peaks[i]`, read-only fields:
  `ok, t_cfd_sample, chi2_per_dof, A, t0, mu, sigma`. `ok=False` (zeros) when no fit was done.
  For `time_algo == 2`: `pk.time == t_cfd_sample * 1e3 / cfg.clk_mhz`.
- `analyze` and `analyze_full` are unchanged.

Build (incremental, existing tree; interpreter = miniforge python 3.13):

    cd build && make prad2py        # ~19 s; writes build/python/prad2py.cpython-313-darwin.so

The build tree only produces the 3.13 module; the 3.14 .so is stale and was not rebuilt.
