flowchart TD
  A["Raw FADC waveform<br/>100 samples x 4 ns"] --> B["Pedestal<br/>first 30 samples, iterative clip: mean, rms"]
  B --> C1["WaveAnalyzer (C++)<br/>smooth, local max, walk-out integral,<br/>linear CFD 0.5 + log-normal time"]
  B --> C2["Matched filter MF (Python, mf_pulse.py)<br/>sampled template bank, up to 4 pulses: A, onset t0"]
  C1 --> D["Pulse list: integral (rawsum)<br/>or A x I_int x f(A) x kappa (MF); time"]
  C2 --> D
  D --> E["Module energy<br/>calibration file to MeV, LMS gain correction,<br/>per-crystal gains"]
  E --> F["Clustering (HyCalCluster)<br/>legacy: largest pulse in pre-window, no time cut<br/>gated: |t - t_seed| <= W"]
  F --> G["Cluster energy"]
  G --> H["Position correction<br/>E to E / R(u,v): bilinear run-type map or hybrid per-module poly8"]
  H --> I["Elastic selection<br/>seed r < 200 mm, seed E cut, row/col exclusion,<br/>elastic window, seed-indexed gain refinement"]
  I --> J["Resolution fit<br/>Gaussian (iterative +-2 sigma) and Crystal Ball"]
  J --> K["Per-split sigma_E/E"]
  K --> L["Per-module medians, run-type means,<br/>2-term fit a, b (headline)"]
