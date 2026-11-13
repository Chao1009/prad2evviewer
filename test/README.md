# Test & Diagnostic Tools

Generic EVIO diagnostic command-line tools.  GEM-specific tools live
in [`gem/`](../gem/README.md).

## Installed tools

### evio_dump

EVIO file structure inspector.

```bash
evio_dump <file> [-m mode] [-n N] [-D daq_config.json] [-v]
```

Modes: (default) summary by tag, `tree`, `tags`, `epics`, `vtp`
(VTP banks), `rf` (RF time TDC hits), `event` (detail for event *N*),
`triggers`, `trig-debug` (event tags vs. TI type vs. FP trigger bits),
and `bank-debug` (bank structure / dispatch coverage).  `-v` prints
every event in the `triggers` / `trig-debug` modes.

### ped_calc

Compute HyCal per-channel pedestals from trigger-selected events.

```bash
ped_calc <evio_file> [-D daq_config.json] [-t bit] [-o file.json] [-n N]
```

Default trigger bit is 3 (LMS_Alpha / pedestal for PRad).

## Dev-only tools (not installed) — sources under `test/dev/`

Built alongside the installed tools but kept out of the install tree;
useful from a build tree only.

| Tool | Purpose |
|---|---|
| `ts_dump` | Dump TI timestamp and trigger info per event. |
| `livetime` | DAQ live-time calculator (DSC2 scalers + pulser).  Accepts a single file, a base name (auto-finds `.00000`, `.00001`, …), or a directory. |
| `dsc_scan` | Walks every 0xE115 DSC2 scaler bank, identifies the parent ROC/crate, parses the per-slot 67-word DSC2 layout, and reports per-channel gated/ungated counts plus the implied live time.  The physics-trigger and reference-clock channels both expose a live-time, so the tool prints both and lets the user pick. |

```bash
ts_dump   <file> [-n max_events] [-D daq_config.json]
livetime  <input> [-D daq_config.json] [-f freq_hz] [-t interval_sec]
dsc_scan  <input> [-D daq_config.json] [-N n_events] [--all]
```

### ET-specific dev tools (`-DWITH_ET=ON`)

| Tool | Purpose |
|---|---|
| `evc_scan` | Three-mode smoketest — read evio buffers, scan with detail, or connect to an ET station. |
| `et_feeder` | Replay an evio file into an ET ring at a controlled rate. |
| `evet_diff` | Read an evio file and an ET ring in parallel and diff their raw buffers; pairs with `et_feeder`. |
