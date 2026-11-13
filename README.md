# PRad-II Event Viewer & Monitor

EVIO decoder, event viewer, and online monitor for the PRad-II experiment
at Jefferson Lab.  Decodes FADC250 waveforms (HyCal), SSP/MPD data (GEM),
V1190 TDC (tagger), DSC2 scalers, and EPICS slow control, and serves a
web GUI for waveform inspection, clustering, gain monitoring, beam
status, GEM tracking efficiency, and physics replay.  The original
PRad ADC1881M format remains supported through DAQ configuration.

## Building

### Linux

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake ≥ 3.14, C++17.  The `evio` and `et` dependencies resolve from the
Hall-B CODA installation by default; if not found, CMake fetches them
from GitHub.  Other dependencies (`nlohmann/json`, `websocketpp`, `asio`)
are always fetched automatically.

Optional flags:

- `-DBUILD_ANALYSIS=ON` — ROOT-based replay and analysis tools (requires ROOT 6.0+).
- `-DBUILD_GUI=ON` — Qt remote client `prad2_client` (requires Qt6 or Qt5 WebEngine).
- `-DBUILD_PYTHON=ON` — pybind11 bindings `prad2py` (also enables the Python GEM tools).
- `-DPython_EXECUTABLE=<path>` — pin the interpreter explicitly (needed on systems with multiple `python3` versions in `PATH`, where CMake may otherwise pick the wrong one).
- `-DWITH_ET=OFF` — disable ET (live monitoring); required on Windows.
- `-DEVIO_SOURCE=fetch` — force fetching evio from GitHub (no CODA needed).

### Windows (MSYS2)

```powershell
# MSYS2 packages
pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake \
          mingw-w64-x86_64-ninja mingw-w64-x86_64-expat

# Build (PowerShell with MinGW in PATH)
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DWITH_ET=OFF
cmake --build build
```

A native WSL2 build is also supported.

## Server

`prad2_server` is a unified HTTP server for both file-based viewing and
online ET monitoring.  It supersedes the former `prad2_viewer` and
`prad2_monitor` executables.

```bash
prad2_server [evio_file] [-p port] [-H] [-i] [-f filter.json] \
             [-c monitor_config.json] [-r reconstruction_config.json] \
             [-D daq_config.json] [-d data_dir] [--et]
```

Opens a web GUI at `http://localhost:5051` with tabs for Waveform Data,
Clustering, Gain Monitoring (LMS), EPICS, Physics, and GEM tracking.
See [docs/SERVER_API.md](docs/SERVER_API.md) for the full
HTTP/WebSocket API reference and the interactive CLI commands
available with `-i`.

### File mode (default)

```bash
prad2_server data.evio -H                                          # PRad-II
prad2_server prad.evio -D database/prad1/prad_daq_config.json -H   # PRad
prad2_server -d /data/stage6 -H                                    # file browser
```

### Online mode

```bash
prad2_server --et                       # connect to ET; config from database/monitor_config.json
prad2_server data.evio --et -H          # start with file; "Go Online" button available in UI
```

ET connection settings are read from the `"online"` section of
`monitor_config.json`.

### Mode switching

The server maintains separate data accumulators for file and online
modes; switching modes does not discard data from the other mode.
Use **Clear All** in the UI to reset the active mode's data.

- Loading a file switches to file mode.
- The **Go Online** button in the web UI (or the `--et` flag) switches to online mode.
- The toggle button appears whenever both capabilities are available (compiled with `WITH_ET`).

Test with `et_feeder`:

```bash
et_start -f /tmp/test_et -s 100000 -n 500
./bin/prad2_server --et -D ../database/prad1/prad_daq_config.json \
                   -c ../database/prad1/prad_monitor_config.json \
                   -r ../database/prad1/prad_reconstruction_config.json
./bin/et_feeder prad.evio -f /tmp/test_et -i 50 -n 5000
```

### Remote access (JLab ifarm)

To run the server on a JLab ifarm node and view it in your local browser:

**1. Start the server on ifarm** (pick a port, e.g. 5678):

```bash
ssh username@ifarm2402
cd /path/to/prad2evviewer/build
./bin/prad2_server data.evio -H -p 5678
```

**2. SSH tunnel from your local machine:**

```bash
ssh -L 5678:ifarm2402:5678 -J username@scilogin.jlab.org username@ifarm2402
```

**3. Open in your local browser:**

```
http://localhost:5678
```

The same recipe works for online monitoring — start
`prad2_server --et -p 5678` on the ifarm node where the ET system is
accessible.

## Native event-by-event GUIs

For interactive per-event inspection, use the PyQt6 tools in `scripts/`
and `gem/`:

- `hycal_event_viewer` — FADC waveform browser plus live HyCal clustering (Waveform / Cluster tabs).
- `gem_event_viewer` — GEM strip hits, clusters, and live threshold tuning.
- `gem_hycal_match_viewer` — Per-event HyCal↔GEM matching, parametric matching cut, find-next search.

All three run standalone against an EVIO file; no server needed.

```bash
hycal_event_viewer       data.evio
gem_event_viewer         data.evio
gem_hycal_match_viewer   data.evio
```

## Qt thin client (optional)

Build with `-DBUILD_GUI=ON` (uses Qt6 WebEngine by default; falls back
to Qt5 if Qt6 is not available).

`prad2_client` is a Qt WebEngine wrapper that connects to a remote
`prad2_server` instance — useful when the web dashboard should live
in a native window rather than a browser tab:

```bash
prad2_client                            # connect to localhost:5051
prad2_client -H clonpc19 -p 8080        # connect to a remote server
```

## Tools

- [test/README.md](test/README.md) — generic EVIO diagnostics (`evio_dump`, `ped_calc`, plus dev tools in `test/dev/`).
- [gem/README.md](gem/README.md) — GEM tracker tools (`gem_dump`, `gem_event_viewer`, `gem_cluster_view`, `gem_layout`, …) and detector reference notes.
- [analysis/README.md](analysis/README.md) — ROOT-based replay, slow-control filtering, and physics analysis (binaries are installed with a `prad2ana_` prefix).
- [calibration/README.md](calibration/README.md) — HyCal gain-equaliser, snake-scan, and PMT response model.
- [calibration/cosmic_gain/README.md](calibration/cosmic_gain/README.md) — cosmic-ray HV iteration macros.
- [scripts/README.md](scripts/README.md) — HyCal scaler / pedestal / gain / coincidence / map-builder GUIs, HyCal event viewer, GEM↔HyCal matching viewer, tagger viewer, and `daq_tool/` (dev-only DAQ editors) plus `shell/` (operator scripts).
- [docs/technical_notes/README.md](docs/technical_notes/README.md) — citation-quality writeups of the waveform, HyCal, and GEM clustering algorithms, each with a Python reference implementation.

## Installation

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/prad2
cmake --build build -j$(nproc)
cmake --install build
```

Set up the environment (sets `PATH`, `LD_LIBRARY_PATH`, `PYTHONPATH`,
`PRAD2_DATABASE_DIR`, `PRAD2_RESOURCE_DIR`):

```bash
source /opt/prad2/bin/prad2_setup.sh    # bash / zsh
source /opt/prad2/bin/prad2_setup.csh   # csh / tcsh
```

`prad2_setup.csh` bakes the install prefix in at configure time, so
sourcing it from inside another script (e.g. a JLab farm wrapper that
also `module load`s ROOT) works correctly.  If the install tree is
later moved, `setenv PRAD2_DIR <new-prefix>` before sourcing to
override.

### Python tools (optional)

The GEM Python tools (`gem_event_viewer`, `gem_cluster_view`, …)
require `matplotlib`, `numpy`, and `PyQt6`.  On shared systems where
`pip install` into site-packages is not permitted, create a venv
next to the install and point the wrapper at it:

```bash
python3 -m venv /opt/prad2/venv
source /opt/prad2/venv/bin/activate     # or activate.csh for tcsh
pip install -r /path/to/prad2evviewer/scripts/requirements.txt
```

Source the venv's `activate` *before* sourcing `prad2_setup.csh` in
your personal wrapper script.

**RHEL 9 caveat.** CMake's `find_package(Python)` may pick up
`/usr/bin/python3.12` over the venv's 3.9.  Pass
`-DPython_EXECUTABLE=$(which python3)` at configure time to pin the
interpreter; otherwise `prad2py.cpython-312-*.so` will not be
importable from the 3.9 venv.

## Configuration

The server reads three top-level configuration files from `database/`:

- `monitor_config.json` (`-c`) — GUI/online server settings: waveform, hycal_hist, LMS, EPICS, livetime, beam status, physics display cuts, GEM diagnostics, elog, colour ranges, online/ET.
- `reconstruction_config.json` (`-r`) — `runinfo` pointer plus cluster/hit reconstruction knobs (HyCal clustering, per-detector GEM `ClusterConfig` with `default` + per-id overrides, matching σ parameters).
- `daq_config.json` (`-D`) — DAQ + raw decoding: event tags, bank tags, ROC layout, sync format, optional pulse-template pile-up deconvolution (`nnls_deconv` block; Levenberg–Marquardt fit), file pointers (`hycal_map_file`, `gem_map_file`, `pedestal_file`).

PRad support: use the mirrored `database/prad1/` set
(`prad_daq_config.json`, `prad_monitor_config.json`,
`prad_reconstruction_config.json`).

## Project structure

```
prad2dec/           libprad2dec — EVIO/ET reader; FADC250/SSP/ADC1881M/TDC/VTP/DSC2 decoders; EPICS store
prad2det/           libprad2det — HyCal/GEM clustering, reconstruction, PipelineBuilder
python/             pybind11 bindings (prad2py)
src/                Server, Qt GUI, data-source layer
resources/          Web frontend (HTML / CSS / JS)
database/           DAQ config, channel maps, calibration constants, runinfo
test/               Generic EVIO diagnostic CLI tools (test/dev/ = not installed)
gem/                GEM tracker: gem_dump binary + Python tools + reference notes
analysis/           ROOT-based replay and physics analysis (optional)
calibration/        HyCal calibration scan tools and cosmic-gain HV iteration
scripts/            HyCal / tagger / GEM-HyCal matching Python utilities (installed)
scripts/daq_tool/   Dev-only DAQ-config editors (not installed)
scripts/dev_tool/   Dev-only one-shot generators (not installed)
scripts/shell/      Operator shell scripts (not installed)
cmake/              Build helpers (PradHelpers, WebDeps, bin-wrapper template)
docs/               Technical notes, ROL references, API documentation
```

## Contributors

Chao Peng — Argonne National Laboratory\
Weizhi Xiong, Yuan Li, Mingyu Li — Shandong University
