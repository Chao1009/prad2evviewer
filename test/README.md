# Test & Diagnostic Tools

## evio_dump

EVIO file structure inspector.

```bash
evio_dump <file> [-m mode] [-n N] [-D daq_config.json]
```

Modes: (default) summary by tag, `tree`, `tags`, `epics`, `event` (detail for event N), `triggers`.

## livetime

DAQ live time calculator. Two methods: DSC2 gated/ungated scalers from SYNC events, and pulser event counting. Accepts a single file, a base name (auto-discovers `.00000`, `.00001`, ...), or a directory.

```bash
livetime <input> [-D daq_config.json] [-f freq_hz] [-t interval_sec]
```

| Option | Description |
|--------|-------------|
| `-f <hz>` | Pulser frequency (default: 100 Hz) |
| `-t <sec>` | Periodic report interval (default: 10 sec) |

Outputs a periodic table (DSC2 and pulser live time vs elapsed time) followed by a summary with per-channel DSC2 breakdowns.

## ts_dump

Dump TI timestamp and trigger information for each event.

```bash
ts_dump <file> [-n max_events] [-D daq_config.json]
```

## ped_calc

Compute HyCal per-channel pedestals from trigger-selected events.

```bash
ped_calc <evio_file> -D <daq_config.json> [-t bit] [-o file.json] [-n N]
```

Default trigger bit 3 (LMS_Alpha/pedestal for PRad).

## gem_dump

GEM data diagnostic tool. Decodes SSP/MPD banks, runs GEM reconstruction, prints diagnostics or dumps events to JSON for visualization.

```bash
gem_dump <evio_file> [options]
```

| Option | Description |
|--------|-------------|
| `-m <mode>` | `summary` (default), `raw`, `hits`, `clusters`, `ped`, `evdump` |
| `-D <file>` | DAQ config (auto-searches `database/daq_config.json`) |
| `-G <file>` | GEM map (default: `gem_map.json` next to DAQ config) |
| `-P <file>` | Pedestal file (required for `hits`/`clusters`/`evdump`) |
| `-n <N>` | Max events (default: 10, 0=all) |
| `-e <N>` | Single event N (1-based) |
| `-t <bit>` | Trigger bit filter (-1=all) |
| `-o <file>` | Output file (`ped`/`evdump` modes) |
| `-z <sigma>` | Override zero-suppression threshold |
| `-f <filter>` | APV filter for `evdump` (see below) |

The `evdump` mode outputs JSON with three pipeline layers: `raw_apvs`, `zs_apvs`, and per-detector `clusters` + `hits_2d`. Use with `gem/gem_cluster_view.py`.

Filter syntax: `-f field=val[,val]:min_dets` (fields: `pos`, `plane`, `match`, `orient`, `det`).

## Dev-only tools (not installed) — sources under `test/dev/`

**ts_dump** — dump trigger timestamps from an EVIO file.

**livetime** — compute DAQ livetime from pulser events.

### ET-specific (requires `-DWITH_ET=ON`)

**evc_scan** — Smoke-test: read EVIO buffers, scan events, or connect to ET.

**et_feeder** — Replay EVIO into ET at configurable rate: `et_feeder data.evio -f /tmp/et -i 50`

**evet_diff** — Compare ET-streamed vs disk-read events word-by-word (pairs with `et_feeder`).
