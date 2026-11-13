# prad2dec

Static library for reading CODA EVIO data and decoding the PRad / PRad-II
front-end electronics.

## Components

| Header | Role |
|---|---|
| `EvChannel` | EVIO file reader — bank-tree scan, lazy per-product decoding, sequential and random-access modes. |
| `EtChannel` | `EvChannel` subclass for live ET-system reading (built when `-DWITH_ET=ON`). |
| `Fadc250Decoder` | FADC250 composite (online-summed) waveform decoding. |
| `Fadc250RawDecoder` | FADC250 raw-mode waveform decoding. |
| `Fadc250FwAnalyzer` | Bit-faithful firmware Mode 1/2/3 emulator (TET / NSB / NSA / NSAT / NPED / MAXPED). |
| `WaveAnalyzer` | Software pedestal (median + MAD), peak finding, integration, and per-module-type pulse-template Levenberg–Marquardt pile-up deconvolution. |
| `Adc1881mDecoder` | ADC1881M decoding (legacy PRad). |
| `SspDecoder` | SSP/MPD/APV fiber decoder (GEM readout). |
| `TdcDecoder` | V1190 multi-hit TDC decoder (tagger crate). |
| `VtpDecoder` | VTP trigger-bit / cluster-summary decoder. |
| `Dsc2Decoder` | DSC2 (0xE115) scaler bank decoder, used for live-time accounting. |
| `EpicsData`, `EpicsStore` | EPICS slow-control records — per-event POD plus a run-scoped accumulator with channel-id registry, value persistence, and event-indexed lookup. |
| `PulseTemplateStore` | Per-module-type pulse-shape templates loaded from JSON for the deconvolution backend. |

The waveform algorithms (soft analyzer, firmware emulator, deconvolution)
are documented in detail in
[`docs/technical_notes/waveform_analysis/wave_analysis.md`](../docs/technical_notes/waveform_analysis/wave_analysis.md).

## Usage

### Sequential read

```cpp
#include "EvChannel.h"

evc::EvChannel ch;
ch.SetConfig(cfg);                    // evc::DaqConfig
ch.OpenAuto("data.evio");             // random-access if the file supports
                                      // it, sequential otherwise; the
                                      // public API is identical in either
                                      // mode.

while (ch.Read() == evc::status::success) {
    if (!ch.Scan()) continue;
    if (ch.GetEventType() != evc::EventType::Physics) continue;

    for (int i = 0; i < ch.GetNEvents(); ++i) {
        ch.SelectEvent(i);                   // set sub-event index, clear cache
        const auto &info = ch.Info();        // TI / trigger metadata
        const auto &fadc = ch.Fadc();        // decoded on first call,
                                             // cached thereafter
        // also: ch.Gem(), ch.Tdc(), ch.Vtp(), ch.Dsc()
    }
}
ch.Close();
```

`Info()` / `Fadc()` / `Gem()` / `Tdc()` / `Vtp()` / `Dsc()` each decode on
the first call after `SelectEvent()` and return a cached reference on
subsequent calls for the same sub-event — requesting two of them on one
event costs no more than requesting one.

### Random access

`OpenRandomAccess()` mmaps the file and asks evio to build an internal
event-pointer table at open time; subsequent `ReadEventByIndex(i)` calls
jump to any event in O(1), in either direction, without close/reopen.

```cpp
evc::EvChannel ch;
ch.SetConfig(cfg);
ch.OpenRandomAccess("data.evio");

int n = ch.GetRandomAccessEventCount();      // total evio events
for (int i : {0, n / 2, n - 1}) {
    if (ch.ReadEventByIndex(i) != evc::status::success) continue;
    if (!ch.Scan()) continue;
    ch.SelectEvent(0);
    auto &fadc = ch.Fadc();
    // ...
}
```

Random-access events are *evio* events (blocks); a CODA built-trigger
block can hold several physics sub-events, so after `ReadEventByIndex` +
`Scan` use `GetNEvents()` + `SelectEvent(i)` as in the sequential path.
(The decoders do not currently use the selected index — every sub-event of
a block decodes the same banks; `SelectEvent` only clears the cache when
the index changes.)
See `src/evio_data_source.cpp` for the two-pass pattern the server
uses (one Scan-only indexing pass to record `{evio_event, sub_event}`
pairs; subsequent random access via the index).

### EPICS accumulator

```cpp
#include "EpicsStore.h"

epics::EpicsStore store;
while (ch.Read() == evc::status::success) {
    ch.Scan();
    if (ch.GetEventType() == evc::EventType::Epics) {
        const auto &rec = ch.Epics();        // EpicsRecord (per-event POD)
        store.Feed(rec.event_number_at_arrival,
                   /*ti timestamp*/ 0, ch.ExtractEpicsText());
    }
}
float current;
store.GetValue(/*event_number=*/ 12345, "hallb_IPM2C21A_CUR", current);
```

### Legacy API

`DecodeEvent(i, event, ssp=nullptr, vtp=nullptr, tdc=nullptr)` is retained
as a thin compatibility wrapper that writes directly into caller-owned
structs without touching the lazy cache. New code should prefer
`SelectEvent()` together with the typed accessors.

## Dependencies

- [evio](https://github.com/JeffersonLab/evio) — required (the evio-6 or
  evio-4 C API; the fetch fallback uses `evio-6.0`).
- [et](https://github.com/JeffersonLab/et) — optional, for `EtChannel`.
- nlohmann/json — required; the `nlohmann_json::nlohmann_json` target is
  defined by the top-level build (`cmake/WebDeps.cmake`), so `prad2dec`
  does not configure on its own.

Both resolve from the Hall-B CODA installation by default; if not found,
CMake fetches from GitHub. Override with `-DEVIO_SOURCE=fetch` or
`-DET_SOURCE=fetch`.
