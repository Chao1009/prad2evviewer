#pragma once
//=============================================================================
// GemEventJson.h — JSON view of one processed GEM event
//
// Serializes a GemSystem after ProcessEvent() (+ Reconstruct() for the
// detector block) into the "zs_apvs" / "detectors" schema that
// `gem_dump -m evdump` writes and gem_event_viewer draws.
//
// round = true rounds ADC values / charges to 0.1 and positions to 0.01 mm
// (file dumps); false keeps full float precision (live view).
//=============================================================================

#include <nlohmann/json.hpp>

namespace ssp { struct SspEventData; }

namespace gem
{

class GemSystem;

// One entry per APV with at least one channel that survived zero
// suppression, in GemSystem index order, or in the readout order of `evt`
// when given (the event ProcessEvent() ran on):
//   {"crate", "mpd", "adc",
//    "channels": {"<ch>": {"charge", "max_timebin", "cross_talk", "ts_adc"[6]}}}
// charge / max_timebin are the first maximum over the processed samples.
nlohmann::json ZsApvsToJson(const GemSystem &sys, bool round = false,
                            const ssp::SspEventData *evt = nullptr);

// One entry per detector:
//   {"id", "name", "x_pitch", "y_pitch", "x_strips", "y_strips",
//    "x_clusters" / "y_clusters": [{"position", "peak_charge", "total_charge",
//        "max_timebin", "cross_talk", "size", "hit_strips"}],
//    "hits_2d": [{"x", "y", "x_charge", "y_charge", "x_peak", "y_peak",
//        "x_size", "y_size"}]}
nlohmann::json DetectorsToJson(const GemSystem &sys, bool round = false);

} // namespace gem
