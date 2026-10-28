#pragma once
//=============================================================================
// HyCalTimeCalib.h — per-module raw-time calibration offsets for HyCal.
//
// Applies a time-offset table to HyCal modules in place by setting
// Module::time_offset, so downstream code with a module pointer uses
//   calib_time = raw_time - mod->time_offset
// PipelineBuilder loads the file named by reconstruction_config.json
// "hycal": {"time_calib_file": ...}; standalone tools call
// LoadHyCalTimeCalib(path, hycal) after hycal.Init().
//
// File format (under database/, e.g. hycal_time_offsets/25308.json; row
// lookup and "default" fallback per HyCalModuleTable.h):
//
//     {
//       "default": 0.0,
//       "modules": [
//         {"name": "W735", "offset_ns": -0.74},
//         {"name": "W736", "offset_ns":  1.21}
//       ]
//     }
//=============================================================================

#include "HyCalModuleTable.h"

#include <iostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalTimeCalibSummary {
    float default_off = 0.f;
    int   n_overrides = 0;
    int   n_unknown   = 0;             // rows naming no known module
};

// Sets every module's time_offset (to `def_off` when `path` is empty or
// unreadable).  Warnings go to *log (if non-null).
inline HyCalTimeCalibSummary LoadHyCalTimeCalib(const std::string &path,
                                                fdec::HyCalSystem &hycal,
                                                float def_off = 0.f,
                                                std::ostream *log = &std::cerr)
{
    HyCalTimeCalibSummary summary;
    summary.default_off = def_off;
    std::vector<float> off;
    const auto stats = LoadHyCalModuleTable(path, hycal, "offset_ns",
        ParseHyCalTableNumber, summary.default_off, off, log);
    for (int i = 0; i < hycal.module_count(); ++i)
        hycal.module(i).time_offset = off[i];
    summary.n_overrides = stats.n_overrides;
    summary.n_unknown   = stats.n_unknown;
    return summary;
}

} // namespace prad2
