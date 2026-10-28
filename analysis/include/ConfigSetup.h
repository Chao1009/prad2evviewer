#pragma once
//=============================================================================
// ConfigSetup.h — analysis-side helpers around RunInfo configs.
//
// The RunConfig struct, LoadRunConfig() and WriteRunConfig() live in
// prad2det/include/RunInfoConfig.h so they can be reused by the viewer,
// Python bindings and ROOT scripts. This header adds:
//   - the analysis-only gRunConfig global + backward-compat aliases
//   - BuildLabTransforms() — bridge from RunConfig to DetectorTransform
//   - ApplyToLab() / ApplyToLocal() / ApplyToHyCal() per-hit helpers
//   - get_run_str() / get_run_int() filename parsers
//
// "Detector-frame → lab-frame" is owned by prad2det's DetectorTransform —
// build the per-detector poses once via BuildLabTransforms(), then call
// xform.toLab(x, y[, z]) per hit.
//=============================================================================

#include "DetectorTransform.h"
#include "EvioFiles.h"
#include "RunInfoConfig.h"

#include <array>
#include <iostream>
#include <string>

namespace analysis {

using RunConfig = ::prad2::RunConfig;
using ::prad2::LoadRunConfig;
using ::prad2::WriteRunConfig;

// Global geometry config for single-run tools.
// Multi-run code should capture LoadRunConfig()'s return value into a
// local RunConfig instead of relying on this global.
inline RunConfig gRunConfig;

// HyCal + 4-GEM DetectorTransform bundle, ready for per-hit `toLab` calls.
// Build once per run (rotation matrices precomputed inside) and reuse.
struct LabTransforms {
    DetectorTransform hycal;
    std::array<DetectorTransform, 4> gem;
};

inline LabTransforms BuildLabTransforms(const RunConfig &geo = gRunConfig)
{
    LabTransforms t;
    t.hycal.set(geo.hycal_x, geo.hycal_y, geo.hycal_z,
                geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        t.gem[d].set(geo.gem_x[d], geo.gem_y[d], geo.gem_z[d],
                     geo.gem_tilt_x[d], geo.gem_tilt_y[d], geo.gem_tilt_z[d]);
    }
    return t;
}

// Apply a DetectorTransform to a hit in-place (lab = R*[x,y,z] + [tx,ty,tz]).
template <typename Hit>
inline void ApplyToLab(const DetectorTransform &xform, Hit &h)
{
    float lx, ly, lz;
    xform.toLab(h.x, h.y, h.z, lx, ly, lz);
    h.x = lx; h.y = ly; h.z = lz;
}
template <typename Hit>
inline void ApplyToLocal(const DetectorTransform &xform, Hit &h)
{
    float dx, dy, dz;
    xform.labToLocal(h.x, h.y, h.z, dx, dy, dz);
    h.x = dx; h.y = dy; h.z = dz;
}
// Transform a hit position to the HyCal coordinate system: a shift by
// (target_x, target_y).  Project the hit to the HyCal surface first.
template <typename Hit>
inline void ApplyToHyCal(Hit &h, const RunConfig &geo = gRunConfig)
{
    float dx = geo.target_x, dy = geo.target_y;
    h.x += dx;
    h.y += dy;
}

// --- run number utilities ---------------------------------------------------
// Run number embedded in a file name, as parsed by prad2::run_number_from_path
// (".../prad_<digits>..." or "run_<digits>").  Returns "unknown" / -1 on failure.
inline int get_run_int(const std::string &file_name)
{
    const int run = prad2::run_number_from_path(file_name);
    if (run < 0)
        std::cerr << "Warning: cannot extract run number from file name " << file_name << ", using 'unknown'.\n";
    return run;
}

inline std::string get_run_str(const std::string &file_name)
{
    const int run = get_run_int(file_name);
    return run < 0 ? "unknown" : std::to_string(run);
}

} // namespace analysis
