#pragma once
//=============================================================================
// HyCalTimeCalib.h — per-module raw-time calibration offsets for HyCal.
//
// Built once per run by PipelineBuilder (or directly by analysis tools via
// LoadHyCalTimeCalib).  Applies a time-offset table to HyCal modules in-place
// by setting Module::time_offset, so downstream code can use
//   calib_time = raw_time - mod.time_offset
// consistently.
//
// Typical usage via PipelineBuilder:
//
//   auto p = prad2::PipelineBuilder()
//       .set_recon_config("reconstruction_config.json")
//       .build();
//
//   // If reconstruction_config.json contains:
//   //   "hycal": { "time_calib_file": "hycal_time_offsets/foo.json" }
//   // then p.hycal.module(i).time_offset is already populated here.
//
// Direct usage in a standalone tool:
//
//   fdec::HyCalSystem hycal;
//   hycal.Init("database/hycal_map.json");
//   prad2::LoadHyCalTimeCalib("database/hycal_time_offsets/foo.json", hycal);
//
// Referencing a module's offset during event processing:
//
//   const auto *mod = hycal.module_by_id(ev.module_id[j]);
//   if (!mod) return;
//   const double raw_time = ev.peak_time[j][0];
//   const double calib_time = raw_time - mod->time_offset;
//
// Since the offsets are stored directly on fdec::Module, callers can access
// mod->time_offset anywhere they already have a module pointer/reference.
//
// File format (under database/, e.g. hycal_time_offsets/25308.json):
//
//     {
//       "default": 0.0,
//       "modules": [
//         {"name": "W735", "offset_ns": -0.74},
//         {"name": "W736", "offset_ns":  1.21}
//       ]
//     }
//
// Resolution per module:
//   * Module listed in `modules`  -> per-module offset.
//   * Else if `default` in file   -> file's default offset.
//   * Else                        -> `def_off` passed by the caller
//                                    (typically 0.0).
//
// Keyed primarily by module name to match other HyCal calibration tables.
//=============================================================================

#include "HyCalSystem.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

namespace prad2 {

struct HyCalTimeCalibSummary {
	float default_off = 0.f;
	int   n_overrides = 0;
};

inline HyCalTimeCalibSummary LoadHyCalTimeCalib(const std::string &path,
												fdec::HyCalSystem &hycal,
												float def_off = 0.f)
{
	HyCalTimeCalibSummary summary;
	summary.default_off = def_off;

	const int n = hycal.module_count();
	for (int i = 0; i < n; ++i) {
		hycal.module(i).time_offset = def_off;
	}

	if (path.empty()) return summary;

	std::ifstream f(path);
	if (!f) {
		std::cerr << "Warning: cannot open HyCal time-calib file " << path
				  << ", using uniform " << def_off << " ns.\n";
		return summary;
	}

	auto j = nlohmann::json::parse(f, nullptr, false, true);
	if (j.is_discarded()) {
		std::cerr << "Warning: failed to parse " << path
				  << ", using uniform " << def_off << " ns.\n";
		return summary;
	}

	if (j.contains("default") && j["default"].is_number()) {
		const float file_default = j["default"].get<float>();
		summary.default_off = file_default;
		for (int i = 0; i < n; ++i) {
			hycal.module(i).time_offset = file_default;
		}
	}

	int unknown = 0;
	if (j.contains("modules") && j["modules"].is_array()) {
		for (const auto &m : j["modules"]) {
			if (!m.contains("offset_ns") || !m["offset_ns"].is_number()) continue;

			const fdec::Module *mod = nullptr;
			if (m.contains("name") && m["name"].is_string()) {
				mod = hycal.module_by_name(m["name"].get<std::string>());
			}
			if (!mod && m.contains("module_id") && m["module_id"].is_number_integer()) {
				mod = hycal.module_by_id(m["module_id"].get<int>());
			}
			if (!mod) {
				++unknown;
				continue;
			}

			hycal.module(mod->index).time_offset = m["offset_ns"].get<float>();
			++summary.n_overrides;
		}
	}

	std::cerr << "HyCal time calib: " << summary.n_overrides << " module overrides";
	if (unknown) std::cerr << " (" << unknown << " unknown modules skipped)";
	std::cerr << ", default=" << summary.default_off << " ns from " << path
			  << "\n";
	return summary;
}

} // namespace prad2
