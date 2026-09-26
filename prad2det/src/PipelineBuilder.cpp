#include "PipelineBuilder.h"

#include "EvioFiles.h"
#include "InstallPaths.h"
#include "JsonUtil.h"
#include "SspData.h"
#include "load_daq_config.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <utility>

using nlohmann::json;

namespace prad2 {

namespace {

bool read_json_bool(const json &obj, const char *key, bool def)
{
    if (!obj.contains(key)) return def;
    const auto &v = obj[key];
    if (v.is_boolean()) return v.get<bool>();
    if (v.is_number_integer()) return v.get<int>() != 0;
    if (v.is_string()) {
        std::string s = v.get<std::string>();
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (s == "true" || s == "1" || s == "yes" || s == "on") return true;
        if (s == "false" || s == "0" || s == "no" || s == "off") return false;
    }
    return def;
}

void apply_gem_cluster_overrides(const json &j, gem::ClusterConfig &cfg)
{
    if (j.contains("min_cluster_hits"))    cfg.min_cluster_hits    = j["min_cluster_hits"];
    if (j.contains("max_cluster_hits"))    cfg.max_cluster_hits    = j["max_cluster_hits"];
    if (j.contains("consecutive_thres"))   cfg.consecutive_thres   = j["consecutive_thres"];
    if (j.contains("split_thres"))         cfg.split_thres         = j["split_thres"];
    if (j.contains("cross_talk_width"))    cfg.cross_talk_width    = j["cross_talk_width"];
    if (j.contains("cross_talk_peak_ratio_max"))
        cfg.cross_talk_peak_ratio_max = j["cross_talk_peak_ratio_max"];
    if (j.contains("charac_dists") && j["charac_dists"].is_array()) {
        cfg.charac_dists.clear();
        for (auto &v : j["charac_dists"])
            cfg.charac_dists.push_back(v.get<float>());
    }
    if (j.contains("match_mode"))          cfg.match_mode          = j["match_mode"];
    if (j.contains("match_adc_asymmetry")) cfg.match_adc_asymmetry = j["match_adc_asymmetry"];
    if (j.contains("match_time_diff"))     cfg.match_time_diff     = j["match_time_diff"];
    if (j.contains("match_ts_period"))     cfg.ts_period           = j["match_ts_period"];

    // SBS-style (mpd_gem_view_ssp Cuts) quality cuts — all off by default.
    // "strip_mean_time_range": [lo, hi] ns; [] or null disables (±inf).
    // Any other shape is ignored (keeps the current value, no throw).
    if (j.contains("strip_mean_time_range")) {
        const auto &r = j["strip_mean_time_range"];
        if (r.is_null() || (r.is_array() && r.empty())) {
            cfg.strip_time_min = -std::numeric_limits<float>::infinity();
            cfg.strip_time_max =  std::numeric_limits<float>::infinity();
        } else if (r.is_array() && r.size() == 2 &&
                   r[0].is_number() && r[1].is_number()) {
            cfg.strip_time_min = r[0].get<float>();
            cfg.strip_time_max = r[1].get<float>();
        }
    }
    if (j.contains("strip_unimodal_shape"))
        cfg.strip_unimodal = read_json_bool(j, "strip_unimodal_shape", cfg.strip_unimodal);
    if (j.contains("seed_min_peak_adc"))    cfg.seed_min_peak_adc    = j["seed_min_peak_adc"];
    if (j.contains("seed_min_sum_adc"))     cfg.seed_min_sum_adc     = j["seed_min_sum_adc"];
    if (j.contains("strip_time_agreement")) cfg.strip_time_agreement = j["strip_time_agreement"];
    if (j.contains("strip_ts_corr_min"))    cfg.strip_ts_corr_min    = j["strip_ts_corr_min"];
}

void apply_hycal_cluster_overrides(const json &j, fdec::ClusterConfig &cfg)
{
    if (j.contains("min_module_energy"))  cfg.min_module_energy  = j["min_module_energy"];
    if (j.contains("min_center_energy"))  cfg.min_center_energy  = j["min_center_energy"];
    if (j.contains("min_cluster_energy")) cfg.min_cluster_energy = j["min_cluster_energy"];
    if (j.contains("min_cluster_size"))   cfg.min_cluster_size   = j["min_cluster_size"];
    if (j.contains("corner_conn"))        cfg.corner_conn        = j["corner_conn"];
    if (j.contains("split_iter"))         cfg.split_iter         = j["split_iter"];
    if (j.contains("least_split"))        cfg.least_split        = j["least_split"];
    if (j.contains("log_weight_thres"))   cfg.log_weight_thres   = j["log_weight_thres"];
    if (j.contains("seed_time_window"))   cfg.seed_time_window   = j["seed_time_window"];
    if (j.contains("non_linear_corr"))    cfg.non_linear_corr    = j["non_linear_corr"];
    if (j.contains("energy_bias_correction"))
        cfg.energy_bias_correction = read_json_bool(
            j, "energy_bias_correction", cfg.energy_bias_correction);
    if (j.contains("leakage_correction")) cfg.leakage_correction = j["leakage_correction"];
    if (j.contains("leakage_iterations")) cfg.leakage_iterations = j["leakage_iterations"];
    if (j.contains("least_leakage_fraction"))
        cfg.least_leakage_fraction = j["least_leakage_fraction"];
    if (j.contains("max_leakage_fraction"))
        cfg.max_leakage_fraction = j["max_leakage_fraction"];
    if (j.contains("leakage_convergence_rel"))
        cfg.leakage_convergence_rel = j["leakage_convergence_rel"];
}

} // namespace

PipelineBuilder &PipelineBuilder::set_database_dir(std::string p)      { database_dir_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_daq_config(std::string p)        { daq_config_path_      = std::move(p); return *this; }

PipelineBuilder &PipelineBuilder::set_loaded_daq_config(evc::DaqConfig cfg)
{
    daq_config_loaded_      = std::move(cfg);
    have_loaded_daq_config_ = true;
    return *this;
}
PipelineBuilder &PipelineBuilder::set_recon_config(std::string p)      { recon_config_path_    = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_runinfo(std::string p)           { runinfo_path_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_map(std::string p)         { hycal_map_path_       = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_map(std::string p)           { gem_map_path_         = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_calib(std::string p)       { hycal_calib_path_     = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_time_calib(std::string p)  { hycal_time_calib_path_ = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_time_cut(std::string p)    { hycal_time_cut_path_  = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_hycal_rf_offset(std::string p)   { hycal_rf_offset_path_ = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_pedestal(std::string p)      { gem_pedestal_path_    = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_gem_common_mode(std::string p)   { gem_common_mode_path_ = std::move(p); return *this; }
PipelineBuilder &PipelineBuilder::set_run_number(int n)                { run_number_ = n; return *this; }

PipelineBuilder &PipelineBuilder::set_run_number_from_evio(const std::string &p)
{
    int n = prad2::run_number_from_path(p);
    if (n > 0) run_number_ = n;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_log_stream(std::ostream *s)
{
    log_stream_         = s;
    log_stream_default_ = false;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_log_pedestal_checksum(bool b)
{
    log_pedestal_checksum_ = b;
    return *this;
}

PipelineBuilder &PipelineBuilder::set_path_resolver(
    std::function<std::string(const std::string &)> resolver)
{
    path_resolver_ = std::move(resolver);
    return *this;
}

Pipeline PipelineBuilder::build()
{
    // --- resolve effective log target -------------------------------------
    std::ostream *log = log_stream_default_ ? &std::cerr : log_stream_;
    // One line per call; a fresh stream keeps manipulators from leaking.
    auto LOG = [log](const auto &...parts) {
        if (!log) return;
        std::ostringstream oss;
        (oss << ... << parts);
        (*log) << oss.str() << '\n';
    };

    // --- resolve database dir + path resolver -----------------------------
    const std::string db_dir = database_dir_.empty() ? prad2::database_dir() : database_dir_;
    auto resolve = [&](const std::string &p) -> std::string {
        if (p.empty() || prad2::is_absolute_path(p)) return p;
        if (path_resolver_) return path_resolver_(p);
        return prad2::resolve_db_path(p, db_dir);
    };
    // The builder override if set, else the default.
    auto pick = [&](const std::string &over, const std::string &def) {
        return resolve(over.empty() ? def : over);
    };

    Pipeline out;

    // --- 1. DAQ config ----------------------------------------------------
    if (have_loaded_daq_config_) {
        out.daq_cfg = std::move(daq_config_loaded_);
        // If the caller knows the source path, surface it; otherwise leave
        // empty so logs make it obvious nothing was loaded by the builder.
        out.daq_config_path = daq_config_path_;
        LOG("[setup] DAQ config : (caller-supplied)",
            daq_config_path_.empty() ? "" : " " + daq_config_path_);
    } else {
        const std::string daq_path = pick(daq_config_path_, "daq_config.json");
        if (daq_path.empty() || !evc::load_daq_config(daq_path, out.daq_cfg)) {
            throw std::runtime_error(
                "PipelineBuilder: cannot load DAQ config '" + daq_path + "'");
        }
        out.daq_config_path = daq_path;
        LOG("[setup] DAQ config : ", daq_path);
    }

    // --- 2. GEM crate remap from daq_cfg.roc_tags ------------------------
    for (const auto &re : out.daq_cfg.roc_tags)
        if (re.type == "gem")
            out.gem_crate_remap[(int)re.tag] = re.crate;

    // --- 3. recon config (soft — falls back to library defaults) ----------
    const std::string recon_path = pick(recon_config_path_, "reconstruction_config.json");
    json recon = json::object();
    std::string recon_err;
    prad2::read_json_file(recon_path, recon, &recon_err);
    if (recon.empty()) {
        LOG("[WARN] PipelineBuilder: cannot load reconstruction_config '",
            recon_path, "'", recon_err.empty() ? "" : " (" + recon_err + ")",
            " — proceeding with library defaults.");
    } else {
        out.recon_config_path = recon_path;
    }

    // --- 4. runinfo (soft — defaults to RunConfig{} if absent) ------------
    std::string recon_runinfo;
    if (recon.contains("runinfo") && recon["runinfo"].is_string())
        recon_runinfo = recon["runinfo"].get<std::string>();
    const std::string ri_path = pick(runinfo_path_, recon_runinfo);

    out.run_number = run_number_;
    if (!ri_path.empty()) {
        if (run_number_ > 0)
            LOG("[setup] Run number : ", run_number_);
        else
            LOG("[setup] Run number : (latest entry — no run-number override)");
        out.run_cfg      = prad2::LoadRunConfig(ri_path, run_number_);
        out.runinfo_path = ri_path;
        LOG("[setup] RunInfo    : ", ri_path,
            "  beam=", static_cast<int>(out.run_cfg.Ebeam),
            " MeV  hycal_z=", std::fixed, std::setprecision(1),
            out.run_cfg.hycal_z, " mm");
    } else {
        LOG("[WARN] PipelineBuilder: no runinfo path resolved — using default "
            "RunConfig (geometry/calibration paths empty).");
    }

    // --- 5. detector transforms ------------------------------------------
    out.hycal_transform.set(
        out.run_cfg.hycal_x, out.run_cfg.hycal_y, out.run_cfg.hycal_z,
        out.run_cfg.hycal_tilt_x, out.run_cfg.hycal_tilt_y, out.run_cfg.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        out.gem_transforms[d].set(
            out.run_cfg.gem_x[d], out.run_cfg.gem_y[d], out.run_cfg.gem_z[d],
            out.run_cfg.gem_tilt_x[d], out.run_cfg.gem_tilt_y[d], out.run_cfg.gem_tilt_z[d]);
    }

    // --- 6. HyCal map (soft — leaves hycal default-constructed on failure) -
    const std::string hc_map = pick(hycal_map_path_, "hycal_map.json");
    if (hc_map.empty()) {
        LOG("[WARN] PipelineBuilder: no HyCal map resolved.");
    } else if (!out.hycal.Init(hc_map)) {
        LOG("[WARN] PipelineBuilder: HyCalSystem.Init failed for '", hc_map, "'.");
    } else {
        out.hycal_map_path = hc_map;
    }

    if (recon.contains("hycal") && recon["hycal"].is_object())
        prad2::read_json_array(recon["hycal"], "energy_resolution",
            out.hycal_energy_res[0], out.hycal_energy_res[1], out.hycal_energy_res[2]);
    out.hycal.SetEnergyResolutionParams(
        out.hycal_energy_res[0], out.hycal_energy_res[1], out.hycal_energy_res[2]);
    LOG("[setup] HC sigma_E : E*sqrt(",
        std::fixed, std::setprecision(3), out.hycal_energy_res[0],
        "^2/E_GeV+", out.hycal_energy_res[2],
        "^2+", out.hycal_energy_res[1], "^2/E_GeV^2)/100");

    // --- 6a. Shared cluster profile --------------------------------------
    const std::string pwo_profile = resolve("cluster_profiles/prof_pwo.dat");
    const std::string glass_profile = resolve("cluster_profiles/prof_lg.dat");
    auto table_profile = std::make_shared<fdec::Geant4Profile>();
    if (table_profile->Load(fdec::ModuleType::PbWO4, pwo_profile) &&
        table_profile->Load(fdec::ModuleType::PbGlass, glass_profile)) {
        out.hycal_profile = std::move(table_profile);
        LOG("[setup] HC profile : loaded PWO and PbGlass tables");
    } else {
        out.hycal_profile = std::make_shared<fdec::SimpleProfile>();
        LOG("[WARN] HC profile : failed to load PWO/LG tables; using SimpleProfile");
    }

    // --- 6b. Per-run HyCal dead modules -----------------------------------
    {
        const auto dead = prad2::ApplyHyCalDeadModules(
            out.run_cfg.hycal_dead_modules, out.hycal);
        std::ostringstream oss;
        oss << "[setup] HC dead    : modules=" << dead.n_dead
            << "  neighbors=" << dead.n_dead_neighbors;
        if (dead.n_unknown > 0)
            oss << "  unknown=" << dead.n_unknown;
        LOG(oss.str());
    }

    // --- 7. HyCal calibration --------------------------------------------
    const std::string hc_calib = pick(hycal_calib_path_, out.run_cfg.energy_calib_file);
    if (!hc_calib.empty()) {
        int n = out.hycal.LoadCalibration(hc_calib);
        out.hycal_calib_path = hc_calib;
        LOG("[setup] HC calib   : ", hc_calib, " (", n, " modules)");
    } else {
        LOG("[WARN] no HyCal calibration file — energies will be wrong.");
    }

    // --- 7b. HyCal per-module raw-time calibration table -----------------
    // Applies offsets directly into HyCalSystem::Module::time_offset so
    // downstream code can use calib_time = raw_time - mod.time_offset.
    // This path is reconstruction-wide, so it comes from recon config
    // (or an explicit builder override), not runinfo.
    {
        const bool has_hycal_cfg = recon.contains("hycal") && recon["hycal"].is_object();
        const bool enable_time_calib = has_hycal_cfg
            ? read_json_bool(recon["hycal"], "time_calib", true)
            : true;

        if (!enable_time_calib) {
            (void)prad2::LoadHyCalTimeCalib("", out.hycal, 0.f, log);
            out.hycal_time_calib_path.clear();
            LOG("[setup] HC t calib: disabled by recon config (hycal.time_calib=false), using 0 ns offsets.");
        } else {
            std::string recon_file;
            if (has_hycal_cfg && recon["hycal"].contains("time_calib_file")
                    && recon["hycal"]["time_calib_file"].is_string())
                recon_file = recon["hycal"]["time_calib_file"].get<std::string>();
            const std::string hc_time_calib_path = pick(hycal_time_calib_path_, recon_file);
            const auto time_calib = prad2::LoadHyCalTimeCalib(
                hc_time_calib_path, out.hycal, 0.f, log);
            out.hycal_time_calib_path = hc_time_calib_path;

            std::ostringstream oss;
            oss << "[setup] HC t calib: default=" << time_calib.default_off << " ns";
            if (time_calib.n_overrides > 0) {
                oss << "  per-module=" << time_calib.n_overrides
                    << " (" << hc_time_calib_path << ")";
            }
            if (time_calib.n_unknown > 0)
                oss << "  unknown=" << time_calib.n_unknown;
            LOG(oss.str());
        }
    }

    // --- 7c. HyCal per-module time-cut table -----------------------------
    // Always populate `out.hycal_time_cuts` (uniform default when no file)
    // so per-event callers can use a single code path.  The path comes
    // from runinfo's `time_cuts.hycal_module_file` unless overridden.
    {
        const std::string hc_time_path = pick(hycal_time_cut_path_, out.run_cfg.hycal_time_cut_file);
        out.hycal_time_cuts = prad2::LoadHyCalTimeCuts(
            hc_time_path, out.hycal,
            out.run_cfg.hc_time_win_lo, out.run_cfg.hc_time_win_hi, log);
        out.hycal_time_cut_path = hc_time_path;

        std::ostringstream oss;
        oss << "[setup] HC time   : default=["
            << out.hycal_time_cuts.default_lo << ", "
            << out.hycal_time_cuts.default_hi << "] ns";
        if (out.hycal_time_cuts.n_overrides > 0) {
            oss << "  per-module=" << out.hycal_time_cuts.n_overrides
                << " (" << hc_time_path << ")";
        }
        if (out.hycal_time_cuts.n_unknown > 0)
            oss << "  unknown=" << out.hycal_time_cuts.n_unknown;
        LOG(oss.str());
    }

    // --- 7d. HyCal per-module HyCal→RF offset table ----------------------
    // Always populate `out.hycal_rf_offsets` (uniform 0 ns when no file)
    // so the per-event Δt fill uses a single call.  Path comes from
    // runinfo's `time_cuts.hycal_rf_offsets` unless overridden.
    {
        const std::string rf_off_path = pick(hycal_rf_offset_path_, out.run_cfg.hycal_rf_offset_file);
        out.hycal_rf_offsets = prad2::LoadHyCalRfOffsets(
            rf_off_path, out.hycal, 0.f, log);
        out.hycal_rf_offset_path = rf_off_path;

        std::ostringstream oss;
        oss << "[setup] HC RF off : default=" << out.hycal_rf_offsets.default_off
            << " ns";
        if (out.hycal_rf_offsets.n_overrides > 0) {
            oss << "  per-module=" << out.hycal_rf_offsets.n_overrides
                << " (" << rf_off_path << ")";
        }
        if (out.hycal_rf_offsets.n_unknown > 0)
            oss << "  unknown=" << out.hycal_rf_offsets.n_unknown;
        LOG(oss.str());
    }

    // --- 8. matching (HyCal sigma + GEM sigma + target sigma) ------------
    if (recon.contains("matching")) {
        const auto &m = recon["matching"];
        prad2::read_json_array(m, "hycal_pos_res",
            out.hycal_pos_res[0], out.hycal_pos_res[1], out.hycal_pos_res[2]);
        if (m.contains("gem_pos_res") && m["gem_pos_res"].is_array()) {
            out.gem_pos_res.clear();
            for (auto &v : m["gem_pos_res"])
                out.gem_pos_res.push_back(v.get<float>());
        }
        prad2::read_json_array(m, "target_pos_res",
            out.target_pos_res[0], out.target_pos_res[1], out.target_pos_res[2]);
        if (m.contains("match_method") && m["match_method"].is_number_integer()) {
            out.match_method = m["match_method"].get<int>();
        }
    }
    out.hycal.SetPositionResolutionParams(
        out.hycal_pos_res[0], out.hycal_pos_res[1], out.hycal_pos_res[2]);
    LOG("[setup] HC sigma(E)= sqrt((",
        std::fixed, std::setprecision(3), out.hycal_pos_res[0],
        "/sqrt(E_GeV))^2+(", out.hycal_pos_res[1],
        "/E_GeV)^2+", out.hycal_pos_res[2], "^2) mm");
    {
        std::ostringstream oss;
        oss << "[setup] GEM sigma  : [";
        for (size_t i = 0; i < out.gem_pos_res.size(); ++i) {
            if (i) oss << ",";
            oss << std::fixed << std::setprecision(3) << out.gem_pos_res[i];
        }
        oss << "] mm";
        LOG(oss.str());
    }

    // --- 9. HyCal cluster config -----------------------------------------
    if (recon.contains("hycal") && recon["hycal"].is_object())
        apply_hycal_cluster_overrides(recon["hycal"], out.hycal_cluster_cfg);
    out.hycal_cluster_cfg.profile = out.hycal_profile;

    if (out.hycal_cluster_cfg.energy_bias_correction) {
        if (out.run_cfg.Ebeam > 0.f) {
            const auto set = fdec::SelectHyCalEnergyBiasSet(out.run_cfg.Ebeam);
            const std::string base = std::string("energy_bias/") + set.file_prefix;
            out.hycal_energy_bias_nominal = set.nominal_mev;
            out.hycal_energy_bias_ee_path = resolve(base + "_ee.json");
            out.hycal_energy_bias_ep_path = resolve(base + "_ep.json");
            out.hycal_energy_bias = fdec::LoadHyCalEnergyBias(
                out.hycal_energy_bias_ee_path,
                out.hycal_energy_bias_ep_path,
                out.hycal, out.run_cfg.Ebeam);
            out.hycal_cluster_cfg.energy_bias = out.hycal_energy_bias;

            LOG("[setup] HC E bias : beam=", out.run_cfg.Ebeam,
                " MeV  set=", set.nominal_mev,
                " MeV  ee_cells=", out.hycal_energy_bias->ee_cells_loaded,
                "  ep_cells=", out.hycal_energy_bias->ep_cells_loaded);
            if (out.hycal_energy_bias->ee_cells_loaded == 0 ||
                out.hycal_energy_bias->ep_cells_loaded == 0) {
                LOG("[WARN] HC E bias : one or both parameter files loaded no cells; zero-bias fallback is active.");
            }
        } else {
            LOG("[WARN] HC E bias : enabled but beam energy is invalid; correction is inactive.");
        }
    } else {
        LOG("[setup] HC E bias : disabled");
    }
    LOG("[setup] HC cluster : min_mod_E=", out.hycal_cluster_cfg.min_module_energy,
        "  min_ctr_E=", out.hycal_cluster_cfg.min_center_energy,
        "  min_cl_E=", out.hycal_cluster_cfg.min_cluster_energy,
        "  split_iter=", out.hycal_cluster_cfg.split_iter,
        "  nonlin=", out.hycal_cluster_cfg.non_linear_corr ? "on" : "off",
        "  E_bias=", out.hycal_cluster_cfg.energy_bias_correction ? "on" : "off",
        "  seed_t_win=", out.hycal_cluster_cfg.seed_time_window, "ns",
        out.hycal_cluster_cfg.seed_time_window > 0.f ? " (gated)" : " (off)");

    // --- 10. GEM map (soft — gem stays default-constructed on failure) ---
    const std::string gem_map = pick(gem_map_path_, "gem_map.json");
    if (gem_map.empty()) {
        LOG("[WARN] PipelineBuilder: no GEM map resolved — GEM disabled.");
    } else {
        out.gem.Init(gem_map);
        out.gem_map_path = gem_map;
        LOG("[setup] GEM map    : ", gem_map,
            "  (", out.gem.GetNDetectors(), " detectors)");
    }

    if (!out.gem_crate_remap.empty()) {
        std::ostringstream oss;
        oss << "[setup] GEM crate remap: {";
        bool first = true;
        for (const auto &kv : out.gem_crate_remap) {
            if (!first) oss << ", ";
            oss << kv.first << ": " << kv.second;
            first = false;
        }
        oss << "}";
        LOG(oss.str());
    }

    // --- 11. GEM pedestals -----------------------------------------------
    const std::string ped_path = pick(gem_pedestal_path_, out.run_cfg.gem_pedestal_file);
    if (!ped_path.empty()) {
        out.gem.LoadPedestals(ped_path, out.gem_crate_remap);
        out.gem_pedestal_path = ped_path;
        LOG("[setup] GEM peds   : ", ped_path);
    } else {
        LOG("[WARN] no GEM pedestal file — full-readout data reconstructs empty.");
    }

    // --- 12. GEM common-mode range ---------------------------------------
    const std::string cm_path = pick(gem_common_mode_path_, out.run_cfg.gem_common_mode_file);
    if (!cm_path.empty()) {
        out.gem.LoadCommonModeRange(cm_path, out.gem_crate_remap);
        out.gem_common_mode_path = cm_path;
        LOG("[setup] GEM CM     : ", cm_path);
    }

    // --- 12b. GEM strip-level cuts (sourced from reconstruction_config)
    if (recon.contains("gem") && recon["gem"].is_object()) {
        const auto &gemr = recon["gem"];
        if (gemr.contains("default") && gemr["default"].is_object()) {
            const auto &def = gemr["default"];
            if (def.contains("reject_first_timebin") && def["reject_first_timebin"].is_boolean())
                out.gem.SetRejectFirstTimebin(def["reject_first_timebin"].get<bool>());
            if (def.contains("reject_last_timebin") && def["reject_last_timebin"].is_boolean())
                out.gem.SetRejectLastTimebin(def["reject_last_timebin"].get<bool>());
            if (def.contains("min_peak_adc") && def["min_peak_adc"].is_number())
                out.gem.SetMinPeakAdc(def["min_peak_adc"].get<float>());
            if (def.contains("min_sum_adc") && def["min_sum_adc"].is_number())
                out.gem.SetMinSumAdc(def["min_sum_adc"].get<float>());
        }
    }

    // --- pedestal checksum (matches Python audit's [PEDSUM] line) --------
    if (log_pedestal_checksum_) {
        int n_apvs = out.gem.GetNApvs();
        double sum_noise = 0.0, sum_off = 0.0;
        long n_strips = 0;
        for (int ai = 0; ai < n_apvs; ++ai) {
            const auto &apv = out.gem.GetApvConfig(ai);
            for (int ch = 0; ch < ssp::APV_STRIP_SIZE; ++ch) {
                sum_noise += apv.pedestal[ch].noise;
                sum_off   += apv.pedestal[ch].offset;
                ++n_strips;
            }
        }
        LOG("[PEDSUM] n_apvs=", n_apvs,
            " n_strips=", n_strips,
            " sum_noise=", std::fixed, std::setprecision(6), sum_noise,
            " sum_offset=", sum_off);
    }

    // --- 13. GEMSYS dump (post-Init globals) -----------------------------
    LOG("[GEMSYS] common_mode_thr=", out.gem.GetCommonModeThreshold(),
        " zero_sup_thr=",    out.gem.GetZeroSupThreshold(),
        " cross_talk_thr=",  out.gem.GetCrossTalkThreshold(),
        " min_peak=",        out.gem.GetMinPeakAdc(),
        " min_sum=",         out.gem.GetMinSumAdc(),
        " rej_first=",       (int)out.gem.GetRejectFirstTimebin(),
        " rej_last=",        (int)out.gem.GetRejectLastTimebin());

    // --- 14. GEM per-detector cluster configs ----------------------------
    {
        gem::ClusterConfig def;
        if (recon.contains("gem") && recon["gem"].is_object()) {
            const auto &gemr = recon["gem"];
            if (gemr.contains("default") && gemr["default"].is_object())
                apply_gem_cluster_overrides(gemr["default"], def);
            std::vector<gem::ClusterConfig> per(out.gem.GetNDetectors(), def);
            for (int d = 0; d < out.gem.GetNDetectors(); ++d) {
                std::string key = std::to_string(d);
                if (gemr.contains(key) && gemr[key].is_object())
                    apply_gem_cluster_overrides(gemr[key], per[d]);
            }
            // [GEMCFG] dump (matches Python audit byte-for-byte).
            for (int d = 0; d < (int)per.size(); ++d) {
                const auto &c = per[d];
                LOG("[GEMCFG] d", d,
                    " min_hits=", c.min_cluster_hits,
                    " max_hits=", c.max_cluster_hits,
                    " consec=",   c.consecutive_thres,
                    " split=",    c.split_thres,
                    " xtalk=",    c.cross_talk_width,
                    " xtalk_peak_ratio_max=", c.cross_talk_peak_ratio_max,
                    " match_mode=", c.match_mode,
                    " asym=",     c.match_adc_asymmetry,
                    " tdiff=",    c.match_time_diff,
                    " tperiod=",  c.ts_period,
                    " strip_t=[", c.strip_time_min, ",", c.strip_time_max, "]",
                    " unimodal=", (int)c.strip_unimodal,
                    " seed_peak=", c.seed_min_peak_adc,
                    " seed_sum=", c.seed_min_sum_adc,
                    " strip_dt=", c.strip_time_agreement,
                    " ts_corr=",  c.strip_ts_corr_min);
            }
            out.gem.SetReconConfigs(std::move(per));
        }
    }

    return out;
}

} // namespace prad2
