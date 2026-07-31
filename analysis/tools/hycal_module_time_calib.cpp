// hycal_module_time_calib
//
// Build neighbor-pair HyCal timing delta-t histograms from raw events,
// fit them with Gaussians, solve a connected-module timing-offset network
// with weighted least squares, and write both the histograms and the
// resulting offsets to disk.
//
// Usage:
//   prad2ana_hycal_module_time_calib <input_raw.root|dir> [more...] \
//       -o <output_prefix> [-n max_events] [-v(validation mode)] [-V validation_file]
//
// The output files are derived from the given -o value:
//   - <output_prefix>.root   : ROOT file containing all pair histograms
//   - <output_prefix>.json   : JSON file containing resolved offsets
//
// Global offsets network resolution algorithm developed by Mingyu Li, implemented and optimized by Yuan Li

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "HyCalTimeCalib.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2Poly.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TText.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <memory>
#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <cmath>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#include <nlohmann/json.hpp>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── File collection helper ───────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_raw.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_path_name, daq_config_file;
    int  max_events  = -1;
    int  num_threads = 4;
    bool validation = false;
    std::string validation_file;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:vV:")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'v': validation       = true; break;
            case 'V': validation_file  = optarg; break;
        }
    }

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_module_time_calib <input_raw.root|dir> [more...] "
                     "[-o ./output_name] [-n max_events] "
                     "[-v(validation mode)] [-V validation_file]\n";
        return 1;
    }

    if (output_path_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }

    const fs::path output_prefix(output_path_name);
    const std::string root_output_path =
        (output_prefix.extension() == ".root")
            ? output_prefix.string()
            : output_prefix.string() + ".root";
    const std::string json_output_path =
        (output_prefix.extension() == ".json")
            ? output_prefix.string()
            : output_prefix.string() + ".json";

    TChain tree("events");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }

    EventVars ev;
    prad2::SetRawReadBranches(&tree, ev);

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    std::string calib_file = db_dir + "/" + gRunConfig.energy_calib_file;

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    int nmatched = hycal.LoadCalibration(calib_file);
    std::cerr << "Main: calibration loaded (" << nmatched << " modules)\n";
    if (nmatched <= 0) {
        std::cerr << "Failed to load calibration from " << calib_file
                  << "; aborting before fitting.\n";
        return 1;
    }
    if (validation_file.empty()) prad2::LoadHyCalTimeCalib(db_dir + "/" + "hycal_time_offsets/test.json", hycal);
    else prad2::LoadHyCalTimeCalib(validation_file, hycal);
    analysis::PhysicsTools physics(hycal);
    fdec::HyCalCluster clusterer(hycal);
    fdec::ClusterConfig cl_cfg;
    clusterer.SetConfig(cl_cfg);

    // Create one delta-t histogram per neighboring module pair.
    // Histograms are detached from any ROOT file to avoid ownership issues
    // while we keep them in memory for later filling/fitting.
    using PairKey = std::pair<int, int>;
    auto ordered_pair = [](int a, int b) -> PairKey {
        return (a < b) ? PairKey(a, b) : PairKey(b, a);
    };

    struct PairResult {
        PairKey key;
        int id1 = -1;
        int id2 = -1;
        std::string name;
        std::unique_ptr<TH1F> hist;
        long long entries = 0;
        double mean = 0.0;
        double mean_error = 0.0;
        double sigma = 0.0;
        double chi2_ndf = 0.0;
        bool fit_ok = false;
        bool valid = false;
        bool used = false;
    };

    std::vector<PairResult> pair_results;
    std::map<PairKey, std::size_t> pair_index;
    std::map<PairKey, TH1F*> delta_t_hists;
    const int nbins = 200;
    const float lo = -10.f;
    const float hi = 10.f;

    const int nmod = hycal.module_count();
    for (int m = 0; m < nmod; ++m) {
        auto &mod = hycal.module(m);
        if (!mod.is_pwo4()) continue;

        hycal.for_each_neighbor(m, true, [&](int neighbor_idx) {
            auto &nbr = hycal.module(neighbor_idx);
            if (!nbr.is_pwo4()) return;

            const PairKey key = ordered_pair(mod.id, nbr.id);
            if (delta_t_hists.count(key)) return;

            const std::string name = "h_dt_" + mod.name + "_vs_" + nbr.name;
            const std::string title = mod.name + " - " + nbr.name +
                ";#Delta t (ns);entries";

            auto h = std::make_unique<TH1F>(name.c_str(), title.c_str(), nbins, lo, hi);
            h->SetDirectory(nullptr);

            pair_results.push_back(PairResult{key, mod.id, nbr.id, name, std::move(h)});
            pair_index[key] = pair_results.size() - 1;
            delta_t_hists[key] = pair_results.back().hist.get();
        });
    }

    std::cout << "Created " << delta_t_hists.size() << " neighbor-pair delta-t histograms.\n";

    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev.nch > 70) continue;

        // Reconstruct clusters for this event.
        clusterer.Clear();
        for (int j = 0; j < ev.nch; ++j) {
            const auto *mod = hycal.module_by_id(ev.module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            float adc      = 0.f;
            int bestIdx = -1;
            float bestHeight = -1.f;
            for(int p = 0; p < ev.npeaks[j]; ++p){
                if(ev.peak_time[j][p] > gRunConfig.hc_time_win_lo &&
                    ev.peak_time[j][p] < gRunConfig.hc_time_win_hi) {
                    if(ev.peak_integral[j][p] > bestHeight) {
                        bestHeight = ev.peak_integral[j][p];
                        bestIdx = p;
                    }
                }
            }
            if (bestIdx < 0) continue;
            adc = ev.peak_integral[j][bestIdx] * ev.gain_factor[j]; // apply gain factor
            float energy = (mod->cal_factor > 0)
                ? static_cast<float>(mod->energize(adc)) : 0.f;
            clusterer.AddHit(mod->index, energy, ev.peak_time[j][bestIdx]);
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        // select single cluster Mott events
        if (hits.size() != 1 || hits[0].nblocks <= 3) continue;
        if (hits[0].energy < 0.8 * gRunConfig.Ebeam) continue;

        // Fill all neighboring-module-pair delta-t histograms for this event.
        // For every pair of modules that are both present in the event and are
        // neighbors, compute dt = t_i - t_j and fill the corresponding histogram.
        std::map<int, float> event_times;
        for (int j = 0; j < ev.nch; ++j) {
            const auto *mod = hycal.module_by_id(ev.module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;
            if (ev.npeaks[j] != 1) continue;
            if (ev.peak_height[j][0] < 100) continue;
            if (validation) event_times[mod->index] = ev.peak_time[j][0] - mod->time_offset;
            else event_times[mod->index] = ev.peak_time[j][0];
        }

        for (const auto &[idx1, t1] : event_times) {
            const auto *mod1 = &hycal.module(idx1);
            if (!mod1 || !mod1->is_pwo4()) continue;

            for (const auto &[idx2, t2] : event_times) {
                if (idx2 <= idx1) continue;

                const auto *mod2 = &hycal.module(idx2);
                if (!mod2 || !mod2->is_pwo4()) continue;
                if (!mod1->is_neighbor(*mod2)) continue;

                const PairKey key = ordered_pair(mod1->id, mod2->id);
                auto hist_it = delta_t_hists.find(key);
                if (hist_it == delta_t_hists.end()) continue;

                const float dt = t1 - t2;
                hist_it->second->Fill(dt);
            }
        }
    }

    // Gaussian fit for each neighbor-pair delta-t histogram.
    // Store the fit summary in the in-memory vector for later analysis.
    constexpr double mean_error_floor = 0.01;
    for (auto &row : pair_results) {
        auto *hist = row.hist.get();
        if (!hist) continue;

        row.entries = static_cast<long long>(hist->GetEntries());
        if (row.entries <= 100) {
            continue;
        }

        const double mean = hist->GetMean();
        const double max  = hist->GetMaximum();
        if (!std::isfinite(mean) || !std::isfinite(max) || max <= 0.0) {
            continue;
        }

        const int max_bin = hist->GetMaximumBin();
        const double peak_center = hist->GetXaxis()->GetBinCenter(max_bin);
        const double x_lo = hist->GetXaxis()->GetXmin();
        const double x_hi = hist->GetXaxis()->GetXmax();
        const double fit_window = std::max(1.0, 1.0 * std::max(0.5, hist->GetRMS()));

        const std::vector<std::pair<double, double>> candidate_ranges = {
            {std::max(x_lo, peak_center - fit_window), std::min(x_hi, peak_center + fit_window)},
            {x_lo, x_hi}
        };

        bool fitted = false;
        for (const auto &[fit_lo, fit_hi] : candidate_ranges) {
            if (!std::isfinite(fit_lo) || !std::isfinite(fit_hi) || fit_hi <= fit_lo) {
                continue;
            }

            auto fit_result = hist->Fit("gaus", "QRS", "", fit_lo, fit_hi);
            const int fit_status = fit_result;
            if (fit_status != 0 || !fit_result.Get()) {
                continue;
            }

            const double mu = fit_result->Parameter(1);
            const double sigma = fit_result->Parameter(2);
            const double ndf = static_cast<double>(fit_result->Ndf());
            const double chi2_ndf = (ndf > 0.0) ? fit_result->Chi2() / ndf : std::numeric_limits<double>::infinity();

            if (!std::isfinite(mu) || !std::isfinite(sigma) || sigma <= 0.0) {
                continue;
            }

            row.fit_ok = true;
            row.mean = mu;
            row.sigma = sigma;
            row.mean_error = std::isfinite(hist->GetMeanError())
                ? hist->GetMeanError()
                : mean_error_floor;
            row.mean_error = std::max(row.mean_error, mean_error_floor);
            row.chi2_ndf = chi2_ndf;
            row.valid = std::isfinite(row.chi2_ndf);
            fitted = true;
            break;
        }

        if (!fitted) {
            continue;
        }
    }

    // cout the pair results valid number out of the total pair number
    std::size_t valid_count = 0;
    for (const auto &row : pair_results) {
        if (row.valid) ++valid_count;
    }
    std::cout << "Valid pair results: " << valid_count << " / " << pair_results.size() << "\n";

    // Build the valid adjacency graph from the fitted pair means.
    std::map<int, std::vector<int>> adjacency;
    for (const auto &row : pair_results) {
        if (!row.valid) continue;

        adjacency[row.id1].push_back(row.id2);
        adjacency[row.id2].push_back(row.id1);
    }

    int reference_module = 1495;
    /*std::size_t max_degree = 0;
    for (const auto &kv : adjacency) {
        if (kv.second.size() > max_degree) {
            max_degree = kv.second.size();
            reference_module = kv.first;
        }
    }*/

    if (reference_module < 0) {
        std::cerr << "No valid timing network could be built from the fitted pair means.\n";
        return 1;
    }

    std::set<int> connected;
    std::queue<int> q;
    connected.insert(reference_module);
    q.push(reference_module);

    while (!q.empty()) {
        const int current = q.front();
        q.pop();

        const auto it = adjacency.find(current);
        if (it == adjacency.end()) continue;

        for (const int next : it->second) {
            if (connected.insert(next).second) {
                q.push(next);
            }
        }
    }

    std::cout << "Connected timing network size = " << connected.size() << "\n";

    std::map<int, double> module_offsets;
    std::map<int, double> module_offset_errors;
    module_offsets[reference_module] = 0.0;
    module_offset_errors[reference_module] = 0.0;

    if (connected.size() < 2) {
        std::cerr << "Reference module has no usable connected timing network; writing default zero offsets.\n";
    } else {
        std::map<int, int> index_by_module;
        std::vector<int> modules_in_network;
        for (const int module_id : connected) {
            if (module_id == reference_module) continue;
            index_by_module[module_id] = static_cast<int>(modules_in_network.size());
            modules_in_network.push_back(module_id);
        }

        const int n_unknown = static_cast<int>(modules_in_network.size());
        TMatrixD normal(n_unknown, n_unknown);
        TVectorD rhs(n_unknown);
        normal.Zero();
        rhs.Zero();

        long long n_used_equations = 0;
        for (const auto &row : pair_results) {
            if (!row.valid) continue;
            if (!connected.count(row.id1) || !connected.count(row.id2)) continue;

            const double sigma = std::max(row.mean_error, mean_error_floor);
            const double weight = 1.0 / (sigma * sigma);

            const bool is_ref1 = (row.id1 == reference_module);
            const bool is_ref2 = (row.id2 == reference_module);

            const int i = is_ref1 ? -1 : index_by_module[row.id1];
            const int j = is_ref2 ? -1 : index_by_module[row.id2];

            if (!is_ref1) {
                normal(i, i) += weight;
                rhs(i) += weight * row.mean;
            }

            if (!is_ref2) {
                normal(j, j) += weight;
                rhs(j) -= weight * row.mean;
            }

            if (!is_ref1 && !is_ref2) {
                normal(i, j) -= weight;
                normal(j, i) -= weight;
            }

            ++n_used_equations;
        }

        std::cout << "Weighted least-squares equations = " << n_used_equations << "\n";
        std::cout << "Unknown offsets = " << n_unknown << "\n";

        if (n_unknown > 0) {
            TDecompSVD svd(normal);
            Bool_t solve_ok = kFALSE;
            TVectorD solution = svd.Solve(rhs, solve_ok);
            if (solve_ok) {
                Bool_t inverse_ok = kFALSE;
                TMatrixD covariance = svd.Invert(inverse_ok);

                for (int i = 0; i < n_unknown; ++i) {
                    const int module_id = modules_in_network[i];
                    module_offsets[module_id] = solution(i);

                    if (inverse_ok && covariance(i, i) >= 0.0) {
                        module_offset_errors[module_id] = std::sqrt(covariance(i, i));
                    } else {
                        module_offset_errors[module_id] = -1.0;
                    }
                }
            } else {
                std::cerr << "Weighted least-squares solve failed; keeping default zero offsets.\n";
            }
        }
    }

    std::cout << "Resolved offsets relative to module " << reference_module << "\n";
    for (const auto &kv : module_offsets) {
        std::cout << "  module " << kv.first
                  << " -> offset = " << std::fixed << std::setprecision(4)
                  << kv.second << " ns"
                  << " (err = " << module_offset_errors[kv.first] << ")\n";
    }

    if (!json_output_path.empty()) {
        const fs::path out_path(json_output_path);
        if (!out_path.parent_path().empty()) {
            fs::create_directories(out_path.parent_path());
        }

        nlohmann::ordered_json out_json = nlohmann::ordered_json::object();
        out_json["reference_module"] = reference_module;
        out_json["reference_name"] = hycal.module_by_id(reference_module)
            ? hycal.module_by_id(reference_module)->name
            : "";
        out_json["offset_units"] = "ns";
        out_json["default"] = 0.0;
        out_json["modules"] = nlohmann::ordered_json::array();

        for (int m = 0; m < hycal.module_count(); ++m) {
            auto &mod = hycal.module(m);
            if (!mod.is_pwo4()) continue;

            const int module_id = mod.id;
            const auto *mod_ptr = hycal.module_by_id(module_id);
            const bool has_offset = module_offsets.count(module_id) > 0;
            const double offset_value = has_offset ? module_offsets.at(module_id) : 0.0;
            const double error_value = has_offset ? module_offset_errors.at(module_id) : 0.0;

            nlohmann::ordered_json entry = nlohmann::ordered_json::object();
            entry["name"] = mod_ptr ? mod_ptr->name : "";
            entry["module_id"] = module_id;
            entry["offset_ns"] = offset_value;
            entry["offset_error_ns"] = error_value;
            entry["connected"] = connected.count(module_id) > 0;
            entry["resolved"] = has_offset;
            out_json["modules"].push_back(entry);
        }

        std::ofstream json_out(out_path);
        if (json_out.is_open()) {
            json_out << std::setw(2) << out_json << "\n";
            json_out.close();
            std::cout << "Wrote offsets JSON to " << out_path.string() << "\n";
        } else {
            std::cerr << "Failed to write offsets JSON to " << out_path.string() << "\n";
        }
    }

    if (!root_output_path.empty()) {
        const fs::path root_path(root_output_path);
        if (!root_path.parent_path().empty()) {
            fs::create_directories(root_path.parent_path());
        }

        TFile root_file(root_path.string().c_str(), "RECREATE");
        if (root_file.IsOpen()) {
            constexpr double map_min = -8.0;
            constexpr double map_max = 8.0;
            auto *h_offset_map = new TH2Poly(
                "h_offset_map",
                "HyCal timing offsets; x (mm); y (mm); offset (ns)",
                -500.0,
                500.0,
                -500.0,
                500.0
            );
            h_offset_map->SetDirectory(&root_file);
            h_offset_map->SetContour(255);
            h_offset_map->SetStats(false);
            h_offset_map->SetMinimum(map_min);
            h_offset_map->SetMaximum(map_max);
            h_offset_map->SetOption("colz");
            h_offset_map->SetTitleOffset(1.2, "Y");
            h_offset_map->SetTitleOffset(1.2, "X");
            h_offset_map->GetXaxis()->SetNoExponent(true);
            h_offset_map->GetYaxis()->SetNoExponent(true);
            h_offset_map->GetXaxis()->SetNdivisions(505);
            h_offset_map->GetYaxis()->SetNdivisions(505);
            h_offset_map->SetMarkerSize(0.0);
            h_offset_map->SetMarkerColor(kBlack);
            h_offset_map->SetLineColor(kBlack);
            h_offset_map->SetLineWidth(1);

            std::map<int, int> offset_bin_map;
            int bad_poly_bins = 0;
            for (int m = 0; m < hycal.module_count(); ++m) {
                auto &mod = hycal.module(m);
                if (!mod.is_pwo4()) continue;

                const double x1 = mod.x - 0.5 * mod.size_x;
                const double x2 = mod.x + 0.5 * mod.size_x;
                const double y1 = mod.y - 0.5 * mod.size_y;
                const double y2 = mod.y + 0.5 * mod.size_y;

                const int bin_id = h_offset_map->AddBin(x1, y1, x2, y2);
                offset_bin_map[mod.id] = bin_id;
                if (bin_id <= 0) {
                    ++bad_poly_bins;
                }
            }

            int clamped_low = 0;
            int clamped_high = 0;
            for (int m = 0; m < hycal.module_count(); ++m) {
                auto &mod = hycal.module(m);
                if (!mod.is_pwo4()) continue;

                const auto it = offset_bin_map.find(mod.id);
                if (it == offset_bin_map.end() || it->second <= 0) continue;

                const bool has_offset = module_offsets.count(mod.id) > 0;
                const double offset_value = has_offset ? module_offsets.at(mod.id) : 0.0;
                if (offset_value < map_min) ++clamped_low;
                if (offset_value > map_max) ++clamped_high;
                h_offset_map->SetBinContent(it->second, std::clamp(offset_value, map_min, map_max));
            }

            if (bad_poly_bins > 0) {
                std::cerr << "Warning: TH2Poly failed to create " << bad_poly_bins
                          << " HyCal map bins.\n";
            }
            if (clamped_low > 0 || clamped_high > 0) {
                std::cerr << "HyCal map: clamped " << clamped_low << " bins below "
                          << map_min << " and " << clamped_high << " bins above "
                          << map_max << ".\n";
            }

            h_offset_map->Write();

            for (const auto &row : pair_results) {
                if (!row.hist) continue;
                row.hist->Write();
            }

            root_file.Close();
            std::cout << "Wrote offset map and pair histograms to " << root_path.string() << "\n";
        } else {
            std::cerr << "Failed to write ROOT histograms to " << root_path.string() << "\n";
        }
    }

    return 0;
}