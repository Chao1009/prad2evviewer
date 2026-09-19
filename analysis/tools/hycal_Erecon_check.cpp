
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TClass.h>
#include <TLorentzVector.h>

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <limits>
#include <thread>
#include <getopt.h>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

//anagles bin edges
const int Nbins = 33;
const Double_t binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};
const int energy_bins = 500; const double energy_min = 0., energy_max = 5000.;

struct HistResult {
    std::unique_ptr<TH2F> h2_hit_module_hycal;
    std::unique_ptr<TH2F> h2_hit_module_gem;
    std::unique_ptr<TH2F> h2_hit_hycal;
    std::unique_ptr<TH2F> h2_hit_gem;
    std::unique_ptr<TH1F> h1_yield_theta;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta_center;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta_edge;
    std::vector<std::unique_ptr<TH1F>> h1_E_modules;
    long long events_processed = 0;
};

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix);
static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result);

static void mergeHistResult(HistResult &destination, const HistResult &source)
{
    destination.h2_hit_module_hycal->Add(source.h2_hit_module_hycal.get());
    destination.h2_hit_module_gem->Add(source.h2_hit_module_gem.get());
    destination.h2_hit_hycal->Add(source.h2_hit_hycal.get());
    destination.h2_hit_gem->Add(source.h2_hit_gem.get());
    destination.h1_yield_theta->Add(source.h1_yield_theta.get());
    for (int i = 0; i < Nbins; ++i) {
        destination.h1_Erecon_theta[i]->Add(source.h1_Erecon_theta[i].get());
        destination.h1_Erecon_theta_center[i]->Add(source.h1_Erecon_theta_center[i].get());
        destination.h1_Erecon_theta_edge[i]->Add(source.h1_Erecon_theta_edge[i].get());
    }
    for (int i = 0; i < 1156; ++i) {
        destination.h1_E_modules[i]->Add(source.h1_E_modules[i].get());
    }
    destination.events_processed += source.events_processed;
}

static std::vector<std::string> collectRootFiles(const std::string &path);

static std::string shell_quote(const std::string &value)
{
    std::string quoted = "'";
    for (char ch : value) {
        if (ch == '\'') quoted += "'\\''";
        else quoted += ch;
    }
    return quoted + "'";
}

static std::string outputFileName(const std::string &output_name, bool corr = false)
{
    const fs::path output_path(output_name);
    const std::string file_name = output_path.filename().string()
        + (corr ? ".corr" : "") + ".root";
    return (output_path.parent_path() / file_name).string();
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            std::string name = entry.path().filename().string();
            if (entry.is_regular_file() &&
                name.find("_recon") != std::string::npos &&
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}


bool inHyCal(float xmm, float ymm) {
    const float module = 20.75; // mm
    return (fabs(xmm) > module * 2.5 || fabs(ymm) > module * 2.5)
        && (fabs(xmm) < module * 15. && fabs(ymm) < module * 15.);
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
    std::string output_name;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;
    bool worker_mode = false;
    bool corr = false;

    static option long_options[] = {
        {"corr", no_argument, nullptr, 'c'},
        {nullptr, 0, nullptr, 0}
    };
    int opt;
    while ((opt = getopt_long_only(argc, argv, "o:n:f:j:",
                                   long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
        }
    }

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        if (num_files > 0) {
            int remaining = num_files - static_cast<int>(root_files.size());
            if (remaining <= 0) break;
            int take = std::min(remaining, static_cast<int>(f.size()));
            root_files.insert(root_files.end(), f.begin(), f.begin() + take);
            if (static_cast<int>(root_files.size()) >= num_files) break;
        } else {
            root_files.insert(root_files.end(), f.begin(), f.end());
        }
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_Erecon_check <input_recon.root|dir> [more...] "
                 "-o <output_name> [-n max_events] [-f nfiles] [-j threads] [-corr]\n";
        return 1;
    }

    if (output_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TH1F");
    TClass::GetClass("TH2F");

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    std::vector<long long> file_limits(root_files.size(), -1);
    if (max_events >= 0) {
        long long remaining = max_events;
        for (size_t i = 0; i < root_files.size(); ++i) {
            TFile input_file(root_files[i].c_str(), "READ");
            auto *tree = input_file.IsOpen()
                ? dynamic_cast<TTree *>(input_file.Get("recon")) : nullptr;
            const long long entries = tree ? tree->GetEntries() : 0;
            file_limits[i] = std::max(0LL, std::min(entries, remaining));
            remaining -= file_limits[i];
        }
    }

    auto merged = makeHistResult("");
    const int threads_count = std::max(1, std::min(num_threads,
        static_cast<int>(root_files.size())));
    const int rounds = (static_cast<int>(root_files.size()) + threads_count - 1)
        / threads_count;
    std::mutex io_mutex;
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << threads_count << " thread(s), " << rounds << " round(s)\n";

    for (int round = 0; round < rounds; ++round) {
        const int first = round * threads_count;
        const int last = std::min(first + threads_count,
                                  static_cast<int>(root_files.size()));
        std::vector<std::unique_ptr<HistResult>> results(last - first);
        std::vector<std::thread> workers;
        workers.reserve(last - first);

        for (int file_index = first; file_index < last; ++file_index) {
            workers.emplace_back([&, file_index, first]() {
                auto result = makeHistResult(Form("worker%d", file_index - first));
                const long long limit = max_events >= 0 ? file_limits[file_index] : -1;
                const bool ok = processRootFile(root_files[file_index], gRunConfig,
                                                db_dir, limit, result.get());
                results[file_index - first] = std::move(result);
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cout << "[worker " << (file_index - first) << "] file "
                          << file_index << " / " << (root_files.size() - 1)
                          << ": " << root_files[file_index] << " -> "
                          << (ok ? "OK" : "FAILED") << "\n";
            });
        }
        for (auto &worker : workers) worker.join();
        for (const auto &result : results) {
            if (result) mergeHistResult(*merged, *result);
        }
    }

    // Fit the merged histograms, fill new histograms with the fit results
    auto h1_deltaE_module = std::make_unique<TH1F>("h1_deltaE_module", "Delta E (Module);#Delta E [MeV];Modules",
        1000, -50, 50);
    auto h2_deltaE_module_map = std::make_unique<TH2Poly>(
        "h2_deltaE_module_map", "Delta E (Module) Map;X [mm];Y [mm];#Delta E [MeV]",
        -360., 360., -360., 360.);

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    analysis::PhysicsTools physics(hycal);

    std::vector<int> deltaE_polygon_bins(1156, -1);
    for (int i = 0; i < 1156; ++i) {
        const auto *mod = hycal.module_by_id(i + 1001);
        if (!mod || mod->size_x <= 0. || mod->size_y <= 0.) continue;
        deltaE_polygon_bins[i] = h2_deltaE_module_map->AddBin(
            mod->x - 0.5 * mod->size_x, mod->y - 0.5 * mod->size_y,
            mod->x + 0.5 * mod->size_x, mod->y + 0.5 * mod->size_y);
    }
    h2_deltaE_module_map->GetZaxis()->SetTitle("#Delta E [MeV]");
    h2_deltaE_module_map->SetMinimum(-10.);
    h2_deltaE_module_map->SetMaximum(10.);
    h2_deltaE_module_map->SetOption("colz");

    for (int i = 0; i < 1156; i++) {
        const auto *mod = hycal.module_by_id(i + 1001);
        if (!mod) continue;

        TH1F *energy_histogram = merged->h1_E_modules[i].get();
        if (energy_histogram->GetEntries() < 200) continue;

        const float theta = std::atan2(
            std::sqrt(mod->x * mod->x + mod->y * mod->y),
            gRunConfig.hycal_z) * 180.0f / M_PI;
        const float expected_energy = analysis::PhysicsTools::ExpectedEnergy(
            theta, gRunConfig.Ebeam, "ep");
        auto fit_result = physics.fitPeak(energy_histogram, expected_energy, true);
        const float peak = static_cast<float>(fit_result[0]);
        const float sigma = static_cast<float>(fit_result[1]);
        const float chi2 = static_cast<float>(fit_result[2]);
        const float peak_error = static_cast<float>(fit_result[3]);
        const float sigma_error = static_cast<float>(fit_result[4]);
        const float expected_sigma = peak > 0.f
            ? 0.03f * peak / std::sqrt(peak / 1000.f) : 0.f;
        const bool fit_good = peak > 0.f && expected_sigma > 0.f
            && sigma > 0.5f * expected_sigma
            && sigma < 1.5f * expected_sigma
            && chi2 < 2.5f;
        if (fit_good) {
            const float delta_energy = peak - expected_energy;
            h1_deltaE_module->Fill(delta_energy);

            // Each module has its own polygon bin; assign rather than accumulate.
            const int polygon_bin = deltaE_polygon_bins[i];
            if (polygon_bin > 0)
                h2_deltaE_module_map->SetBinContent(polygon_bin, delta_energy);
        }
    }
    std::vector<double> angle_points;
    std::vector<double> angle_center_points;
    std::vector<double> angle_edge_points;
    std::vector<double> expected_angle_points;
    std::vector<double> energy_points;
    std::vector<double> energy_center_points;
    std::vector<double> energy_edge_points;
    std::vector<double> expect_energy_points;
    std::vector<double> sigma_points;
    std::vector<double> sigma_center_points;
    std::vector<double> sigma_edge_points;
    std::vector<double> peak_error_points;
    std::vector<double> peak_error_center_points;
    std::vector<double> peak_error_edge_points;
    std::vector<double> sigma_error_points;
    std::vector<double> sigma_error_center_points;
    std::vector<double> sigma_error_edge_points;
    std::vector<double> resolution_points;
    std::vector<double> resolution_center_points;
    std::vector<double> resolution_edge_points;
    std::vector<double> resolution_error_points;
    std::vector<double> resolution_error_center_points;
    std::vector<double> resolution_error_edge_points;
    for (int i = 0; i < Nbins; ++i) {
        const double angle = 0.5 * (binEdge[i] + binEdge[i + 1]);
        const double expected_energy = analysis::PhysicsTools::ExpectedEnergy(
            angle, gRunConfig.Ebeam, "ep");
        const auto fit_histogram = [&](TH1F *hist) {
            return hist->GetEntries() >= 200
                ? physics.fitPeak(hist, static_cast<float>(expected_energy), true)
                : std::array<double, 5>{0., 0., 0., 0., 0.};
        };
        const auto fit_all = fit_histogram(merged->h1_Erecon_theta[i].get());
        const auto fit_center = fit_histogram(merged->h1_Erecon_theta_center[i].get());
        const auto fit_edge = fit_histogram(merged->h1_Erecon_theta_edge[i].get());

        const auto add_fit = [&](const std::array<double, 5> &fit,
                                 std::vector<double> &energy,
                                 std::vector<double> &sigma,
                                 std::vector<double> &peak_error,
                                 std::vector<double> &sigma_error,
                                 std::vector<double> &resolution,
                                 std::vector<double> &resolution_error,
                                 std::vector<double> &angles) {
            if (fit[0] <= 0. || fit[1] <= 0.) return false;
            angles.push_back(angle);
            energy.push_back(fit[0]);
            sigma.push_back(fit[1]);
            peak_error.push_back(fit[3]);
            sigma_error.push_back(fit[4]);
            resolution.push_back(100.0 * fit[1] / fit[0] * std::sqrt(fit[0] / 1000.0));
            resolution_error.push_back(100.0 * fit[4] / fit[0] * std::sqrt(fit[0] / 1000.0));
            return true;
        };

        const bool has_all = add_fit(fit_all, energy_points, sigma_points,
                                     peak_error_points, sigma_error_points,
                                                                         resolution_points, resolution_error_points,
                                                                         angle_points);
        const bool has_center = add_fit(fit_center, energy_center_points, sigma_center_points,
                                        peak_error_center_points, sigma_error_center_points,
                                                                                resolution_center_points, resolution_error_center_points,
                                                                                angle_center_points);
        const bool has_edge = add_fit(fit_edge, energy_edge_points, sigma_edge_points,
                                      peak_error_edge_points, sigma_error_edge_points,
                                                                            resolution_edge_points, resolution_error_edge_points,
                                                                            angle_edge_points);
        if (has_all || has_center || has_edge) {
            expected_angle_points.push_back(angle);
            expect_energy_points.push_back(expected_energy);
        }
    }
    TCanvas c_reconE("c_reconE", "Reconstructed energy and resolution vs angle", 900, 800);
    TPad upper_pad("upper_pad", "Reconstructed energy", 0.0, 0.35, 1.0, 1.0);
    TPad lower_pad("lower_pad", "Resolution", 0.0, 0.0, 1.0, 0.35);
    upper_pad.SetTopMargin(0.08);
    upper_pad.SetBottomMargin(0.0);
    lower_pad.SetTopMargin(0.0);
    lower_pad.SetBottomMargin(0.20);
    upper_pad.Draw();
    lower_pad.Draw();

    TGraphErrors gErecon_vs_theta(static_cast<int>(energy_points.size()));
    TGraphErrors gErecon_vs_theta_center(static_cast<int>(energy_center_points.size()));
    TGraphErrors gErecon_vs_theta_edge(static_cast<int>(energy_edge_points.size()));
    TGraph gExpectedE_vs_theta(static_cast<int>(expected_angle_points.size()));
    TGraphErrors gResolution_vs_theta(static_cast<int>(resolution_points.size()));
    TGraphErrors gResolution_vs_theta_center(static_cast<int>(resolution_center_points.size()));
    TGraphErrors gResolution_vs_theta_edge(static_cast<int>(resolution_edge_points.size()));
    for (size_t i = 0; i < expected_angle_points.size(); ++i) {
        gExpectedE_vs_theta.SetPoint(i, expected_angle_points[i], expect_energy_points[i]);
    }
    for (size_t i = 0; i < energy_points.size(); ++i) {
        gErecon_vs_theta.SetPoint(i, angle_points[i], energy_points[i]);
        gErecon_vs_theta.SetPointError(i, 0.0, peak_error_points[i]);
        gResolution_vs_theta.SetPoint(i, angle_points[i], resolution_points[i]);
        gResolution_vs_theta.SetPointError(i, 0.0, resolution_error_points[i]);
    }
    for (size_t i = 0; i < energy_center_points.size(); ++i) {
        gErecon_vs_theta_center.SetPoint(i, angle_center_points[i], energy_center_points[i]);
        gErecon_vs_theta_center.SetPointError(i, 0.0, peak_error_center_points[i]);
        gResolution_vs_theta_center.SetPoint(i, angle_center_points[i], resolution_center_points[i]);
        gResolution_vs_theta_center.SetPointError(i, 0.0, resolution_error_center_points[i]);
    }
    for (size_t i = 0; i < energy_edge_points.size(); ++i) {
        gErecon_vs_theta_edge.SetPoint(i, angle_edge_points[i], energy_edge_points[i]);
        gErecon_vs_theta_edge.SetPointError(i, 0.0, peak_error_edge_points[i]);
        gResolution_vs_theta_edge.SetPoint(i, angle_edge_points[i], resolution_edge_points[i]);
        gResolution_vs_theta_edge.SetPointError(i, 0.0, resolution_error_edge_points[i]);
    }

    const auto axis_range = [](const std::vector<double> &values,
                               const std::vector<double> &errors,
                               double fallback_min, double fallback_max) {
        double min_value = std::numeric_limits<double>::max();
        double max_value = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < values.size(); ++i) {
            const double error = i < errors.size() ? std::abs(errors[i]) : 0.0;
            min_value = std::min(min_value, values[i] - error);
            max_value = std::max(max_value, values[i] + error);
        }
        if (!(max_value > min_value)) return std::pair<double, double>{fallback_min, fallback_max};
        const double padding = 0.08 * (max_value - min_value);
        return std::pair<double, double>{min_value - padding, max_value + padding};
    };

    double energy_min_value = std::numeric_limits<double>::max();
    double energy_max_value = std::numeric_limits<double>::lowest();
    const auto include_energy_range = [&](const std::vector<double> &values,
                                          const std::vector<double> &errors) {
        const auto range = axis_range(values, errors, 0.0, 1.0);
        energy_min_value = std::min(energy_min_value, range.first);
        energy_max_value = std::max(energy_max_value, range.second);
    };
    include_energy_range(energy_points, peak_error_points);
    include_energy_range(energy_center_points, peak_error_center_points);
    include_energy_range(energy_edge_points, peak_error_edge_points);
    include_energy_range(expect_energy_points, {});
    if (!(energy_max_value > energy_min_value)) {
        energy_min_value = 0.0;
        energy_max_value = 1.0;
    }

    double resolution_min_value = std::numeric_limits<double>::max();
    double resolution_max_value = std::numeric_limits<double>::lowest();
    const auto include_resolution_range = [&](const std::vector<double> &values,
                                               const std::vector<double> &errors) {
        const auto range = axis_range(values, errors, 0.0, 1.0);
        resolution_min_value = std::min(resolution_min_value, range.first);
        resolution_max_value = std::max(resolution_max_value, range.second);
    };
    include_resolution_range(resolution_points, resolution_error_points);
    include_resolution_range(resolution_center_points, resolution_error_center_points);
    include_resolution_range(resolution_edge_points, resolution_error_edge_points);
    if (!(resolution_max_value > resolution_min_value)) {
        resolution_min_value = 0.0;
        resolution_max_value = 1.0;
    }

    upper_pad.cd();
    gErecon_vs_theta.SetTitle("Reconstructed energy vs angle; ;E_{recon} [MeV]");
    gErecon_vs_theta.SetMarkerStyle(20);
    gErecon_vs_theta.SetMarkerColor(kBlue + 1);
    gErecon_vs_theta.SetLineColor(kBlue + 1);
    gErecon_vs_theta.GetXaxis()->SetLimits(binEdge[0], binEdge[Nbins]);
    gErecon_vs_theta.GetXaxis()->SetLabelSize(0.0);
    gErecon_vs_theta.GetXaxis()->SetTitleSize(0.0);
    gErecon_vs_theta.GetYaxis()->SetLabelSize(0.045);
    gErecon_vs_theta.GetYaxis()->SetTitleSize(0.045);
    gErecon_vs_theta.GetYaxis()->SetTitleOffset(0.81);
    gErecon_vs_theta.GetYaxis()->CenterTitle();
    gErecon_vs_theta.SetMinimum(energy_min_value);
    gErecon_vs_theta.SetMaximum(energy_max_value);
    gErecon_vs_theta.Draw("AP");
    gErecon_vs_theta_center.SetMarkerStyle(21);
    gErecon_vs_theta_center.SetMarkerColor(kGreen + 2);
    gErecon_vs_theta_center.SetLineColor(kGreen + 2);
    gErecon_vs_theta_center.Draw("P SAME");
    gErecon_vs_theta_edge.SetMarkerStyle(22);
    gErecon_vs_theta_edge.SetMarkerColor(kOrange + 1);
    gErecon_vs_theta_edge.SetLineColor(kOrange + 1);
    gErecon_vs_theta_edge.Draw("P SAME");
    gExpectedE_vs_theta.SetLineColor(kRed + 1);
    gExpectedE_vs_theta.SetLineWidth(2);
    gExpectedE_vs_theta.Draw("L SAME");
    TLegend energy_legend(0.62, 0.78, 0.93, 0.92);
    energy_legend.AddEntry(&gErecon_vs_theta, "All", "ep");
    energy_legend.AddEntry(&gErecon_vs_theta_center, "Center", "p");
    energy_legend.AddEntry(&gErecon_vs_theta_edge, "Edge", "p");
    energy_legend.AddEntry(&gExpectedE_vs_theta, "Expected E", "l");
    energy_legend.Draw();

    lower_pad.cd();
    gResolution_vs_theta.SetTitle("Energy resolution vs angle;Scattering angle [deg];#sigma / E * #sqrt{E[GeV]} [%]");
    gResolution_vs_theta.SetMarkerStyle(20);
    gResolution_vs_theta.SetMarkerColor(kBlue + 1);
    gResolution_vs_theta.SetLineColor(kBlue + 1);
    gResolution_vs_theta.GetXaxis()->SetLimits(binEdge[0], binEdge[Nbins]);
    gResolution_vs_theta.GetXaxis()->SetLabelSize(0.060);
    gResolution_vs_theta.GetXaxis()->SetTitleSize(0.070);
    gResolution_vs_theta.GetXaxis()->SetTitleOffset(1.05);
    gResolution_vs_theta.GetYaxis()->SetLabelSize(0.08);
    gResolution_vs_theta.GetYaxis()->SetTitleSize(0.08);
    gResolution_vs_theta.GetYaxis()->SetTitleOffset(0.58);
    gResolution_vs_theta.GetYaxis()->CenterTitle();
    gResolution_vs_theta.SetMinimum(resolution_min_value);
    gResolution_vs_theta.SetMaximum(resolution_max_value);
    gResolution_vs_theta.Draw("AP");
    gResolution_vs_theta_center.SetMarkerStyle(21);
    gResolution_vs_theta_center.SetMarkerColor(kGreen + 2);
    gResolution_vs_theta_center.SetLineColor(kGreen + 2);
    gResolution_vs_theta_center.Draw("P SAME");
    gResolution_vs_theta_edge.SetMarkerStyle(22);
    gResolution_vs_theta_edge.SetMarkerColor(kOrange + 1);
    gResolution_vs_theta_edge.SetLineColor(kOrange + 1);
    gResolution_vs_theta_edge.Draw("P SAME");
    c_reconE.Update();


    const std::string output_file_name = outputFileName(output_name, corr);
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Cannot create output file " << output_file_name << "\n";
        return 1;
    }
    merged->h2_hit_hycal->Write();
    merged->h2_hit_gem->Write();
    merged->h2_hit_module_hycal->Write();
    merged->h2_hit_module_gem->Write();
    merged->h1_yield_theta->Write();
    h1_deltaE_module->Write();
    h2_deltaE_module_map->Write();
    c_reconE.Write("recon_energy_resolution_vs_angle");
    output_file.cd();
    output_file.mkdir("Erecon_theta");
    output_file.cd("Erecon_theta");
    for (int i = 0; i < Nbins; ++i) {
        merged->h1_Erecon_theta[i]->Write();
        merged->h1_Erecon_theta_center[i]->Write();
        merged->h1_Erecon_theta_edge[i]->Write();
    }
    output_file.cd();
    output_file.mkdir("E_modules");
    output_file.cd("E_modules");
    for (int i = 0; i < 1156; ++i) {
        merged->h1_E_modules[i]->Write();
    }
    output_file.cd();
    output_file.Close();

}

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix)
{
    auto result = std::make_unique<HistResult>();
    const std::string name_suffix = suffix.empty() ? "" : "_" + suffix;
    result->h2_hit_module_hycal = std::make_unique<TH2F>(
        Form("h2_hit_module_hycal%s", name_suffix.c_str()),
        "HyCal Hit Distribution (Module);(X_{hycal}-X_{cell center})/d_{cell size};(Y_{hycal}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    result->h2_hit_module_gem = std::make_unique<TH2F>(
        Form("h2_hit_module_gem%s", name_suffix.c_str()),
        "GEM Hit Distribution (Module);(X_{gem}-X_{cell center})/d_{cell size};(Y_{gem}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    result->h2_hit_hycal = std::make_unique<TH2F>(
        Form("h2_hit_hycal%s", name_suffix.c_str()),
        "HyCal Hit Distribution;X [mm];Y [mm]", 720, -360, 360, 720, -360, 360);
    result->h2_hit_gem = std::make_unique<TH2F>(
        Form("h2_hit_gem%s", name_suffix.c_str()),
        "GEM Hit Distribution;X [mm];Y [mm]", 720, -360, 360, 720, -360, 360);
    result->h1_yield_theta = std::make_unique<TH1F>(
        Form("h1_yield_theta%s", name_suffix.c_str()),
        "Yield vs Theta;#theta [deg];Yield/binWidth", Nbins, binEdge);

    result->h1_Erecon_theta.reserve(Nbins);
    result->h1_Erecon_theta_center.reserve(Nbins);
    result->h1_Erecon_theta_edge.reserve(Nbins);
    for (int i = 0; i < Nbins; ++i) {
        result->h1_Erecon_theta.push_back(std::make_unique<TH1F>(
            Form("h1_Erecon_theta_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg];E_{recon} [MeV];Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
        result->h1_Erecon_theta_center.push_back(std::make_unique<TH1F>(
            Form("h1_Erecon_theta_center_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg] (Center);E_{recon} [MeV]; Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
        result->h1_Erecon_theta_edge.push_back(std::make_unique<TH1F>(
            Form("h1_Erecon_theta_edge_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg] (Edge);E_{recon} [MeV]; Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
    }

    result->h1_E_modules.reserve(1156);
    for (int i = 0; i < 1156; ++i) {
        const int mod_id = i + 1001;
        result->h1_E_modules.push_back(std::make_unique<TH1F>(
            Form("h1_E_mod_%d%s", mod_id, name_suffix.c_str()),
            Form("Module W%d cluster energy;E (MeV);Counts", mod_id - 1000),
            energy_bins, energy_min, energy_max));
    }
    return result;
}

static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result)
{
    TFile input(input_file.c_str(), "READ");
    if (input.IsZombie() || !input.IsOpen()) {
        std::cerr << "Cannot open input file " << input_file << "\n";
        return false;
    }
    auto *tree = dynamic_cast<TTree *>(input.Get("recon"));
    if (!tree) {
        std::cerr << "Input file has no 'recon' tree: " << input_file << "\n";
        return false;
    }

    EventVars_Recon event;
    prad2::SetReconReadBranches(tree, event);
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    Long64_t entries = tree->GetEntries();
    if (max_events >= 0) entries = std::min(entries, max_events);
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if ((event.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (event.n_clusters != 1) continue;
        //if (event.matchNum != 1) continue;
        if (event.cl_nblocks[0] < 3) continue;

        HCHit hc_hit;
        GEMHit gem_hit;
        hc_hit.x = event.cl_x[0];
        hc_hit.y = event.cl_y[0];
        hc_hit.z = event.cl_z[0];
        hc_hit.energy = event.cl_energy[0];
        gem_hit.x = event.mHit_gx[0][0];
        gem_hit.y = event.mHit_gy[0][0];
        gem_hit.z = event.mHit_gz[0][0];

        if (gem_hit.z == 0.f) continue;
        const float scale = hc_hit.z / gem_hit.z;
        gem_hit.x *= scale;
        gem_hit.y *= scale;
        gem_hit.z *= scale;
        ApplyToHyCal(gem_hit, run_config);
        ApplyToHyCal(hc_hit, run_config);
        if (!inHyCal(hc_hit.x, hc_hit.y)) continue;

        const auto *mod = hycal.module_by_id(event.cl_center[0]);
        if (!mod) continue;
        const float dx = hc_hit.x - gem_hit.x;
        const float dy = hc_hit.y - gem_hit.y;
        float xd_hycal = (hc_hit.x - mod->x) / mod->size_x;
        float yd_hycal = (hc_hit.y - mod->y) / mod->size_y;
        if (xd_hycal < -0.5f) xd_hycal += 1.0f;
        if (xd_hycal >  0.5f) xd_hycal -= 1.0f;
        if (yd_hycal < -0.5f) yd_hycal += 1.0f;
        if (yd_hycal >  0.5f) yd_hycal -= 1.0f;
        float xd_gem = (gem_hit.x - mod->x) / mod->size_x;
        float yd_gem = (gem_hit.y - mod->y) / mod->size_y;
        if (xd_gem < -0.5f) xd_gem += 1.0f;
        if (xd_gem >  0.5f) xd_gem -= 1.0f;
        if (yd_gem < -0.5f) yd_gem += 1.0f;
        if (yd_gem >  0.5f) yd_gem -= 1.0f;

        float theta = std::atan2(
            std::sqrt(gem_hit.x * gem_hit.x + gem_hit.y * gem_hit.y), gem_hit.z)
            * 180.0f / M_PI;
        if (event.matchNum != 1)
            theta = std::atan2(
                std::sqrt(hc_hit.x * hc_hit.x + hc_hit.y * hc_hit.y), hc_hit.z)
                * 180.0f / M_PI;
        const int theta_bin = result->h1_yield_theta->GetXaxis()->FindBin(theta);
        if (theta_bin >= 1 && theta_bin <= Nbins) {
            result->h1_yield_theta->Fill(
                theta, 1.0 / result->h1_yield_theta->GetXaxis()->GetBinWidth(theta_bin));
            result->h1_Erecon_theta[theta_bin - 1]->Fill(event.cl_energy[0]);
            if (std::fabs(xd_hycal) < 0.3f && std::fabs(yd_hycal) < 0.3f) {
                result->h1_Erecon_theta_center[theta_bin - 1]->Fill(event.cl_energy[0]);
            } else {
                result->h1_Erecon_theta_edge[theta_bin - 1]->Fill(event.cl_energy[0]);
            }
        }

        const int module_index = mod->id - 1001;
        if (module_index < 0 || module_index >= 1156) continue;
        if (std::fabs(xd_hycal) < 0.3f && std::fabs(yd_hycal) < 0.3f)
            result->h1_E_modules[module_index]->Fill(event.cl_energy[0]);
        result->h2_hit_hycal->Fill(hc_hit.x, hc_hit.y);
        result->h2_hit_gem->Fill(gem_hit.x, gem_hit.y);
        result->h2_hit_module_hycal->Fill(xd_hycal, yd_hycal);
        result->h2_hit_module_gem->Fill(xd_gem, yd_gem);
        ++result->events_processed;
    }
    return true;
}