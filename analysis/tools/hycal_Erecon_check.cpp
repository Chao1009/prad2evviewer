
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2Poly.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>

#include <iostream>
#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <limits>
#include <getopt.h>

using namespace analysis;

using EventVars_Recon = prad2::ReconEventData;

// scattering-angle bin edges (deg)
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
    std::unique_ptr<TH2F> h2_Erecon_theta;
    std::unique_ptr<TH2F> h2_Erecon_theta_center;
    std::unique_ptr<TH2F> h2_Erecon_theta_edge;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta_center;
    std::vector<std::unique_ptr<TH1F>> h1_Erecon_theta_edge;
    std::vector<std::unique_ptr<TH1F>> h1_E_modules;
    HistList all;
    long long events_processed = 0;
};

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix);
static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result);

static void mergeHistResult(HistResult &destination, const HistResult &source)
{
    AddAll(destination.all, source.all);
    destination.events_processed += source.events_processed;
}

int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_name;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;

    // getopt_long_only reports an unknown word option such as '-xyz' as one
    // unrecognized option instead of splitting it into short options.
    static option long_options[] = {
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

    std::vector<std::string> root_files =
        CollectInputs(argc, argv, optind, IsReconRootName, num_files);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_Erecon_check <input_recon.root|dir> [more...] "
                 "-o <output_name> [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (output_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }
    InitRootThreading();

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    const auto file_limits = DistributeEventBudget(root_files, "recon", max_events);

    auto merged = makeHistResult("");
    std::vector<std::unique_ptr<HistResult>> results(root_files.size());
    RunFilesInRounds(root_files, num_threads,
        [&](int idx, int slot) {
            auto result = makeHistResult(Form("worker%d", slot));
            const bool ok = processRootFile(root_files[idx], gRunConfig, db_dir,
                                            file_limits[idx], result.get());
            results[idx] = std::move(result);
            return ok;
        },
        [&](int first, int last) {
            for (int i = first; i < last; ++i) {
                if (!results[i]) continue;
                mergeHistResult(*merged, *results[i]);
                results[i].reset();
            }
        });

    // Fit the merged histograms, fill new histograms with the fit results
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    // Both maps get the same polygon bins, so one bin table serves both.
    std::vector<int> map_bins;
    auto h1_deltaE_module = std::make_unique<TH1F>("h1_deltaE_module", "Delta E (Module);#Delta E [MeV];Modules",
        1000, -50, 50);
    std::unique_ptr<TH2Poly> h2_deltaE_module_map(analysis::PhysicsTools::MakeModuleMap(
        hycal, "h2_deltaE_module_map", "Delta E (Module) Map;X [mm];Y [mm];#Delta E [MeV]",
        360., map_bins));
    auto h1_resolution_module = std::make_unique<TH1F>("h1_resolution_module", "Resolution (Module);#sigma_{E}/E*#sqrt{E} [%];Modules",
        50, 2.5, 3.5);
    std::unique_ptr<TH2Poly> h2_resolution_module_map(analysis::PhysicsTools::MakeModuleMap(
        hycal, "h2_resolution_module_map", "Resolution (Module) Map;X [mm];Y [mm];#sigma_{E}/E*#sqrt{E} [%]",
        360., map_bins));
    analysis::PhysicsTools physics(hycal);

    h2_deltaE_module_map->GetZaxis()->SetTitle("#Delta E [MeV]");
    h2_deltaE_module_map->SetMinimum(-10.);
    h2_deltaE_module_map->SetMaximum(10.);
    h2_deltaE_module_map->SetOption("colz");
    h2_resolution_module_map->GetZaxis()->SetTitle("#sigma_{E}/E*#sqrt{E} [%]");
    h2_resolution_module_map->SetMinimum(2.5);
    h2_resolution_module_map->SetMaximum(3.5);
    h2_resolution_module_map->SetOption("colz text");

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
            && sigma < 1.5f * expected_sigma;
        if (fit_good) {
            const float delta_energy = peak - expected_energy;
            float resolution = sigma / peak * std::sqrt(peak / 1000.f) * 100.f;
            h1_deltaE_module->Fill(delta_energy);
            h1_resolution_module->Fill(resolution);
            resolution = std::round(resolution * 100.f) / 100.f;

            // Each module has its own polygon bin; assign rather than accumulate.
            const int polygon_bin = map_bins[mod->index];
            if (polygon_bin > 0) {
                h2_deltaE_module_map->SetBinContent(polygon_bin, delta_energy);
                h2_resolution_module_map->SetBinContent(polygon_bin, resolution);
            }
        }
    }
    // Per angle bin fits of all clusters (0), center hits (1) and edge hits (2).
    struct FitSeries {
        std::vector<double> angle, energy, energy_error, resolution, resolution_error;
    };
    std::array<FitSeries, 3> series;
    std::vector<double> expected_angle_points;
    std::vector<double> expect_energy_points;
    for (int i = 0; i < Nbins; ++i) {
        const double angle = 0.5 * (binEdge[i] + binEdge[i + 1]);
        const double expected_energy = analysis::PhysicsTools::ExpectedEnergy(
            angle, gRunConfig.Ebeam, "ep");
        TH1F *const hists[3] = {merged->h1_Erecon_theta[i].get(),
                                merged->h1_Erecon_theta_center[i].get(),
                                merged->h1_Erecon_theta_edge[i].get()};
        bool has_fit = false;
        for (int k = 0; k < 3; ++k) {
            if (hists[k]->GetEntries() < 200) continue;
            const auto fit = physics.fitPeak(hists[k], static_cast<float>(expected_energy), true);
            if (fit[0] <= 0. || fit[1] <= 0.) continue;
            FitSeries &s = series[k];
            s.angle.push_back(angle);
            s.energy.push_back(fit[0]);
            s.energy_error.push_back(fit[3]);
            s.resolution.push_back(100.0 * fit[1] / fit[0] * std::sqrt(fit[0] / 1000.0));
            s.resolution_error.push_back(100.0 * fit[4] / fit[0] * std::sqrt(fit[0] / 1000.0));
            has_fit = true;
        }
        if (has_fit) {
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

    const auto sized_graph = [&](int k) {
        return TGraphErrors(static_cast<int>(series[k].angle.size()));
    };
    TGraphErrors gErecon_vs_theta[3] = {sized_graph(0), sized_graph(1), sized_graph(2)};
    TGraphErrors gResolution_vs_theta[3] = {sized_graph(0), sized_graph(1), sized_graph(2)};
    TGraph gExpectedE_vs_theta(static_cast<int>(expected_angle_points.size()));
    for (size_t i = 0; i < expected_angle_points.size(); ++i) {
        gExpectedE_vs_theta.SetPoint(i, expected_angle_points[i], expect_energy_points[i]);
    }
    for (int k = 0; k < 3; ++k) {
        const FitSeries &s = series[k];
        for (size_t i = 0; i < s.angle.size(); ++i) {
            gErecon_vs_theta[k].SetPoint(i, s.angle[i], s.energy[i]);
            gErecon_vs_theta[k].SetPointError(i, 0.0, s.energy_error[i]);
            gResolution_vs_theta[k].SetPoint(i, s.angle[i], s.resolution[i]);
            gResolution_vs_theta[k].SetPointError(i, 0.0, s.resolution_error[i]);
        }
    }

    // Extend [lo, hi] over values +- errors padded by 8% of their span, or
    // over [0, 1] when they span nothing.
    const auto include_range = [](double &lo, double &hi,
                                  const std::vector<double> &values,
                                  const std::vector<double> &errors) {
        double min_value = std::numeric_limits<double>::max();
        double max_value = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < values.size(); ++i) {
            const double error = i < errors.size() ? std::abs(errors[i]) : 0.0;
            min_value = std::min(min_value, values[i] - error);
            max_value = std::max(max_value, values[i] + error);
        }
        double range_min = 0.0, range_max = 1.0;
        if (max_value > min_value) {
            const double padding = 0.08 * (max_value - min_value);
            range_min = min_value - padding;
            range_max = max_value + padding;
        }
        lo = std::min(lo, range_min);
        hi = std::max(hi, range_max);
    };

    double energy_min_value = std::numeric_limits<double>::max();
    double energy_max_value = std::numeric_limits<double>::lowest();
    double resolution_min_value = std::numeric_limits<double>::max();
    double resolution_max_value = std::numeric_limits<double>::lowest();
    for (const FitSeries &s : series) {
        include_range(energy_min_value, energy_max_value, s.energy, s.energy_error);
        include_range(resolution_min_value, resolution_max_value, s.resolution, s.resolution_error);
    }
    include_range(energy_min_value, energy_max_value, expect_energy_points, {});
    if (!(energy_max_value > energy_min_value)) {
        energy_min_value = 0.0;
        energy_max_value = 1.0;
    }
    if (!(resolution_max_value > resolution_min_value)) {
        resolution_min_value = 0.0;
        resolution_max_value = 1.0;
    }

    struct SeriesStyle {
        Style_t marker;
        Color_t color;
        const char *label;
        const char *legend_option;
    };
    const SeriesStyle series_style[3] = {
        {20, kBlue + 1, "All", "ep"},
        {21, kGreen + 2, "Center", "p"},
        {22, kOrange + 1, "Edge", "p"}};
    // The first graph carries the frame, so it is drawn with the axes.
    const auto draw_series = [&](TGraphErrors *graphs) {
        for (int k = 0; k < 3; ++k) {
            graphs[k].SetMarkerStyle(series_style[k].marker);
            graphs[k].SetMarkerColor(series_style[k].color);
            graphs[k].SetLineColor(series_style[k].color);
            graphs[k].Draw(k == 0 ? "AP" : "P SAME");
        }
    };

    upper_pad.cd();
    TGraphErrors &gErecon_frame = gErecon_vs_theta[0];
    gErecon_frame.SetTitle("Reconstructed energy vs angle; ;E_{recon} [MeV]");
    gErecon_frame.GetXaxis()->SetLimits(binEdge[0], binEdge[Nbins]);
    gErecon_frame.GetXaxis()->SetLabelSize(0.0);
    gErecon_frame.GetXaxis()->SetTitleSize(0.0);
    gErecon_frame.GetYaxis()->SetLabelSize(0.045);
    gErecon_frame.GetYaxis()->SetTitleSize(0.045);
    gErecon_frame.GetYaxis()->SetTitleOffset(0.81);
    gErecon_frame.GetYaxis()->CenterTitle();
    gErecon_frame.SetMinimum(energy_min_value);
    gErecon_frame.SetMaximum(energy_max_value);
    draw_series(gErecon_vs_theta);
    gExpectedE_vs_theta.SetLineColor(kRed + 1);
    gExpectedE_vs_theta.SetLineWidth(2);
    gExpectedE_vs_theta.Draw("L SAME");
    TLegend energy_legend(0.62, 0.78, 0.93, 0.92);
    for (int k = 0; k < 3; ++k)
        energy_legend.AddEntry(&gErecon_vs_theta[k], series_style[k].label,
                               series_style[k].legend_option);
    energy_legend.AddEntry(&gExpectedE_vs_theta, "Expected E", "l");
    energy_legend.Draw();

    lower_pad.cd();
    TGraphErrors &gResolution_frame = gResolution_vs_theta[0];
    gResolution_frame.SetTitle("Energy resolution vs angle;Scattering angle [deg];#sigma / E * #sqrt{E[GeV]} [%]");
    gResolution_frame.GetXaxis()->SetLimits(binEdge[0], binEdge[Nbins]);
    gResolution_frame.GetXaxis()->SetLabelSize(0.060);
    gResolution_frame.GetXaxis()->SetTitleSize(0.070);
    gResolution_frame.GetXaxis()->SetTitleOffset(1.05);
    gResolution_frame.GetYaxis()->SetLabelSize(0.08);
    gResolution_frame.GetYaxis()->SetTitleSize(0.08);
    gResolution_frame.GetYaxis()->SetTitleOffset(0.58);
    gResolution_frame.GetYaxis()->CenterTitle();
    gResolution_frame.SetMinimum(resolution_min_value);
    gResolution_frame.SetMaximum(resolution_max_value);
    draw_series(gResolution_vs_theta);
    c_reconE.Update();

    const std::string output_file_name = output_name + ".root";
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
    h1_resolution_module->Write();
    h2_resolution_module_map->Write();
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
    merged->h2_Erecon_theta->Write();
    merged->h2_Erecon_theta_center->Write();
    merged->h2_Erecon_theta_edge->Write();
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
    result->h2_hit_module_hycal = Book<TH2F>(result->all,
        Form("h2_hit_module_hycal%s", name_suffix.c_str()),
        "HyCal Hit Distribution (Module);(X_{hycal}-X_{cell center})/d_{cell size};(Y_{hycal}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    result->h2_hit_module_gem = Book<TH2F>(result->all,
        Form("h2_hit_module_gem%s", name_suffix.c_str()),
        "GEM Hit Distribution (Module);(X_{gem}-X_{cell center})/d_{cell size};(Y_{gem}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    result->h2_hit_hycal = Book<TH2F>(result->all,
        Form("h2_hit_hycal%s", name_suffix.c_str()),
        "HyCal Hit Distribution;X [mm];Y [mm]", 720, -360, 360, 720, -360, 360);
    result->h2_hit_gem = Book<TH2F>(result->all,
        Form("h2_hit_gem%s", name_suffix.c_str()),
        "GEM Hit Distribution;X [mm];Y [mm]", 720, -360, 360, 720, -360, 360);
    result->h1_yield_theta = Book<TH1F>(result->all,
        Form("h1_yield_theta%s", name_suffix.c_str()),
        "Yield vs Theta;#theta [deg];Yield/binWidth", Nbins, binEdge);

    result->h2_Erecon_theta = Book<TH2F>(result->all,
        Form("h2_Erecon_theta%s", name_suffix.c_str()),
        "Reconstructed Energy vs Scattering Angles;#theta [deg];E_{recon} [MeV]", 320, 0, 8, 10000, 0, 5000);
    result->h2_Erecon_theta_center = Book<TH2F>(result->all,
        Form("h2_Erecon_theta_center%s", name_suffix.c_str()),
        "Reconstructed Energy vs Scattering Angles (Center);#theta [deg];E_{recon} [MeV]", 320, 0, 8, 10000, 0, 5000);
    result->h2_Erecon_theta_edge = Book<TH2F>(result->all,
        Form("h2_Erecon_theta_edge%s", name_suffix.c_str()),
        "Reconstructed Energy vs Scattering Angles (Edge);#theta [deg];E_{recon} [MeV]", 320, 0, 8, 10000, 0, 5000);

    result->h1_Erecon_theta.reserve(Nbins);
    result->h1_Erecon_theta_center.reserve(Nbins);
    result->h1_Erecon_theta_edge.reserve(Nbins);
    for (int i = 0; i < Nbins; ++i) {
        result->h1_Erecon_theta.push_back(Book<TH1F>(result->all,
            Form("h1_Erecon_theta_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg];E_{recon} [MeV];Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
        result->h1_Erecon_theta_center.push_back(Book<TH1F>(result->all,
            Form("h1_Erecon_theta_center_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg] (Center);E_{recon} [MeV]; Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
        result->h1_Erecon_theta_edge.push_back(Book<TH1F>(result->all,
            Form("h1_Erecon_theta_edge_%d%s", i, name_suffix.c_str()),
            Form("Reconstructed Energy in angles [%.3f-%.3f deg] (Edge);E_{recon} [MeV]; Counts", binEdge[i], binEdge[i+1]),
            energy_bins, energy_min, energy_max));
    }

    result->h1_E_modules.reserve(1156);
    for (int i = 0; i < 1156; ++i) {
        const int mod_id = i + 1001;
        result->h1_E_modules.push_back(Book<TH1F>(result->all,
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
        if (event.cl_nblocks[0] < 3) continue;

        HCHit hc_hit;
        GEMHit gem_hit;
        hc_hit.x = event.cl_x[0];
        hc_hit.y = event.cl_y[0];
        hc_hit.z = event.cl_z[0];
        hc_hit.energy = event.cl_energy[0];
        gem_hit.x = event.mHit_gx[0][0];
        gem_hit.y = event.mHit_gy[0][0];

        ApplyToHyCal(hc_hit, run_config);
        if (!InHyCalRing(hc_hit.x, hc_hit.y, 2.5, 15.)) continue;

        const auto *mod = hycal.module_by_id(event.cl_center[0]);
        if (!mod) continue;
        const auto [xd_hycal, yd_hycal] = mod->cell_offset(hc_hit.x, hc_hit.y, true);
        const auto [xd_gem, yd_gem] = mod->cell_offset(gem_hit.x, gem_hit.y, true);

        const float theta = std::atan2(
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

        result->h2_Erecon_theta->Fill(theta, event.cl_energy[0]);
        if (std::fabs(xd_hycal) < 0.3f && std::fabs(yd_hycal) < 0.3f) {
            result->h2_Erecon_theta_center->Fill(theta, event.cl_energy[0]);
        } else {
            result->h2_Erecon_theta_edge->Fill(theta, event.cl_energy[0]);
        }

        const int module_index = mod->id - 1001;
        if (module_index < 0 || module_index >= 1156) continue;
        result->h1_E_modules[module_index]->Fill(event.cl_energy[0]);
        result->h2_hit_hycal->Fill(hc_hit.x, hc_hit.y);
        result->h2_hit_gem->Fill(gem_hit.x, gem_hit.y);
        result->h2_hit_module_hycal->Fill(xd_hycal, yd_hycal);
        result->h2_hit_module_gem->Fill(xd_gem, yd_gem);
        ++result->events_processed;
    }
    return true;
}