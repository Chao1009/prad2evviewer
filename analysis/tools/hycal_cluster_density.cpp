#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
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
#include <TF1.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <cstdio>
#include <memory>
#include <mutex>
#include <getopt.h>

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

static std::string outputFileName(const std::string &output_name, bool corr = false)
{
    const fs::path output_path(output_name);
    const std::string file_name = "hycal_cluster_density"
        + output_path.filename().string()
        + (corr ? ".corr" : "") + ".root";
    return (output_path.parent_path() / file_name).string();
}

static std::unique_ptr<TGraphErrors> extractMainPeakCenters(
    const TH2 &histogram, const char *graph_name)
{
    constexpr int x_bins_per_slice = 1;
    constexpr double fit_half_width = 1.5;
    auto centers = std::make_unique<TGraphErrors>();
    centers->SetName(graph_name);

    const TAxis *x_axis = histogram.GetXaxis();
    const TAxis *y_axis = histogram.GetYaxis();
    for (int first_x_bin = 1; first_x_bin <= x_axis->GetNbins();
         first_x_bin += x_bins_per_slice) {
        const int last_x_bin = first_x_bin;
        std::unique_ptr<TH1D> projection(histogram.ProjectionY(
            "slice_peak_fit", first_x_bin, last_x_bin, "e"));
        if (projection->GetEntries() < 100.0) continue;

        const double peak_position = y_axis->GetBinCenter(projection->GetMaximumBin());
        TF1 peak_fit("slice_peak_fit_function", "gaus",
                     peak_position - fit_half_width, peak_position + fit_half_width);
        peak_fit.SetParameters(projection->GetMaximum(), peak_position, 0.8);
        peak_fit.SetParLimits(2, 0.1, fit_half_width);
        if (projection->Fit(&peak_fit, "QNR") != 0) continue;

        const double center = peak_fit.GetParameter(1);
        const double center_error = peak_fit.GetParError(1);
        if (center_error <= 0.0 || !std::isfinite(center)) continue;
        const int point = centers->GetN();
        const double x_center = x_axis->GetBinCenter(first_x_bin);
        centers->SetPoint(point, x_center, center);
        centers->SetPointError(point, 0.0, center_error);
    }
    return centers;
}

int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

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
    while ((opt = getopt_long_only(argc, argv, "o:n:f:j:wc",
                                   long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
            case 'w': worker_mode      = true; break;
            case 'c': corr             = true; break;
        }
    }

    std::vector<std::string> root_files =
        CollectInputs(argc, argv, optind, IsReconRootName, num_files);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_cluster_density <input_recon.root|dir> [more...] "
                 "-o <output_name> [-n max_events] [-f nfiles] [-j threads] [-corr]\n";
        return 1;
    }

    if (output_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }

    // ROOT TTree branch addresses are not safe to share across worker threads:
    // run one isolated worker process per input file, merge histograms here.
    if (!corr && !worker_mode && root_files.size() > 1) {
        std::vector<std::string> worker_outputs(root_files.size());
        std::mutex io_mutex;
        ParallelFor(root_files.size(), num_threads, [&](size_t i, int) {
            const std::string worker_prefix = output_name + ".worker_" + std::to_string(i);
            worker_outputs[i] = outputFileName(worker_prefix);
            std::vector<std::string> args{argv[0], "-w", "-o", worker_prefix};
            if (max_events > 0) args.insert(args.end(), {"-n", std::to_string(max_events)});
            args.insert(args.end(), {"-j", "1", root_files[i]});
            const int rc = RunCommand(args);
            std::lock_guard<std::mutex> lock(io_mutex);
            if (rc != 0)
                std::cerr << "Worker failed for " << root_files[i]
                          << " (exit code " << rc << ")\n";
        });
        if (!MergeTopLevelHistograms(worker_outputs, outputFileName(output_name, corr)))
            return 1;
        for (const auto &worker_output : worker_outputs)
            std::remove(worker_output.c_str());
        return 0;
    }

    std::unique_ptr<TF1> fit_func_x;
    std::unique_ptr<TF1> fit_func_y;
    std::unique_ptr<TGraphErrors> peak_centers_x;
    std::unique_ptr<TGraphErrors> peak_centers_y;
    std::unique_ptr<TH2> fit_source_dx;
    std::unique_ptr<TH2> fit_source_dy;
    if (corr) {
        const std::string source_name = outputFileName(output_name);
        TFile source_file(source_name.c_str(), "READ");
        if (source_file.IsZombie()) {
            std::cerr << "Cannot open uncorrected result " << source_name
                      << ". Run once without -corr first.\n";
            return 1;
        }

        auto *source_dx = dynamic_cast<TH2 *>(source_file.Get("h2_dist_dx_xd"));
        auto *source_dy = dynamic_cast<TH2 *>(source_file.Get("h2_dist_dy_yd"));
        if (!source_dx || !source_dy) {
            std::cerr << "Missing h2_dist_dx_xd or h2_dist_dy_yd in "
                      << source_name << "\n";
            return 1;
        }
        fit_source_dx.reset(dynamic_cast<TH2 *>(source_dx->Clone("fit_source_dx")));
        fit_source_dy.reset(dynamic_cast<TH2 *>(source_dy->Clone("fit_source_dy")));
        if (!fit_source_dx || !fit_source_dy) {
            std::cerr << "Cannot clone two-dimensional residual histograms from "
                      << source_name << "\n";
            return 1;
        }
        fit_source_dx->SetDirectory(nullptr);
        fit_source_dy->SetDirectory(nullptr);

        peak_centers_x = extractMainPeakCenters(*fit_source_dx, "peak_centers_x");
        peak_centers_y = extractMainPeakCenters(*fit_source_dy, "peak_centers_y");
        if (peak_centers_x->GetN() < 4 || peak_centers_y->GetN() < 4) {
            std::cerr << "Cannot extract enough local peak centers for correction fit.\n";
            return 1;
        }
        const char *position_correction =
            "([0]*x + [1]*x^3 + [2]*x^5 + [3]*x^7)"
            "*(x^2 - 0.25)";
        fit_func_x = std::make_unique<TF1>("fit_func_x", position_correction,
                                           -0.5, 0.5);
        fit_func_y = std::make_unique<TF1>("fit_func_y", position_correction,
                                           -0.5, 0.5);
        fit_func_x->SetParNames("c0", "c1", "c2", "c3");
        fit_func_y->SetParNames("c0", "c1", "c2", "c3");
        peak_centers_x->Fit(fit_func_x.get(), "QR0");
        peak_centers_y->Fit(fit_func_y.get(), "QR0");
    }

    TChain tree("recon");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }
    EventVars_Recon ev;
    prad2::SetReconReadBranches(&tree, ev);

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    // Histograms
    TH2F *h2_hit_hycal = new TH2F("h2_hit_hycal", "HyCal Hit Distribution;(X_{hycal}-X_{cell center})/d_{cell size};(Y_{hycal}-Y_{cell center})/d_{cell size}", 100, -0.5, 0.5, 100, -0.5, 0.5);
    TH2F *h2_hit_gem = new TH2F("h2_hit_gem", "GEM Hit Distribution;(X_{gem}-X_{cell center})/d_{cell size};(Y_{gem}-Y_{cell center})/d_{cell size}", 100, -0.5, 0.5, 100, -0.5, 0.5);
    TH2F *h2_dist_dx_xd = new TH2F("h2_dist_dx_xd", "dx vs xd_hycal;relative x to cell center;x_hycal-x_gem [mm]", 100, -0.5, 0.5, 200, -8.0, 8.0);
    TH2F *h2_dist_dy_yd = new TH2F("h2_dist_dy_yd", "dy vs yd_hycal;relative y to cell center;y_hycal-y_gem [mm]", 100, -0.5, 0.5, 200, -8.0, 8.0);
    TH2F *h2_dist_dx_xd_gem = new TH2F("h2_dist_dx_xd_gem", "dx vs xd_gem;relative x to cell center;x_gem-x_hycal [mm]", 100, -0.5, 0.5, 200, -8.0, 8.0);
    TH2F *h2_dist_dy_yd_gem = new TH2F("h2_dist_dy_yd_gem", "dy vs yd_gem;relative y to cell center;y_gem-y_hycal [mm]", 100, -0.5, 0.5, 200, -8.0, 8.0);
    TH1F *h1_density_hycal_xd = new TH1F("h1_density_hycal_xd", "Density of HyCal xd;X_{hycal}-X_{cell center}/d_{cell size};Counts", 100, -0.5, 0.5);
    TH1F *h1_density_hycal_yd = new TH1F("h1_density_hycal_yd", "Density of HyCal yd;Y_{hycal}-Y_{cell center}/d_{cell size};Counts", 100, -0.5, 0.5);
    TH1F *h1_density_gem_xd = new TH1F("h1_density_gem_xd", "Density of GEM xd;X_{gem}-X_{cell center}/d_{cell size};Counts", 100, -0.5, 0.5);
    TH1F *h1_density_gem_yd = new TH1F("h1_density_gem_yd", "Density of GEM yd;Y_{gem}-Y_{cell center}/d_{cell size};Counts", 100, -0.5, 0.5);
    TH2F *h2_npos_xd_hycal = new TH2F("h2_npos_xd_hycal", "Number of blocks to recon pos vs xd;relative x to cell center;Number of blocks", 100, -0.5, 0.5, 9, 0.5, 9.5);
    TH2F *h2_npos_yd_hycal = new TH2F("h2_npos_yd_hycal", "Number of blocks to recon pos vs yd;relative y to cell center;Number of blocks", 100, -0.5, 0.5, 9, 0.5, 9.5);

    // check the residual between HyCal and GEM hits
    TH1F *h1_residual_dx = new TH1F("h1_residual_dx", "Residual in x between HyCal and GEM hits;dx [mm];Counts", 800, -8.0, 8.0);
    TH1F *h1_residual_dy = new TH1F("h1_residual_dy", "Residual in y between HyCal and GEM hits;dy [mm];Counts", 800, -8.0, 8.0);
    TH2F *h2_residual_dx_dy = new TH2F("h2_residual_dx_dy", "Residuals in x vs y between HyCal and GEM hits;dx [mm];dy [mm]", 800, -8.0, 8.0, 800, -8.0, 8.0);

    Long64_t n = tree.GetEntries();
    if (max_events >= 0 && max_events < n) n = max_events;

    for (Long64_t i = 0; i < n; i++) {
        tree.GetEntry(i);
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << n << " entries.\r" << std::flush;

        // trigger selection
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_sum) continue;

        //Event selection, single cluster e-p events, no "kSplit" flag
        if (ev.n_clusters != 1 || ev.matchNum != 1) continue;
        if (ev.cl_nblocks[0] < 2) continue;
        if (fdec::test_bit(ev.cl_flag[0], fdec::kSplit)) continue;
        if (std::fabs(ev.cl_energy[0] - gRunConfig.Ebeam) > 3.0f * 0.03f * std::sqrt(gRunConfig.Ebeam * 1000.f)) continue;

        HCHit hc_hit;
        GEMHit g_hit;

        hc_hit.x = ev.cl_x[0];
        hc_hit.y = ev.cl_y[0];
        hc_hit.z = ev.cl_z[0];
        hc_hit.energy = ev.cl_energy[0];

        // [0][0]: the matched hit of the downstream GEM pair (GEM1/GEM2)
        g_hit.x = ev.mHit_gx[0][0];
        g_hit.y = ev.mHit_gy[0][0];
        g_hit.z = ev.mHit_gz[0][0];
        GetProjection(g_hit, hc_hit.z);

        ApplyToHyCal(g_hit, gRunConfig);
        ApplyToHyCal(hc_hit, gRunConfig);

        const auto &mod = hycal.module_by_id(ev.cl_center[0]);
        if (!InHyCalRing(hc_hit.x, hc_hit.y, 3.0, 15.)) continue;
        if (!mod) continue;

        if (corr) {
            const auto [xd_hycal, yd_hycal] = mod->cell_offset(hc_hit.x, hc_hit.y, true);
            float corr_x = fit_func_x->Eval(xd_hycal);
            float corr_y = fit_func_y->Eval(yd_hycal);
            hc_hit.x -= corr_x;
            hc_hit.y -= corr_y;
        }

        float dx = hc_hit.x - g_hit.x;
        float dy = hc_hit.y - g_hit.y;

        h1_residual_dx->Fill(dx);
        h1_residual_dy->Fill(dy);
        h2_residual_dx_dy->Fill(dx, dy);
        
        const auto [xd_hycal, yd_hycal] = mod->cell_offset(hc_hit.x, hc_hit.y, true);
        const auto [xd_gem, yd_gem] = mod->cell_offset(g_hit.x, g_hit.y, true);

        h2_hit_hycal->Fill(xd_hycal, yd_hycal);
        h2_hit_gem->Fill(xd_gem, yd_gem);
        h2_dist_dx_xd->Fill(xd_hycal, dx);
        h2_dist_dy_yd->Fill(yd_hycal, dy);
        h1_density_hycal_xd->Fill(xd_hycal);
        h1_density_hycal_yd->Fill(yd_hycal);
        h1_density_gem_xd->Fill(xd_gem);
        h1_density_gem_yd->Fill(yd_gem);
        h2_npos_xd_hycal->Fill(xd_hycal, ev.cl_npos[0]);
        h2_npos_yd_hycal->Fill(yd_hycal, ev.cl_npos[0]);
        h2_dist_dx_xd_gem->Fill(xd_gem, dx);
        h2_dist_dy_yd_gem->Fill(yd_gem, dy);

    }

    const std::string output_file_name = outputFileName(output_name, corr);
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Cannot create output file " << output_file_name << "\n";
        return 1;
    }
    h2_hit_hycal->Write();
    h2_hit_gem->Write();
    h2_dist_dx_xd->Write();
    h2_dist_dy_yd->Write();
    h2_dist_dx_xd_gem->Write();
    h2_dist_dy_yd_gem->Write();
    h1_density_hycal_xd->Write();
    h1_density_hycal_yd->Write();
    h1_density_gem_xd->Write();
    h1_density_gem_yd->Write();
    h2_npos_xd_hycal->Write();
    h2_npos_yd_hycal->Write();
    h1_residual_dx->Write();
    h1_residual_dy->Write();
    h2_residual_dx_dy->Write();

    if (corr) {
        fit_func_x->SetLineColor(kRed);
        fit_func_x->SetLineWidth(2);
        fit_func_y->SetLineColor(kRed);
        fit_func_y->SetLineWidth(2);

        TCanvas c_pos_fit("c_pos_fit", "Position Correction Fits (X & Y)",
                          1400, 600);
        c_pos_fit.Divide(2, 1);
        c_pos_fit.cd(1);
        gPad->SetGrid();
        fit_source_dx->Draw("COLZ");
        fit_func_x->Draw("SAME");
        TLegend leg_x(0.55, 0.72, 0.88, 0.88);
        leg_x.SetFillStyle(0);
        leg_x.AddEntry(fit_source_dx.get(), "2D distribution", "f");
        leg_x.AddEntry(fit_func_x.get(), "Peak-center polynomial fit", "l");
        leg_x.Draw();

        c_pos_fit.cd(2);
        gPad->SetGrid();
        fit_source_dy->Draw("COLZ");
        fit_func_y->Draw("SAME");
        TLegend leg_y(0.55, 0.72, 0.88, 0.88);
        leg_y.SetFillStyle(0);
        leg_y.AddEntry(fit_source_dy.get(), "2D distribution", "f");
        leg_y.AddEntry(fit_func_y.get(), "Peak-center polynomial fit", "l");
        leg_y.Draw();

        output_file.cd();
        fit_func_x->Write("fit_func_x");
        fit_func_y->Write("fit_func_y");
        peak_centers_x->Write();
        peak_centers_y->Write();
        c_pos_fit.Write();
    }

    output_file.Close();

}