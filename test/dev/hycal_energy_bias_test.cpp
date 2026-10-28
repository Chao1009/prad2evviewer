#include "HyCalCluster.h"
#include "HyCalEnergyBias.h"
#include "HyCalSystem.h"
#include "PipelineBuilder.h"
#include "test_util.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;
using nlohmann::json;
using namespace testutil;

namespace {

json make_grid(float center, float upper_left)
{
    json module = json::object();
    for (int row = 0; row < fdec::HyCalEnergyBias::GRID_SIZE; ++row) {
        json values = json::object();
        for (int column = 0; column < fdec::HyCalEnergyBias::GRID_SIZE; ++column)
            values["x" + std::to_string(column)] = 0.f;
        module["y" + std::to_string(row)] = std::move(values);
    }
    module["y2"]["x2"] = center;
    module["y4"]["x0"] = upper_left;
    return module;
}

void write_fixture(const fs::path &path, float center, float upper_left)
{
    json root = {
        {"comment", "test fixture"},
        {"W1", make_grid(center, upper_left)}
    };
    std::ofstream output(path);
    output << root.dump(2) << '\n';
}

} // namespace

int main()
{
    fdec::HyCalSystem hycal;
    check(hycal.Init(std::string(DATABASE_DIR) + "/hycal_map.json"),
          "load HyCal map");
    const auto *w1 = hycal.module_by_name("W1");
    const auto *g1 = hycal.module_by_name("G1");
    check(w1 != nullptr, "find W1");
    check(g1 != nullptr, "find G1");
    if (!w1 || !g1) return 1;

    const fs::path fixture_dir = make_temp_dir("prad2_energy_bias_");
    const fs::path ee_path = fixture_dir / "ee.json";
    const fs::path ep_path = fixture_dir / "ep.json";
    write_fixture(ee_path, 0.10f, 0.25f);
    write_fixture(ep_path, -0.20f, -0.30f);

    auto bias = fdec::LoadHyCalEnergyBias(
        ee_path.string(), ep_path.string(), hycal, 2200.f);
    check(bias->ee_cells_loaded == 25, "load 25 ee cells");
    check(bias->ep_cells_loaded == 25, "load 25 ep cells");

    const float threshold = bias->ep_threshold();
    check(!bias->is_ep(threshold), "threshold equality selects ee");
    check(bias->is_ep(threshold + 0.01f), "energy above threshold selects ep");

    check(close_to(bias->bias(*w1, static_cast<float>(w1->x),
                              static_cast<float>(w1->y), threshold), 0.10f),
          "center ee cell lookup");
    check(close_to(bias->bias(*w1,
                              static_cast<float>(w1->x - 0.49 * w1->size_x),
                              static_cast<float>(w1->y + 0.49 * w1->size_y),
                              threshold), 0.25f),
          "columns increase with +x and rows increase with +y");
    check(close_to(bias->bias(*w1,
                              static_cast<float>(w1->x - 2.0 * w1->size_x),
                              static_cast<float>(w1->y + 2.0 * w1->size_y),
                              threshold), 0.25f),
          "out-of-cell positions clamp to edge bin");
    check(close_to(bias->correction_factor(*g1, static_cast<float>(g1->x),
                                           static_cast<float>(g1->y), 2100.f), 1.f),
          "PbGlass uses identity correction");
    check(close_to(bias->correction_factor(*w1, static_cast<float>(w1->x),
                                           static_cast<float>(w1->y), 2100.f), 1.25f),
          "ep bias converts to reciprocal correction factor");

    check(close_to(fdec::SelectHyCalEnergyBiasSet(728.9f).nominal_mev, 700.f),
          "select nearest 0.7 GeV parameter set");
    check(close_to(fdec::SelectHyCalEnergyBiasSet(2108.8f).nominal_mev, 2200.f),
          "select nearest 2.2 GeV parameter set");
    check(close_to(fdec::SelectHyCalEnergyBiasSet(3488.43f).nominal_mev, 3500.f),
          "select nearest 3.5 GeV parameter set");

    fdec::ClusterConfig config;
    config.min_module_energy = 0.f;
    config.min_center_energy = 10.f;
    config.min_cluster_energy = 50.f;
    config.non_linear_corr = false;
    config.energy_bias_correction = true;
    config.energy_bias = bias;

    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(config);
    clusterer.AddHit(w1->index, 2100.f, 0.f);
    clusterer.FormClusters();
    std::vector<fdec::ClusterHit> hits;
    clusterer.ReconstructHits(hits);
    check(hits.size() == 1, "reconstruct one W1 cluster");
    if (hits.size() == 1) {
          check(close_to(hits[0].linear_corr, 1.f),
              "energy bias does not change non-linearity factor");
          check(close_to(hits[0].bias_corr, 1.25f),
              "cluster reports energy-bias correction separately");
        check(close_to(hits[0].energy, 2625.f, 0.01f),
              "cluster applies energy-bias correction before output");
    }

    config.energy_bias_correction = false;
    clusterer.SetConfig(config);
    clusterer.Clear();
    clusterer.AddHit(w1->index, 2100.f, 0.f);
    clusterer.FormClusters();
    clusterer.ReconstructHits(hits);
        check(hits.size() == 1 && close_to(hits[0].bias_corr, 1.f) &&
            close_to(hits[0].energy, 2100.f),
          "disabled correction preserves cluster energy");

        const fs::path recon_path = fixture_dir / "reconstruction_config.json";
        {
          std::ofstream output(recon_path);
          output << json({
            {"runinfo", "runinfo/general.json"},
            {"hycal", {{"energy_bias_correction", true}}}
          }).dump(2) << '\n';
        }
        auto pipeline = prad2::PipelineBuilder()
          .set_database_dir(DATABASE_DIR)
          .set_recon_config(recon_path.string())
          .set_run_number(24560)
          .set_log_stream(nullptr)
          .build();
        check(close_to(pipeline.run_cfg.Ebeam, 728.9f),
            "PipelineBuilder loads representative 0.7 GeV run");
        check(close_to(pipeline.hycal_energy_bias_nominal, 700.f),
            "PipelineBuilder selects 0.7 GeV energy-bias files");
        check(pipeline.hycal_energy_bias != nullptr &&
            pipeline.hycal_energy_bias->ee_cells_loaded > 0 &&
            pipeline.hycal_energy_bias->ep_cells_loaded > 0,
            "PipelineBuilder loads both ee and ep tables");
        check(pipeline.hycal_cluster_cfg.energy_bias == pipeline.hycal_energy_bias,
            "PipelineBuilder propagates table into ClusterConfig");
        check(pipeline.hycal_energy_bias_ee_path.find("0p7GeV_ee.json") !=
              std::string::npos &&
            pipeline.hycal_energy_bias_ep_path.find("0p7GeV_ep.json") !=
              std::string::npos,
            "PipelineBuilder exposes resolved ee and ep paths");

    fs::remove_all(fixture_dir);
    return finish("energy-bias", "HyCal energy-bias");
}
