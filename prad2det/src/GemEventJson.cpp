#include "GemEventJson.h"
#include "GemSystem.h"
#include "SspData.h"

#include <cmath>
#include <string>

using json = nlohmann::json;

namespace gem
{

namespace
{

// Display-only cross-talk flag: a ZS survivor whose peak stays below the
// cross-talk threshold.  StripHit::cross_talk (GemSystem::collectHits)
// uses the one-sided test charge > noise * zs instead, so the two differ.
bool isCrossTalkCandidate(float charge, float noise, float zs_thres, float xt_thres)
{
    return (charge < noise * xt_thres) && (charge > noise * zs_thres);
}

} // namespace

json ZsApvsToJson(const GemSystem &sys, bool round, const ssp::SspEventData *evt)
{
    auto r1 = [round](float v) -> double { return round ? std::round(v * 10.) / 10. : v; };
    const float zs_thres = sys.GetZeroSupThreshold();
    const float xt_thres = sys.GetCrossTalkThreshold();

    json out = json::array();
    auto add_apv = [&](int idx) {
        if (idx < 0 || !sys.HasApvZsHits(idx)) return;
        const auto &cfg = sys.GetApvConfig(idx);

        json channels = json::object();
        for (int ch = 0; ch < ssp::APV_STRIP_SIZE; ++ch) {
            if (!sys.IsChannelHit(idx, ch)) continue;

            float max_charge = -1e9f;
            short max_tb = 0;
            json ts = json::array();
            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
                const float val = sys.GetProcessedAdc(idx, ch, t);
                ts.push_back(r1(val));
                if (val > max_charge) { max_charge = val; max_tb = static_cast<short>(t); }
            }
            const bool xtalk = isCrossTalkCandidate(max_charge, cfg.pedestal[ch].noise,
                                                    zs_thres, xt_thres);
            channels[std::to_string(ch)] = {
                {"charge", r1(max_charge)}, {"max_timebin", max_tb},
                {"cross_talk", xtalk}, {"ts_adc", ts}
            };
        }
        out.push_back({
            {"crate", cfg.crate_id}, {"mpd", cfg.mpd_id}, {"adc", cfg.adc_ch},
            {"channels", channels}
        });
    };

    if (evt)
        evt->forEachApv([&](const ssp::MpdData &mpd, int a, const ssp::ApvData &) {
            add_apv(sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a));
        });
    else
        for (int idx = 0; idx < sys.GetNApvs(); ++idx) add_apv(idx);
    return out;
}

json DetectorsToJson(const GemSystem &sys, bool round)
{
    auto r1 = [round](float v) -> double { return round ? std::round(v * 10.) / 10. : v; };
    auto r2 = [round](float v) -> double { return round ? std::round(v * 100.) / 100. : v; };

    json out = json::array();
    const auto &dets = sys.GetDetectors();
    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        const auto &det = dets[d];
        json dj;
        dj["id"]       = d;
        dj["name"]     = det.name;
        dj["x_pitch"]  = det.planes[0].pitch;
        dj["y_pitch"]  = det.planes[1].pitch;
        dj["x_strips"] = det.planes[0].n_apvs * ssp::APV_STRIP_SIZE;
        dj["y_strips"] = det.planes[1].n_apvs * ssp::APV_STRIP_SIZE;

        for (int p = 0; p < 2; ++p) {
            json cl_arr = json::array();
            for (const auto &cl : sys.GetPlaneClusters(d, p)) {
                json strips = json::array();
                for (const auto &sh : cl.hits) strips.push_back(sh.strip);
                json cj;
                cj["position"]     = r2(cl.position);
                cj["peak_charge"]  = r1(cl.peak_charge);
                cj["total_charge"] = r1(cl.total_charge);
                cj["max_timebin"]  = cl.max_timebin;
                cj["cross_talk"]   = cl.cross_talk;
                cj["size"]         = static_cast<int>(cl.hits.size());
                cj["hit_strips"]   = strips;
                cl_arr.push_back(cj);
            }
            dj[p == 0 ? "x_clusters" : "y_clusters"] = cl_arr;
        }

        json h2d_arr = json::array();
        for (const auto &h : sys.GetHits(d)) {
            json hj;
            hj["x"]        = r2(h.x);
            hj["y"]        = r2(h.y);
            hj["x_charge"] = r1(h.x_charge);
            hj["y_charge"] = r1(h.y_charge);
            hj["x_peak"]   = r1(h.x_peak);
            hj["y_peak"]   = r1(h.y_peak);
            hj["x_size"]   = h.x_size;
            hj["y_size"]   = h.y_size;
            h2d_arr.push_back(hj);
        }
        dj["hits_2d"] = h2d_arr;

        out.push_back(dj);
    }
    return out;
}

} // namespace gem
