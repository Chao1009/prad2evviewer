//=============================================================================
// Replay.cpp — EVIO to ROOT tree conversion
//=============================================================================

#include "Replay.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

namespace analysis {

void Replay::LoadDaqMap(const std::string &json_path)
{
    std::ifstream f(json_path);
    if (!f.is_open()) {
        std::cerr << "Replay: cannot open DAQ map: " << json_path << "\n";
        return;
    }
    auto j = json::parse(f, nullptr, false, true);
    if (j.is_array()) {
        for (auto &entry : j) {
            std::string name = entry.value("name", "");
            int crate   = entry.value("crate", -1);
            int slot    = entry.value("slot", -1);
            int channel = entry.value("channel", -1);
            if (!name.empty() && crate >= 0)
                daq_map_[std::to_string(crate) + "_" + std::to_string(slot) +
                         "_" + std::to_string(channel)] = name;
        }
    }
    std::cerr << "Replay: loaded " << daq_map_.size() << " DAQ map entries\n";
}

std::string Replay::moduleName(int roc, int slot, int ch) const
{
    auto it = daq_map_.find(std::to_string(roc) + "_" + std::to_string(slot) +
                            "_" + std::to_string(ch));
    return (it != daq_map_.end()) ? it->second : "";
}

int Replay::moduleID(int roc, int slot, int ch) const
{
    auto name = moduleName(roc, slot, ch);
    if (name.empty()) return -1;
    if(name[0] == 'G') return std::stoi(name.substr(1));
    else if(name[0] == 'W') return std::stoi(name.substr(1))+1000;
    else return -1; // Unknown module type
}

float Replay::computeIntegral(const fdec::ChannelData &cd, float pedestal) const
{
    float sum = 0.f;
    for (int i = 0; i < cd.nsamples; ++i)
        sum += cd.samples[i] - pedestal;
    return sum;
}

void Replay::clearEvent(EventVars &ev)
{
    ev.event_num = 0;
    ev.trigger = 0;
    ev.timestamp = 0;
    ev.nch = 0;
    ev.gem_nch = 0;
}

void Replay::setupBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    tree->Branch("event_num", &ev.event_num, "event_num/i");
    tree->Branch("trigger",   &ev.trigger,   "trigger/i");
    tree->Branch("timestamp", &ev.timestamp, "timestamp/L");
    tree->Branch("hycal.nch",       &ev.nch,       "nch/I");
    tree->Branch("hycal.crate",     ev.crate,      "crate[nch]/b");
    tree->Branch("hycal.slot",      ev.slot,       "slot[nch]/b");
    tree->Branch("hycal.channel",   ev.channel,    "channel[nch]/b");
    tree->Branch("hycal.module_id", ev.module_id,  "module_id[nch]/I");
    tree->Branch("hycal.nsamples",  ev.nsamples,   "nsamples[nch]/b");
    tree->Branch("hycal.ped_mean",  ev.ped_mean,   "ped_mean[nch]/F");
    tree->Branch("hycal.ped_rms",   ev.ped_rms,    "ped_rms[nch]/F");
    tree->Branch("hycal.integral",  ev.integral,   "integral[nch]/F");
    if (write_peaks) {
        tree->Branch("hycal.npeaks",       ev.npeaks,       "npeaks[nch]/b");
        tree->Branch("hycal.peak_height",  ev.peak_height,  Form("peak_height[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_time",    ev.peak_time,    Form("peak_time[nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_integral",ev.peak_integral, Form("peak_integral[nch][%d]/F", fdec::MAX_PEAKS));
    }
    //GEM part
    tree->Branch("gem.nch",        &ev.gem_nch,   "gem_nch/I");
    tree->Branch("gem.mpd_crate",  ev.mpd_crate,  "mpd_crate[gem_nch]/b");
    tree->Branch("gem.mpd_fiber",  ev.mpd_fiber,  "mpd_fiber[gem_nch]/b");
    tree->Branch("gem.apv",        ev.apv,        "apv[gem_nch]/b");
    tree->Branch("gem.strip",      ev.strip,      "strip[gem_nch]/b");
    tree->Branch("gem.ssp_samples",ev.ssp_samples,Form("ssp_samples[gem_nch][%d]/F", ssp::SSP_TIME_SAMPLES));
}

bool Replay::Process(const std::string &input_evio, const std::string &output_root,
                     int max_events, bool write_peaks , const std::string &daq_config_file)
{
    // build ROC tag → crate index mapping from DAQ config JSON
    std::unordered_map<int, int> roc_to_crate;
    if (!daq_config_file.empty()) {
        std::cout << "Loading DAQ config from " << daq_config_file << "\n";
        std::ifstream dcf(daq_config_file);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("roc_tags") && dcj["roc_tags"].is_array()) {
                for (auto &entry : dcj["roc_tags"]) {
                    int tag   = std::stoi(entry.at("tag").get<std::string>(), nullptr, 16);
                    int crate = entry.at("crate").get<int>();
                    roc_to_crate[tag] = crate;
                }
            }
        }
    }
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.Open(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("events", "PRad2 replay data");
    //EventVars ev;
    auto ev = std::make_unique<EventVars>();
    setupBranches(tree, *ev, write_peaks);

    fdec::EventData event;
    ssp::SspEventData ssp_evt;
    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    int total = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != evc::EventType::Physics) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event.clear();
            ssp_evt.clear();
            if (!ch.DecodeEvent(ie, event, &ssp_evt)) continue;
            if (max_events > 0 && total >= max_events) break;

            clearEvent(*ev);
            ev->event_num = event.info.event_number;
            ev->trigger   = event.info.trigger_bits;
            ev->timestamp = event.info.timestamp;

            // decode HyCal FADC250 data
            int nch = 0;
            for (int r = 0; r < event.nrocs; ++r) {
                auto &roc = event.rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                int crate;
                if (cit == roc_to_crate.end()) crate = roc.tag;
                else crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 16; ++c) {
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0 || nch >= kMaxCh) continue;

                        ev->crate[nch]   = crate;
                        ev->slot[nch]    = s;
                        ev->channel[nch] = c;
                        //ev->module_id[nch] = moduleID(roc.tag, s, c);
                        ev->nsamples[nch] = cd.nsamples;

                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        ev->ped_mean[nch] = wres.ped.mean;
                        ev->ped_rms[nch]  = wres.ped.rms;
                        ev->integral[nch] = computeIntegral(cd, wres.ped.mean);

                        if (write_peaks) {
                            ev->npeaks[nch] = wres.npeaks;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                ev->peak_height[nch][p]   = wres.peaks[p].height;
                                ev->peak_time[nch][p]     = wres.peaks[p].time;
                                ev->peak_integral[nch][p] = wres.peaks[p].integral;
                            }
                        }
                        nch++;
                    }
                }
            }
            ev->nch = nch;

            // decode GEM SSP data
            int gem_ch = 0;
            for (int m = 0; m < ssp_evt.nmpds; ++m) {
                auto &mpd = ssp_evt.mpds[m];
                if (!mpd.present) continue;
                for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
                    auto &apv = mpd.apvs[a];
                    if (!apv.present) continue;
                    int idx = -1; // find APV index in GemSystem if needed
                    for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                        if (!apv.hasStrip(s)) continue;
                        if (gem_ch >= GEMkMaxCH) continue;
                        
                        ev->mpd_crate[gem_ch] = mpd.crate_id;
                        ev->mpd_fiber[gem_ch] = mpd.mpd_id;
                        ev->apv[gem_ch]       = a;
                        ev->strip[gem_ch]     = s;
                        for (int t = 0; t < ssp::SSP_TIME_SAMPLES; t++)
                            ev->ssp_samples[gem_ch][t] = apv.strips[s][t];

                        gem_ch++;
                    }
                }
            }
            ev->gem_nch = gem_ch; // total channels = HyCal + GEM
            tree->Fill();
            total++;

            if (total % 1000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
        if (max_events > 0 && total >= max_events) break;
    }

    std::cerr << "\rReplay: " << total << " events written to " << output_root << "\n";
    tree->Write();
    delete outfile;
    return true;
}

} // namespace analysis
