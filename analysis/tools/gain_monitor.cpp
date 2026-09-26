//=============================================================================
// gain_monitor.cpp — EVIO -> per-channel LMS / alpha amplitude histograms
//
// Usage: gain_monitor -r <run_number> [-s start_file] [-e end_file]
//                     [-o output_dir] [-i evio_dir]
//
// Reads splits start_file..end_file of <evio_dir>/prad_<run>/ and fills the
// amplitude (mean of +-3 samples around the maximum in samples [40,60],
// minus pedestal; 900 bins over 0-4500) of LMS1-3 on alpha triggers and of
// LMS1-3 plus every W/G module on LMS triggers.  The output
// <output_dir>/prad_<run>_LMS_file_<start>_<end>.root is merged by
// scripts/shell/run_gain_monitor.sh and fitted by prad2ana_gain_fitter, whose
// .dat feeds scripts/hycal_gain_monitor.py.
//=============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"
#include "Fadc250Data.h"
#include "EventData.h"
#include "EvioFiles.h"
#include "InstallPaths.h"

#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

std::string EVIODIR = ".";
std::string OUTPUTDIR = ".";

static constexpr int N_LMS_REF = 3;   // LMS1, LMS2, LMS3

struct ChInfo {
    std::string name;
    bool is_lms    = false;
    int  hist_idx  = -1;   // precomputed index into LMSHist / AlphaHist
};

static std::unordered_map<int, ChInfo> gCh;   // key = crate*10000 + slot*100 + ch
std::vector<TH1F*> ObjectContainer;   // container for saving root histogram/graph

static int packAddr(int c, int s, int ch) { return c * 10000 + s * 100 + ch; }
static std::vector<std::string> loadDaqMap(const char *path);
std::vector<std::string> discoverFiles(const char *dir, const std::string &run, unsigned int startFileNum, unsigned int endFileNum);
TH1F* Init1DHist(const std::string & name,  const std::string & title, const int & nbin,
                 const double & min, const double & max,
                 const std::string & xaxis, const std::string & yaxis, const int & color);
void WritePlotToFile(const std::string fileName);

int main(int argc, char *argv[])
{
    std::string run_number = "";
    unsigned int startFileNum = 0;
    unsigned int endFileNum = 100;
    int opt;
    while ((opt = getopt(argc, argv, "r:s:e:o:i:")) != -1) {
        switch (opt) {
            case 'r': run_number    = optarg; break;
            case 's': startFileNum  = std::atoi(optarg); break;
            case 'e': endFileNum    = std::atoi(optarg); break;
            case 'o': OUTPUTDIR     = optarg; break;
            case 'i': EVIODIR     = optarg; break;
        }
    }
    if (endFileNum < startFileNum) endFileNum = startFileNum;

    if (run_number.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " -r <run_number> [-s start_file] [-e end_file]"
                  << " [-o output_dir] [-i evio_dir]\n";
        return 1;
    }

    TString dbDir = prad2::database_dir().c_str();

    TString cfgFile = Form("%s/daq_config.json", dbDir.Data());
    TString mapFile = Form("%s/hycal_map.json", dbDir.Data());

    printf("============================================\n");
    printf(" PRad LMS / Alpha Normalization\n");
    printf(" Run        : %s\n", run_number.data());
    printf(" Data Dir   : %s\n", EVIODIR.data());
    printf(" Output Dir : %s\n", OUTPUTDIR.data());
    printf(" DAQ cfg    : %s\n", cfgFile.Data());
    printf(" HyCal map  : %s\n", mapFile.Data());
    printf("============================================\n");

    // --- load configs ---
    evc::DaqConfig cfg;
    if (!evc::load_daq_config(cfgFile.Data(), cfg)) {
        std::cerr << "FATAL: failed to load " << cfgFile << "\n";
        return 0;
    }

    std::vector<std::string> modules = loadDaqMap(mapFile.Data());
    if (modules.empty()) {
        std::cerr << "FATAL: no HyCal modules found in daq_map\n";
        return 0;
    }
    printf("HyCal modules: %zu\n", modules.size());

    const auto tagToCrate = cfg.roc_crate_map();

    // LMS histogram order: LMS1-3, then W before G by module number.
    // gain_fitter takes the reference channels by key position.
    std::sort(modules.begin(), modules.end(),
              [](const std::string &a, const std::string &b) {
                  if (a[0] != b[0]) return a[0] == 'W';
                  return std::atoi(a.c_str() + 1) < std::atoi(b.c_str() + 1);
              });
    std::vector<std::string> indexToName = {"LMS1", "LMS2", "LMS3"};
    indexToName.insert(indexToName.end(), modules.begin(), modules.end());

    std::unordered_map<std::string, int> nameToIndex;
    nameToIndex.reserve(indexToName.size());
    for (size_t i = 0; i < indexToName.size(); ++i)
        nameToIndex[indexToName[i]] = static_cast<int>(i);

    // Pre-compute hist_idx in each ChInfo to avoid a second map lookup in the hot loop
    for (auto &[addr, ci] : gCh) {
        auto it = nameToIndex.find(ci.name);
        if (it != nameToIndex.end()) ci.hist_idx = it->second;
    }

    std::vector<TH1F*> LMSHist;
    std::vector<TH1F*> AlphaHist;

    for (int i=0; i<N_LMS_REF; ++i){
        AlphaHist.push_back(Init1DHist(Form("LMS%d_Alpha", i+1),
                                       Form("LMS%d_Alpha", i+1),
                                       900, 0, 4500, "ADC", "count", 1));
    }
    for (const auto &name : indexToName){
        LMSHist.push_back(Init1DHist(Form("%s_LMS", name.c_str()),
                                     Form("%s_LMS", name.c_str()),
                                     900, 0, 4500, "ADC", "count", 1));
    }
    std::vector<std::string> InputFiles = discoverFiles(EVIODIR.data(), run_number, startFileNum, endFileNum);
    printf("found EVIO files : %zu\n", InputFiles.size());

    // --- decoder objects ---
    evc::EvChannel reader;
    reader.SetConfig(cfg);

    static fdec::EventData evt;
    fdec::WaveAnalyzer wave(cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    template_store.LoadFromConfig(cfg.wave_cfg, std::string(dbDir.Data()));
    wave.SetTemplateStore(&template_store);

    // Hoist outside all loops: avoids repeated heap allocation/deallocation per event.
    std::unordered_map<int, float> integrals;
    integrals.reserve(2048);

    for (size_t fi = 0; fi < InputFiles.size(); fi++) {
        printf("[%zu/%zu] %s\n", fi + 1, InputFiles.size(), InputFiles[fi].c_str());

        if (reader.OpenAuto(InputFiles[fi]) != evc::status::success) {
            printf("  WARNING: cannot open, skipping\n");
            continue;
        }

        while (reader.Read() == evc::status::success) {
            if (!reader.Scan()) continue;
            if (reader.GetEventType() != evc::EventType::Physics) continue;

            for (int ie = 0; ie < reader.GetNEvents(); ie++) {

                evt.clear();
                if (!reader.DecodeEvent(ie, evt)) continue;

                bool isLMS = (evt.info.trigger_bits & prad2::TBIT_lms);
                bool isAlpha = (evt.info.trigger_bits & prad2::TBIT_alpha);

                if (!(isLMS || isAlpha)) continue;

                // -- compute waveform integrals for all channels ----
                // store as:  packed_addr -> integral
                integrals.clear();

                for (int r = 0; r < fdec::MAX_ROCS; r++) {
                    auto &roc = evt.rocs[r];
                    if (!roc.present) continue;
                    auto ct = tagToCrate.find(roc.tag);
                    if (ct == tagToCrate.end()) continue;
                    int crate = ct->second;

                    fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                        fdec::WaveResult wres;
                        wave.SetChannelKey(roc.tag, s, c);
                        wave.Analyze(cd.samples, cd.nsamples, wres);

                        // Find peak sample index in window [40, 60]
                        int peakSample = 40;
                        float peakVal = cd.samples[40] - wres.ped.mean;
                        int searchEnd = std::min(60, cd.nsamples - 1);
                        for (int k = 41; k <= searchEnd; ++k) {
                            float v = cd.samples[k] - wres.ped.mean;
                            if (v > peakVal) { peakVal = v; peakSample = k; }
                        }

                        // Average ±3 samples around the peak
                        int avgStart = std::max(0, peakSample - 3);
                        int avgEnd   = std::min(cd.nsamples - 1, peakSample + 3);
                        float sum = 0.f;
                        int   cnt = 0;
                        for (int k = avgStart; k <= avgEnd; ++k) {
                            sum += cd.samples[k] - wres.ped.mean;
                            ++cnt;
                        }
                        float integral = (cnt > 0) ? sum / cnt : 0.f;

                        integrals[packAddr(crate, s, c)] = integral;
                    });
                }

                const size_t nIntegrals = integrals.size();
                for (auto &[addr, integ] : integrals) {
                    auto it = gCh.find(addr);
                    if (it == gCh.end()) continue;
                    const auto &ci = it->second;
                    if (ci.hist_idx < 0) continue;

                    if (ci.is_lms && isAlpha && nIntegrals < 500){
                        assert(ci.hist_idx < N_LMS_REF);
                        AlphaHist[ci.hist_idx]->Fill(integ);
                    }
                    if (isLMS && nIntegrals > 500){
                        assert(ci.hist_idx < (int)LMSHist.size());
                        LMSHist[ci.hist_idx]->Fill(integ);
                    }
                }
            }
        }

        reader.Close();
    }

    WritePlotToFile(Form("%s/prad_%s_LMS_file_%d_%d.root", OUTPUTDIR.c_str(), run_number.c_str(), startFileNum, endFileNum));
    return 0;
}

// Load database/hycal_map.json — maps every (crate,slot,channel) to a name.
// W* / G* names are HyCal modules; LMS1/2/3 are reference channels.
// Records without a "daq" block (boosters, PRad-1 V1-V4) are skipped.
// Returns the HyCal module names.
static std::vector<std::string> loadDaqMap(const char *path)
{
    std::vector<std::string> modules;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open hycal_map " << path << "\n";
        return modules;
    }
    auto j = nlohmann::json::parse(f);

    for (auto &e : j) {
        if (!e.contains("daq")) continue;
        std::string nm = e.value("n", "");
        const auto &d  = e["daq"];
        int cr = d.value("crate", -1);
        int sl = d.value("slot", -1);
        int ch = d.value("channel", -1);

        ChInfo ci;
        ci.name = nm;
        ci.is_lms = (nm == "LMS1" || nm == "LMS2" || nm == "LMS3");
        if (!nm.empty() && (nm[0] == 'W' || nm[0] == 'G')) modules.push_back(nm);

        gCh[packAddr(cr, sl, ch)] = ci;
    }
    return modules;
}

// Splits startFileNum..endFileNum of <dir>/prad_<run>/, in split order.
std::vector<std::string> discoverFiles(const char *dir, const std::string &run, unsigned int startFileNum, unsigned int endFileNum)
{
    fs::path subdir = fs::path(dir) / ("prad_" + run);
    if (!fs::is_directory(subdir)) {
        throw std::runtime_error("Subdirectory not found: " + subdir.string());
    }

    std::vector<std::string> result;
    for (auto &file : prad2::discover_split_files(subdir.string())) {
        unsigned long subrun = std::stoul(file.substr(file.rfind('.') + 1));
        if (subrun >= startFileNum && subrun <= endFileNum) result.push_back(std::move(file));
    }
    return result;
}

TH1F* Init1DHist(const std::string & name,  const std::string & title, const int & nbin, const double & min, const double & max,
                 const std::string & xaxis, const std::string & yaxis, const int & color)
{
    TH1F* h = new TH1F(name.c_str(), title.c_str(), nbin, min, max);
    h->GetXaxis()->SetTitle(xaxis.c_str());
    h->GetYaxis()->SetTitle(yaxis.c_str());
    h->SetLineWidth(2);
    h->SetLineColor(color);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    ObjectContainer.push_back(h);
    return h;
}

void WritePlotToFile(const std::string fileName)
{
    TFile* f = new TFile(fileName.c_str(), "RECREATE");
    f->cd();

    for (unsigned int i = 0; i < ObjectContainer.size(); i++)
    ObjectContainer[i]->Write();

    f->Close();
    delete f;
}
