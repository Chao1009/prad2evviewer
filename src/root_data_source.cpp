#ifdef WITH_ROOT

#include "root_data_source.h"
#include "EventData_io.h"
#include "EvioFiles.h"

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <algorithm>

std::unique_ptr<DataSource> createRootDataSource(
    const std::string &path,
    const std::unordered_map<int, uint32_t> &crate_to_roc,
    const fdec::HyCalSystem *hycal)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return nullptr;

    if (f->Get<TTree>("events"))
        return std::make_unique<RootRawDataSource>(crate_to_roc, hycal);
    if (f->Get<TTree>("recon"))
        return std::make_unique<RootReconDataSource>();

    return nullptr;
}

std::string RootRawDataSource::open(const std::string &path)
{
    close();
    file_.reset(TFile::Open(path.c_str(), "READ"));
    if (!file_ || file_->IsZombie()) {
        file_.reset();
        return "cannot open ROOT file";
    }
    tree_ = file_->Get<TTree>("events");
    if (!tree_) { close(); return "no 'events' tree in ROOT file"; }

    n_entries_ = static_cast<int>(tree_->GetEntries());

    auto status = prad2::SetRawReadBranches(tree_, ev_);
    has_peaks_ = status.has_peaks;
    has_gem_   = status.has_gem;

    if (!hycal_) {
        std::cerr << "ROOT raw: warning — no HyCalSystem provided; "
                     "module_id → DAQ mapping unavailable, HyCal channels will be skipped\n";
    }

    std::cerr << "ROOT raw: " << n_entries_ << " events"
              << (has_peaks_ ? ", peaks" : "")
              << (has_gem_   ? ", GEM"   : "") << "\n";
    return "";
}

void RootRawDataSource::close()
{
    tree_ = nullptr;
    file_.reset();
    n_entries_ = 0;
    has_peaks_ = false;
    has_gem_ = false;
}

DataSourceCaps RootRawDataSource::capabilities() const
{
    return {
        true,       // has_waveforms
        has_peaks_, // has_peaks
        true,       // has_pedestals
        true,       // has_clusters (computed)
        has_gem_,   // has_gem_raw
        has_gem_,   // has_gem_hits (computed)
        false,      // has_epics
        false,      // has_sync
        "root_raw"
    };
}

void RootRawDataSource::fillEventData(fdec::EventData &evt) const
{
    evt.clear();
    evt.info.event_number = ev_.event_num;
    evt.info.trigger_type = ev_.trigger_type;
    evt.info.trigger_bits = ev_.trigger_bits;
    evt.info.timestamp = static_cast<uint64_t>(ev_.timestamp);

    // Without the HyCalSystem we cannot reverse module_id → (crate, slot, ch).
    // Downstream code indexes by ROC/slot/channel, so skip HyCal channels.
    if (!hycal_) return;

    for (int i = 0; i < ev_.nch && i < prad2::kMaxChannels; ++i) {
        const fdec::Module *mod = hycal_->module_by_id(ev_.module_id[i]);
        if (!mod) continue;

        int crate = mod->daq.crate;
        int sl    = mod->daq.slot;
        int ch    = mod->daq.channel;
        if (sl < 0 || sl >= fdec::MAX_SLOTS || ch < 0 || ch >= fdec::MAX_CHANNELS) continue;

        uint32_t roc_tag = static_cast<uint32_t>(crate);
        auto it = crate_to_roc_.find(crate);
        if (it != crate_to_roc_.end()) roc_tag = it->second;

        int roc_idx = -1;
        for (int r = 0; r < evt.nrocs; ++r) {
            if (evt.rocs[r].tag == roc_tag) { roc_idx = r; break; }
        }
        if (roc_idx < 0) {
            if (evt.nrocs >= fdec::MAX_ROCS) continue;
            roc_idx = evt.nrocs++;
            evt.rocs[roc_idx].present = true;
            evt.rocs[roc_idx].tag = roc_tag;
        }

        auto &slot = evt.rocs[roc_idx].slots[sl];
        slot.present = true;
        slot.channel_mask |= (1ull << ch);
        auto &cd = slot.channels[ch];
        cd.nsamples = std::min((int)ev_.nsamples[i], fdec::MAX_SAMPLES);
        for (int s = 0; s < cd.nsamples; ++s)
            cd.samples[s] = ev_.samples[i][s];
    }
}

std::string RootRawDataSource::decodeEvent(int index, fdec::EventData &evt,
                                            ssp::SspEventData *ssp)
{
    if (index < 0 || index >= n_entries_) return "event out of range";
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    fillEventData(evt);
    if (ssp) ssp->clear();
    return "";
}

void RootRawDataSource::iterateAll(EventCallback ev_cb, ReconCallback /*recon_cb*/,
                                    ControlCallback /*ctrl_cb*/, EpicsCallback /*epics_cb*/,
                                    DscCallback /*dsc_cb*/, int /*dsc_bank_tag*/)
{
    if (!tree_ || !ev_cb) return;

    std::lock_guard<std::mutex> lk(mtx_);  // block concurrent decodeEvent calls
    auto event_ptr = std::make_unique<fdec::EventData>();
    for (int i = 0; i < n_entries_; ++i) {
        tree_->GetEntry(i);
        fillEventData(*event_ptr);
        ev_cb(i, *event_ptr, nullptr);
    }
}

std::string RootReconDataSource::open(const std::string &path)
{
    close();
    file_.reset(TFile::Open(path.c_str(), "READ"));
    if (!file_ || file_->IsZombie()) {
        file_.reset();
        return "cannot open ROOT file";
    }
    tree_ = file_->Get<TTree>("recon");
    if (!tree_) { close(); return "no 'recon' tree in ROOT file"; }

    n_entries_ = static_cast<int>(tree_->GetEntries());

    prad2::SetReconReadBranches(tree_, ev_);

    // Recon trees don't carry the run number per event; take it from the
    // file name ("prad_NNNNNN.*_recon.root").
    const int run = prad2::run_number_from_path(path);
    run_number_ = run > 0 ? static_cast<uint32_t>(run) : 0;

    std::cerr << "ROOT recon: " << n_entries_ << " events";
    if (run_number_) std::cerr << "  run=" << run_number_;
    std::cerr << "\n";
    return "";
}

void RootReconDataSource::close()
{
    tree_ = nullptr;
    file_.reset();
    n_entries_ = 0;
}

DataSourceCaps RootReconDataSource::capabilities() const
{
    return {
        false,        // has_waveforms
        false,        // has_peaks
        false,        // has_pedestals
        true,         // has_clusters (pre-computed)
        false,        // has_gem_raw
        true,         // has_gem_hits (pre-computed)
        false,        // has_epics
        false,        // has_sync
        "root_recon"
    };
}

std::string RootReconDataSource::decodeEvent(int index, fdec::EventData &evt,
                                              ssp::SspEventData *ssp)
{
    if (index < 0 || index >= n_entries_) return "event out of range";
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    evt.clear();
    evt.info.event_number = ev_.event_num;
    evt.info.trigger_type = ev_.trigger_type;
    evt.info.trigger_bits = ev_.trigger_bits;
    evt.info.run_number = run_number_;
    evt.info.timestamp = static_cast<uint64_t>(ev_.timestamp);
    if (ssp) ssp->clear();
    return "";
}

bool RootReconDataSource::decodeReconEvent(int index, prad2::ReconEventData &recon)
{
    if (index < 0 || index >= n_entries_) return false;
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    recon = ev_;
    return true;
}

void RootReconDataSource::iterateAll(EventCallback /*ev_cb*/, ReconCallback recon_cb,
                                      ControlCallback /*ctrl_cb*/, EpicsCallback /*epics_cb*/,
                                      DscCallback /*dsc_cb*/, int /*dsc_bank_tag*/)
{
    if (!tree_ || !recon_cb) return;

    std::lock_guard<std::mutex> lk(mtx_);  // block concurrent decodeReconEvent calls
    for (int i = 0; i < n_entries_; ++i) {
        tree_->GetEntry(i);
        recon_cb(i, ev_);
    }
}

#endif // WITH_ROOT
