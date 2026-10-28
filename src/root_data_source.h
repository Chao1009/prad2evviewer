#pragma once
// ROOT file data sources (requires WITH_ROOT):
//   RootRawDataSource:   reads "events" tree from replay_rawdata output
//   RootReconDataSource: reads "recon" tree from replay_recon output
//
// Uses shared data structs from EventData.h — any schema change there
// automatically updates both the writer (Replay) and these readers.
//
// The raw source needs a HyCalSystem pointer to translate the stored
// module_id back to DAQ (crate, slot, channel) addressing so downstream
// code (which indexes by ROC/slot/channel) still works.

#ifdef WITH_ROOT

#include "data_source.h"
#include "EventData.h"
#include "HyCalSystem.h"

#include <TFile.h>
#include <TTree.h>

#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

class RootRawDataSource : public DataSource {
public:
    RootRawDataSource(const std::unordered_map<int, uint32_t> &crate_to_roc,
                      const fdec::HyCalSystem *hycal)
        : crate_to_roc_(crate_to_roc), hycal_(hycal) {}

    std::string open(const std::string &path) override;
    void close() override;
    DataSourceCaps capabilities() const override;
    int eventCount() const override { return n_entries_; }
    std::string decodeEvent(int index, fdec::EventData &evt,
                             ssp::SspEventData *ssp = nullptr) override;
    void iterateAll(EventCallback ev_cb, ReconCallback recon_cb,
                    ControlCallback ctrl_cb, EpicsCallback epics_cb,
                    DscCallback dsc_cb, int dsc_bank_tag) override;

private:
    std::unordered_map<int, uint32_t> crate_to_roc_;
    const fdec::HyCalSystem *hycal_ = nullptr;
    std::unique_ptr<TFile> file_;
    TTree *tree_ = nullptr;
    int n_entries_ = 0;
    bool has_peaks_ = false;
    bool has_gem_ = false;
    std::mutex mtx_;

    prad2::RawEventData ev_;     // branch read buffer

    void fillEventData(fdec::EventData &evt) const;
};

class RootReconDataSource : public DataSource {
public:
    std::string open(const std::string &path) override;
    void close() override;
    DataSourceCaps capabilities() const override;
    int eventCount() const override { return n_entries_; }
    std::string decodeEvent(int index, fdec::EventData &evt,
                             ssp::SspEventData *ssp = nullptr) override;
    bool decodeReconEvent(int index, prad2::ReconEventData &recon) override;
    uint32_t runNumber() const override { return run_number_; }
    void iterateAll(EventCallback ev_cb, ReconCallback recon_cb,
                    ControlCallback ctrl_cb, EpicsCallback epics_cb,
                    DscCallback dsc_cb, int dsc_bank_tag) override;

private:
    std::unique_ptr<TFile> file_;
    TTree *tree_ = nullptr;
    int n_entries_ = 0;
    uint32_t run_number_ = 0;   // parsed from filename "prad_NNNNNN.*_recon.root"
    std::mutex mtx_;

    prad2::ReconEventData ev_;   // branch read buffer
};

// Factory helper: detect tree type and return appropriate source.
// hycal is required for "events" tree decoding (module_id → DAQ reverse lookup).
std::unique_ptr<DataSource> createRootDataSource(
    const std::string &path,
    const std::unordered_map<int, uint32_t> &crate_to_roc,
    const fdec::HyCalSystem *hycal);

#endif // WITH_ROOT
