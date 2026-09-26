#include "evio_data_source.h"

#include <memory>

using namespace evc;

std::string EvioDataSource::open(const std::string &path)
{
    close();
    filepath_ = path;

    if (!reopenReader()) {
        filepath_.clear();
        return "cannot open file";
    }

    // Index every non-monitoring event so the EPICS tab + control-event
    // bookkeeping see them.  The viewer's HTTP layer tags non-Physics samples
    // in the JSON response (`event_kind`) so the status bar can label them
    // instead of showing "0 channels, no trigger".
    if (reader_.IsRandomAccess()) {
        int n_evio = reader_.GetRandomAccessEventCount();
        for (int ei = 0; ei < n_evio; ++ei)
            if (reader_.ReadEventByIndex(ei) == status::success) indexRecord(ei);
    } else {
        for (int ei = 0; reader_.Read() == status::success; ++ei) {
            reader_pos_ = ei;
            indexRecord(ei);
        }
    }
    return "";
}

void EvioDataSource::close()
{
    index_.clear();
    invalidateReader();
    filepath_.clear();
}

void EvioDataSource::indexRecord(int ei)
{
    if (!reader_.Scan() || cfg_.is_monitoring(reader_.GetEvHeader().tag)) return;
    const EventType et = reader_.GetEventType();
    for (int si = 0; si < reader_.GetNEvents(); ++si)
        index_.push_back({ei, si, et});
}

DataSourceCaps EvioDataSource::nativeCaps()
{
    return {
        true,   // has_waveforms
        true,   // has_peaks (computed by WaveAnalyzer)
        true,   // has_pedestals
        true,   // has_clusters (computed by HyCalCluster)
        true,   // has_gem_raw
        true,   // has_gem_hits (computed by GemCluster)
        true,   // has_epics
        true,   // has_sync
        "evio"  // source_type
    };
}

void EvioDataSource::invalidateReader()
{
    reader_.Close();
    reader_path_.clear();
    reader_pos_ = -1;
    last_decoded_index_ = -1;
}

bool EvioDataSource::reopenReader()
{
    // OpenAuto closes a still-open handle itself.
    if (reader_.OpenAuto(filepath_) != status::success) {
        invalidateReader();
        return false;
    }
    reader_path_ = filepath_;
    reader_pos_ = -1;
    return true;
}

std::string EvioDataSource::seekTo(int evio_target)
{
    // Target is behind us — reopen and walk forward from the start.
    if (reader_pos_ > evio_target && !reopenReader())
        return "cannot reopen file";
    while (reader_pos_ < evio_target) {
        if (reader_.Read() != status::success) {
            invalidateReader();
            return "read error while seeking";
        }
        ++reader_pos_;
    }
    return "";
}

std::string EvioDataSource::decodeEvent(int index, fdec::EventData &evt,
                                         ssp::SspEventData *ssp)
{
    if (index < 0 || index >= (int)index_.size())
        return "event out of range";

    std::lock_guard<std::mutex> lk(reader_mtx_);

    // Cache-hit path: same event as last decode → reader_'s lazy cache is
    // already valid, skip Scan + decode and just copy out.
    if (index != last_decoded_index_) {
        // Recover the reader if something invalidated it (e.g. prior error).
        if (reader_path_ != filepath_ && !reopenReader())
            return "cannot open file";
        const auto &ei = index_[index];
        if (reader_.IsRandomAccess()) {
            if (reader_.ReadEventByIndex(ei.evio_event) != status::success) {
                last_decoded_index_ = -1;
                return "read error";
            }
        } else {
            std::string err = seekTo(ei.evio_event);
            if (!err.empty()) { last_decoded_index_ = -1; return err; }
        }
        if (!reader_.Scan()) { last_decoded_index_ = -1; return "scan error"; }
        reader_.SelectEvent(ei.sub_event);
        last_decoded_index_ = index;
    }

    evt = reader_.Fadc();                      // lazy-decodes on first access,
    if (ssp) *ssp = reader_.Gem();             // returns cached ref thereafter
    return "";
}

void EvioDataSource::iterateAll(EventCallback ev_cb, ReconCallback /*recon_cb*/,
                                ControlCallback ctrl_cb, EpicsCallback epics_cb,
                                DscCallback dsc_cb, int dsc_bank_tag)
{
    EvChannel ch;
    ch.SetConfig(cfg_);
    if (ch.OpenAuto(filepath_) != status::success) return;

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    auto ssp_ptr = std::make_unique<ssp::SspEventData>();
    auto &ssp_evt = *ssp_ptr;
    uint64_t last_ti_ts = 0;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;

        // control events (sync/prestart/go/end) — absolute unix time lands in
        // ch.Sync() along with run number / type / counter; Sync()'s snapshot
        // persists across events, so gate on event type so physics events
        // don't re-fire the callback.
        if (ctrl_cb) {
            auto et = ch.GetEventType();
            if (et == EventType::Prestart || et == EventType::Go ||
                et == EventType::End      || et == EventType::Sync)
            {
                const auto &s = ch.Sync();
                if (s.unix_time != 0) ctrl_cb(s.unix_time, last_ti_ts);
            }
        }

        if (epics_cb && ch.GetEventType() == EventType::Epics) {
            std::string text = ch.ExtractEpicsText();
            if (!text.empty())
                epics_cb(text, 0, last_ti_ts);
        }

        // DSC2 scaler bank.  Bank-format detection + (source, channel)
        // selection live in EvChannel::Dsc(); fire only when a configured
        // slot matched.
        if (dsc_cb && dsc_bank_tag >= 0) {
            auto et = ch.GetEventType();
            if (et == EventType::Sync || et == EventType::Physics) {
                const auto &dsc = ch.Dsc();
                if (dsc.present) dsc_cb(dsc);
            }
        }

        // physics events (skip monitoring — no waveforms to process)
        if (cfg_.is_monitoring(ch.GetEvHeader().tag)) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            if (!ch.DecodeEvent(i, event, &ssp_evt)) continue;
            last_ti_ts = event.info.timestamp;
            if (ev_cb) ev_cb(i, event, &ssp_evt);
        }
    }
    ch.Close();
}
