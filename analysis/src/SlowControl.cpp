//=============================================================================
// SlowControl.cpp — slow-event row loaders, delta livetime, live charge
//=============================================================================

#include "SlowControl.h"
#include "DscData.h"
#include "EventData_io.h"
#include "Fadc250Data.h"

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <memory>

namespace analysis {

namespace {

// Null, after printing "<tag>: cannot open <path>", when path cannot be read.
std::unique_ptr<TFile> openInput(const std::string &path, const char *tag)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (f && !f->IsZombie()) return f;
    std::cerr << tag << ": cannot open " << path << "\n";
    return nullptr;
}

// Binds the per-row `good` verdict that replay_filter writes, if present.
void bindGood(TTree *t, bool &good, bool *has_good)
{
    if (!t->GetBranch("good")) return;
    t->SetBranchAddress("good", &good);
    if (has_good) *has_good = true;
}

} // anonymous namespace

bool LoadScalerRows(const std::vector<std::string> &files, std::vector<ScalerRow> &out,
                    const char *tag, bool *has_good)
{
    ScalerRow row;
    for (const auto &path : files) {
        const auto f = openInput(path, tag);
        if (!f) return false;
        auto *t = dynamic_cast<TTree *>(f->Get("scalers"));
        if (!t) continue;
        prad2::SetScalerReadBranches(t, row);
        bindGood(t, row.good, has_good);
        const Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            row.good = true;
            t->GetEntry(i);
            out.push_back(row);
        }
    }
    return true;
}

bool LoadEpicsRows(const std::vector<std::string> &files, std::vector<EpicsRow> &out,
                   const char *tag, bool *has_good)
{
    for (const auto &path : files) {
        const auto f = openInput(path, tag);
        if (!f) return false;
        auto *t = dynamic_cast<TTree *>(f->Get("epics"));
        if (!t) continue;
        // Fresh per file: branches missing from older files read as 0.
        prad2::RawEpicsData ep;
        prad2::EpicsVectorBindings vectors;
        prad2::SetEpicsReadBranches(t, ep);
        prad2::BindEpicsVectorBranches(t, ep, vectors);
        bool good = true;
        bindGood(t, good, has_good);
        const Long64_t n = t->GetEntries();
        out.reserve(out.size() + n);
        for (Long64_t i = 0; i < n; ++i) {
            good = true;
            t->GetEntry(i);
            EpicsRow r;
            r.event_number = ep.event_number_at_arrival;
            r.ti_ticks     = ep.ti_ticks_at_arrival;
            r.unix_time    = ep.unix_time;
            r.sync_counter = ep.sync_counter;
            r.run_number   = ep.run_number;
            r.good         = good;
            const size_t k_max = std::min(ep.channel.size(), ep.value.size());
            for (size_t k = 0; k < k_max; ++k)
                r.updates[ep.channel[k]] = ep.value[k];
            out.push_back(std::move(r));
        }
    }
    return true;
}

std::pair<uint32_t, uint32_t> SelectDscPair(const prad2::RawScalerData &row,
                                            const std::string &source, int channel)
{
    if (source == "ref") return {row.ref_gated, row.ref_ungated};
    const int c = std::clamp(channel, 0, prad2::kDscChannels - 1);
    if (source == "trg") return {row.trg_gated[c], row.trg_ungated[c]};
    if (source == "tdc") return {row.tdc_gated[c], row.tdc_ungated[c]};
    return {0, 0};
}

std::vector<double> DeltaLivetime(const std::vector<ScalerRow> &rows,
                                  const std::vector<size_t> &order,
                                  const std::string &source, int channel,
                                  double scale)
{
    std::vector<double> out(rows.size(), -1.0);
    uint32_t prev_g = 0, prev_u = 0;
    for (const size_t i : order) {
        const auto [g, u] = SelectDscPair(rows[i], source, channel);
        const double lf = dsc::delta_live_ratio(g, u, prev_g, prev_u);
        if (lf >= 0) out[i] = lf * scale;
        prev_g = g;
        prev_u = u;
    }
    return out;
}

void ChargeSums::AddPair(int64_t ticks_a, int64_t ticks_b, double live_fraction_b,
                         double current_a, double current_b, bool good)
{
    const bool data_ok = !(ticks_a <= 0 || ticks_b <= 0 || ticks_b <= ticks_a
        || !std::isfinite(live_fraction_b) || live_fraction_b < 0
        || !std::isfinite(current_a) || !std::isfinite(current_b));
    if (!data_ok) {
        if (good) ++n_pairs_skipped;
        ++n_ungated_pairs_skipped;
        return;
    }
    const double dt = (ticks_b - ticks_a) * fdec::TI_TICK_SEC;
    const double I  = 0.5 * (current_a + current_b);
    const double dQ = live_fraction_b * dt * I;
    const double dL = live_fraction_b * dt;
    ungated_value_nC     += dQ;
    ungated_live_seconds += dL;
    ungated_real_seconds += dt;
    ++n_ungated_pairs_integrated;
    if (good) {
        value_nC     += dQ;
        live_seconds += dL;
        real_seconds += dt;
        ++n_pairs_integrated;
    }
}

ChargeSums &ChargeSums::operator+=(const ChargeSums &o)
{
    value_nC                   += o.value_nC;
    live_seconds               += o.live_seconds;
    real_seconds               += o.real_seconds;
    n_pairs_integrated         += o.n_pairs_integrated;
    n_pairs_skipped            += o.n_pairs_skipped;
    ungated_value_nC           += o.ungated_value_nC;
    ungated_live_seconds       += o.ungated_live_seconds;
    ungated_real_seconds       += o.ungated_real_seconds;
    n_ungated_pairs_integrated += o.n_ungated_pairs_integrated;
    n_ungated_pairs_skipped    += o.n_ungated_pairs_skipped;
    return *this;
}

} // namespace analysis
