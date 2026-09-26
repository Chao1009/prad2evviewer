#pragma once
//=============================================================================
// SspData.h — pre-allocated, flat event data for SSP/MPD/APV GEM readout
//
// No heap allocation in the event loop. Fixed-size arrays sized for worst case.
// APV indexed by (crate, mpd/fiber, adc_ch) for O(1) access.
//=============================================================================

#include <cstdint>
#include <cstring>

namespace ssp
{

// --- capacity limits --------------------------------------------------------
static constexpr int APV_STRIP_SIZE     = 128;   // channels per APV25 chip
static constexpr int SSP_TIME_SAMPLES   = 6;     // fixed by SSP firmware
static constexpr int MAX_APVS_PER_MPD   = 16;    // APV slots per MPD
static constexpr int MAX_MPDS           = 64;    // MPDs across all crates

// --- APV address ------------------------------------------------------------
struct ApvAddress {
    int crate_id = -1;
    int mpd_id   = -1;   // fiber ID in SSP readout
    int adc_ch   = -1;   // APV ID within MPD

    bool operator==(const ApvAddress &o) const
    {
        return crate_id == o.crate_id && mpd_id == o.mpd_id && adc_ch == o.adc_ch;
    }
};

// --- per-APV data -----------------------------------------------------------
struct ApvData {
    ApvAddress addr;
    bool     present = false;

    // Raw ADC samples: strips[strip][time_sample]
    // 13-bit sign-extended values from SSP firmware
    int16_t  strips[APV_STRIP_SIZE][SSP_TIME_SAMPLES];

    int      nstrips = 0;           // count of populated strips
    uint64_t strip_mask[2] = {};    // bitmask: bit i set = strip i has data
                                    // strip_mask[0] = strips 0-63, [1] = 64-127

    // SSP firmware flags (from MPD frame header)
    uint32_t flags = 0;

    // Online common mode from firmware debug header (6 values, one per time sample)
    int16_t  online_cm[SSP_TIME_SAMPLES] = {};
    bool     has_online_cm = false;

    void clear()
    {
        present = false;
        nstrips = 0;
        strip_mask[0] = strip_mask[1] = 0;
        flags = 0;
        has_online_cm = false;
        // Zero strips fully, not lazily via strip_mask: some read paths
        // (e.g. the Python `apv.strips` binding) copy the full 128×6
        // buffer without gating on hasStrip(), so stale samples would
        // leak between events.  Costs one ~1.5 KB memset per APV.
        std::memset(strips, 0, sizeof(strips));
    }

    void setStrip(int strip, int ts, int16_t value)
    {
        strips[strip][ts] = value;
        int idx = strip >> 6;   // 0 or 1
        int bit = strip & 63;
        if (!(strip_mask[idx] & (1ULL << bit))) {
            strip_mask[idx] |= (1ULL << bit);
            ++nstrips;
        }
    }

    bool hasStrip(int strip) const
    {
        int idx = strip >> 6;
        int bit = strip & 63;
        return (strip_mask[idx] & (1ULL << bit)) != 0;
    }

    // Full readout: the firmware sent every strip (no online zero
    // suppression).  nstrips is the authoritative signal; has_online_cm is
    // not — the MPD can emit CM debug headers while still sending all strips.
    bool isFullReadout() const { return present && nstrips >= APV_STRIP_SIZE; }
};

// --- per-MPD data -----------------------------------------------------------
struct MpdData {
    int      crate_id = -1;
    int      mpd_id   = -1;
    bool     present  = false;
    int      napvs    = 0;
    ApvData  apvs[MAX_APVS_PER_MPD];   // indexed by APV ID (adc_ch)

    void clear()
    {
        present = false;
        napvs = 0;
        crate_id = -1;
        mpd_id = -1;
        for (int i = 0; i < MAX_APVS_PER_MPD; ++i)
            apvs[i].clear();
    }
};

// --- full SSP event data ----------------------------------------------------
struct SspEventData {
    int      nmpds = 0;
    MpdData  mpds[MAX_MPDS];

    void clear()
    {
        for (int i = 0; i < nmpds; ++i)
            mpds[i].clear();
        nmpds = 0;
    }

    // Find or create MPD entry by (crate, mpd_id).
    // Returns pointer to MpdData, or nullptr if full.
    MpdData *findOrCreateMpd(int crate, int mpd)
    {
        for (int i = 0; i < nmpds; ++i) {
            if (mpds[i].crate_id == crate && mpds[i].mpd_id == mpd)
                return &mpds[i];
        }
        if (nmpds >= MAX_MPDS) return nullptr;
        MpdData &m = mpds[nmpds++];
        m.present = true;
        m.crate_id = crate;
        m.mpd_id = mpd;
        return &m;
    }

    // True if any APV in the event is a full readout (ApvData::isFullReadout).
    bool hasFullReadout() const
    {
        for (int i = 0; i < nmpds; ++i)
            for (const ApvData &apv : mpds[i].apvs)
                if (apv.isFullReadout()) return true;
        return false;
    }

    // Call fn(const MpdData &mpd, int apv_id, const ApvData &apv) for every
    // present APV of every present MPD, in readout order; apv_id is the
    // index into mpd.apvs (= adc_ch).
    template <class Fn>
    void forEachApv(Fn &&fn) const
    {
        for (int i = 0; i < nmpds; ++i) {
            const MpdData &mpd = mpds[i];
            if (!mpd.present) continue;
            for (int a = 0; a < MAX_APVS_PER_MPD; ++a)
                if (mpd.apvs[a].present) fn(mpd, a, mpd.apvs[a]);
        }
    }

    // Find APV by full address. Returns nullptr if not found.
    const ApvData *findApv(int crate, int mpd, int adc) const
    {
        for (int i = 0; i < nmpds; ++i) {
            if (mpds[i].crate_id == crate && mpds[i].mpd_id == mpd) {
                if (adc >= 0 && adc < MAX_APVS_PER_MPD && mpds[i].apvs[adc].present)
                    return &mpds[i].apvs[adc];
                return nullptr;
            }
        }
        return nullptr;
    }
};

} // namespace ssp
