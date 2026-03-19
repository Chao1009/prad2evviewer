#pragma once
//=============================================================================
// DaqConfig.h — configurable DAQ bank tags and event type identification
//
// All tags are configurable to accommodate DAQ format changes.
// Defaults match PRad-II data (prad_023109.evio format).
//
// This is a plain struct — no JSON dependency. Loading from JSON is handled
// by the application layer (see load_daq_config.h).
//=============================================================================

#include <cstdint>
#include <vector>
#include <string>

namespace evc
{

struct DaqConfig
{
    // --- event type identification (top-level bank tag ranges) ---------------

    // physics event tag range
    // PRad uses 0xFE (single-event CODA mode, num = session ID)
    // CODA built-trigger uses 0xFF50-0xFF8F (num = event count)
    uint32_t physics_tag_min  = 0x00FE;
    uint32_t physics_tag_max  = 0x00FE;

    // control event tags
    uint32_t prestart_tag     = 0x11;
    uint32_t go_tag           = 0x12;
    uint32_t end_tag          = 0x14;

    // sync event tag
    uint32_t sync_tag         = 0xC1;

    // EPICS slow control event tag
    uint32_t epics_tag        = 0x1F;

    // --- bank tags within physics events ------------------------------------

    // FADC250 composite data bank tag
    uint32_t fadc_composite_tag = 0xE101;

    // Trigger Interface (TI) data bank tag (present in every ROC)
    uint32_t ti_bank_tag        = 0xE10A;

    // Trigger/event number bank (top-level child)
    uint32_t trigger_bank_tag   = 0xC000;

    // Run info bank (in TI master crate only)
    uint32_t run_info_tag       = 0xE10F;

    // DAQ configuration readback string bank
    uint32_t daq_config_tag     = 0xE10E;

    // EPICS data bank tag (within EPICS events)
    uint32_t epics_bank_tag     = 0xE114;

    // --- TI data format -----------------------------------------------------
    // TI bank layout: word[0]=header, word[1]=trigger#, word[2]=ts_low, word[3]=ts_high
    int ti_trigger_word   = 1;
    int ti_time_low_word  = 2;      // lower 32 bits of 48-bit timestamp
    int ti_time_high_word = 3;      // upper bits of timestamp (shifted)
    uint32_t ti_time_high_mask  = 0xFFFF0000;
    int      ti_time_high_shift = 16;   // right-shift before combining

    // --- trigger bank format (tag 0xC000) -----------------------------------
    int trig_event_number_word = 0;
    int trig_event_type_word   = 1;

    // --- run info bank format (tag 0xE10F, in TI master crate) --------------
    int ri_run_number_word     = 1;
    int ri_event_count_word    = 2;
    int ri_unix_time_word      = 3;

    // --- ROC identification -------------------------------------------------
    struct RocEntry {
        uint32_t    tag;
        std::string name;
    };
    std::vector<RocEntry> roc_tags;

    // TI master crate tag (contains run info bank)
    uint32_t ti_master_tag = 0x27;

    // --- helpers ------------------------------------------------------------
    bool is_physics(uint32_t tag) const
    {
        return tag >= physics_tag_min && tag <= physics_tag_max;
    }

    bool is_control(uint32_t tag) const
    {
        return tag == prestart_tag || tag == go_tag || tag == end_tag;
    }

    bool is_sync(uint32_t tag) const { return tag == sync_tag; }
    bool is_epics(uint32_t tag) const { return tag == epics_tag; }
};

// --- event type enum --------------------------------------------------------
enum class EventType : uint8_t {
    Unknown   = 0,
    Physics   = 1,
    Sync      = 2,
    Epics     = 3,
    Prestart  = 4,
    Go        = 5,
    End       = 6,
    Control   = 7,
};

inline EventType classify_event(uint32_t tag, const DaqConfig &cfg)
{
    if (cfg.is_physics(tag)) return EventType::Physics;
    if (cfg.is_sync(tag))    return EventType::Sync;
    if (cfg.is_epics(tag))   return EventType::Epics;
    if (tag == cfg.prestart_tag) return EventType::Prestart;
    if (tag == cfg.go_tag)       return EventType::Go;
    if (tag == cfg.end_tag)      return EventType::End;
    if (cfg.is_control(tag))     return EventType::Control;
    return EventType::Unknown;
}

} // namespace evc
