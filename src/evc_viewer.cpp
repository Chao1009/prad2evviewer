// src/evc_viewer.cpp — HyCal event viewer
//
// Single C++ binary: evio decoder + waveform analysis + HTTP server.
// Auto-discovers database files (hycal_modules.json, daq_map.json) and
// resources (viewer.html) from compile-time DATABASE_DIR / RESOURCE_DIR.
//
// Usage: evc_viewer <evio_file> [port]

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "Fadc250Decoder.h"
#include "WaveAnalyzer.h"

#include <nlohmann/json.hpp>

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

using json = nlohmann::json;
using WsServer = websocketpp::server<websocketpp::config::asio>;
using namespace evc;

// compile-time defaults (set by CMake)
#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif
#ifndef RESOURCE_DIR
#define RESOURCE_DIR "."
#endif

// -------------------------------------------------------------------------
// Globals
// -------------------------------------------------------------------------
struct EventIndex { int buffer_num, sub_event; };

static std::string g_filepath;
static std::vector<EventIndex> g_index;
static std::string g_viewer_html;
static json g_config;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
static std::string readFile(const std::string &path)
{
    std::ifstream f(path);
    if (!f) return "";
    return {std::istreambuf_iterator<char>(f), {}};
}

// try path as-is, then relative to a base directory
static std::string findFile(const std::string &name, const std::string &base_dir)
{
    // as given
    { std::ifstream f(name); if (f.good()) return name; }
    // under base_dir
    std::string p = base_dir + "/" + name;
    { std::ifstream f(p); if (f.good()) return p; }
    return "";
}

// -------------------------------------------------------------------------
// Index the evio file
// -------------------------------------------------------------------------
static void buildIndex(const std::string &path)
{
    g_filepath = path;
    g_index.clear();

    EvChannel ch;
    if (ch.Open(path) != status::success) {
        std::cerr << "Failed to open " << path << "\n";
        return;
    }

    int buf = 0;
    while (ch.Read() == status::success) {
        ++buf;
        if (!ch.Scan()) continue;
        for (int i = 0; i < ch.GetNEvents(); ++i)
            g_index.push_back({buf, i});
    }
    ch.Close();
    std::cerr << "Indexed " << g_index.size() << " events in " << buf << " buffers\n";
}

// -------------------------------------------------------------------------
// Decode one event → JSON
// -------------------------------------------------------------------------
static json decodeEvent(int ev1)
{
    int idx = ev1 - 1;
    if (idx < 0 || idx >= (int)g_index.size())
        return {{"error", "event out of range"}};

    auto &ei = g_index[idx];
    EvChannel ch;
    if (ch.Open(g_filepath) != status::success)
        return {{"error", "cannot open file"}};

    for (int b = 0; b < ei.buffer_num; ++b)
        if (ch.Read() != status::success) { ch.Close(); return {{"error", "read error"}}; }
    if (!ch.Scan()) { ch.Close(); return {{"error", "scan error"}}; }

    fdec::EventData event;
    if (!ch.DecodeEvent(ei.sub_event, event)) { ch.Close(); return {{"error", "decode error"}}; }
    ch.Close();

    fdec::WaveAnalyzer ana;
    fdec::WaveResult wres;
    json channels = json::object();

    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];

            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1u << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);

                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                json sarr = json::array();
                for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);

                json parr = json::array();
                for (int p = 0; p < wres.npeaks; ++p) {
                    auto &pk = wres.peaks[p];
                    parr.push_back({
                        {"p", pk.pos}, {"t", std::round(pk.time * 10) / 10},
                        {"h", std::round(pk.height * 10) / 10},
                        {"i", std::round(pk.integral * 10) / 10},
                        {"l", pk.left}, {"r", pk.right},
                        {"o", pk.overflow ? 1 : 0},
                    });
                }

                channels[key] = {
                    {"s", sarr},
                    {"pm", std::round(wres.ped.mean * 10) / 10},
                    {"pr", std::round(wres.ped.rms * 10) / 10},
                    {"pk", parr},
                };
            }
        }
    }
    return {{"event", ev1}, {"channels", channels}};
}

// -------------------------------------------------------------------------
// HTTP handler
// -------------------------------------------------------------------------
static void onHttp(WsServer *srv, websocketpp::connection_hdl hdl)
{
    auto con = srv->get_con_from_hdl(hdl);
    std::string uri = con->get_resource();

    if (uri == "/") {
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(g_viewer_html);
        con->append_header("Content-Type", "text/html; charset=utf-8");
        return;
    }

    if (uri == "/api/config") {
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(g_config.dump());
        con->append_header("Content-Type", "application/json");
        return;
    }

    if (uri.rfind("/api/event/", 0) == 0) {
        int evnum = std::atoi(uri.c_str() + 11);
        con->set_status(websocketpp::http::status_code::ok);
        con->set_body(decodeEvent(evnum).dump());
        con->append_header("Content-Type", "application/json");
        return;
    }

    con->set_status(websocketpp::http::status_code::not_found);
    con->set_body("404 Not Found");
}

// -------------------------------------------------------------------------
// Main
// -------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <evio_file> [port]\n";
        return 1;
    }

    std::string evio_file = argv[1];
    int port = (argc >= 3 && argv[2][0] != '-') ? std::atoi(argv[2]) : 5050;

    // index events
    buildIndex(evio_file);

    // auto-discover files from compile-time paths
    std::string db_dir  = DATABASE_DIR;
    std::string res_dir = RESOURCE_DIR;

    std::string html_file = findFile("viewer.html", res_dir);
    std::string mod_file  = findFile("hycal_modules.json", db_dir);
    std::string daq_file  = findFile("daq_map.json", db_dir);

    g_viewer_html = readFile(html_file);
    if (g_viewer_html.empty())
        std::cerr << "Warning: viewer.html not found (tried " << res_dir << ")\n";

    // build config
    json modules_j = json::array(), daq_j = json::array();
    { std::string s = readFile(mod_file);  if (!s.empty()) modules_j = json::parse(s, nullptr, false); }
    { std::string s = readFile(daq_file);  if (!s.empty()) daq_j     = json::parse(s, nullptr, false); }

    g_config = {
        {"modules", modules_j},
        {"daq", daq_j},
        {"crate_roc", {{"0",0x80},{"1",0x82},{"2",0x84},{"3",0x86},{"4",0x88},{"5",0x8a},{"6",0x8c}}},
        {"total_events", (int)g_index.size()},
    };

    std::cerr << "Database  : " << db_dir << " ("
              << modules_j.size() << " modules, " << daq_j.size() << " DAQ channels)\n"
              << "Resources : " << res_dir << "\n";

    // start server
    WsServer server;
    server.set_access_channels(websocketpp::log::alevel::none);
    server.set_error_channels(websocketpp::log::elevel::warn | websocketpp::log::elevel::rerror);
    server.init_asio();
    server.set_reuse_addr(true);
    server.set_http_handler([&server](websocketpp::connection_hdl hdl) { onHttp(&server, hdl); });
    server.listen(port);
    server.start_accept();

    std::cout << "Viewer at http://localhost:" << port << "\n"
              << "  " << g_index.size() << " events, Ctrl+C to stop\n";

    server.run();
    return 0;
}
