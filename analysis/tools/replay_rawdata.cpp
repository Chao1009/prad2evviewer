//=============================================================================
// replay_rawdata — convert multiple EVIO files to ROOT trees (multi-threaded)
//
// Usage: replay_rawdata <evio_file_or_dir> [more files/dirs...]
//                       -o output_dir [-f max_files] [-n max_events] [-p] [-j num_threads]
//                       [-c daq_config.json] [-d hycal_map.json] [-x] [-z threshold]
//                       [--Ecalib] [--noWaveform]
//   -o  output directory (REQUIRED)
//   -f  max files to process (default: all)
//   -n  max events per file (default: all)
//   -p  include peak analysis branches
//   -j  number of threads (default: 4)
//   -c  DAQ configuration file
//   -d  HyCal map file (default: <db>/hycal_map.json)
//   -x  run the X17 reconstruction path
//   -z  override GEM zero-suppression threshold (sigma)
//   --Ecalib  enable Ecalib mode, throw out extra branches to reduce output size (implies -p)
//   --noWaveform  disable saving raw waveform samples (implies -p)
//=============================================================================

#include "Replay.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "GainCorrCompute.h"
#include "ToolUtils.h"

#include <iostream>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <algorithm>
#include <vector>

using namespace analysis;

static std::string makeOutputFile(const std::string &evio_path)
{
    std::string out = std::filesystem::path(evio_path).filename().string();
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_raw.root";
    return out;
}

static const struct option long_options[] = {
    {"Ecalib", no_argument, nullptr, 1000},
    {"noWaveform", no_argument, nullptr, 1001},
    {nullptr, 0, nullptr, 0}
};

int main(int argc, char *argv[])
{
    analysis::InitRootThreading();

    std::string daq_config, daq_map, output_dir;
    int max_events = -1;
    int max_files = -1;
    bool peaks = false;
    bool Ecalib = false;
    bool noWaveform = false;
    int num_threads = 4;
    float zerosup_override = 5.f;
    bool x17 = false;

    std::string db_dir = prad2::database_dir();
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    opterr = 0;
    int opt;
    while ((opt = getopt_long(argc, argv, "o:f:n:c:d:j:z:px", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'f':
                if (!ParseIntOption(optarg, max_files))
                    return InvalidOptionValue("-f", optarg);
                break;
            case 'n':
                if (!ParseIntOption(optarg, max_events))
                    return InvalidOptionValue("-n", optarg);
                break;
            case 'c': daq_config = optarg; break;
            case 'd': daq_map = optarg; break;
            case 'j':
                if (!ParseIntOption(optarg, num_threads) || num_threads <= 0)
                    return InvalidOptionValue("-j", optarg);
                break;
            case 'p': peaks = true; break;
            case 'x': x17 = true; break;
            case 'z':
                if (!ParseFloatOption(optarg, zerosup_override)
                        || zerosup_override < 0.f)
                    return InvalidOptionValue("-z", optarg);
                break;
            case 1000: Ecalib = true; peaks = true; break;
            case 1001: noWaveform = true; peaks = true; break;
            case '?':
            default:
                return InvalidOption(argv, optind, optopt);
        }
    }

    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> evio_files = CollectInputs(argc, argv, optind, IsEvioName);

    if (evio_files.empty() || output_dir.empty()) {
        std::cerr << "Usage: replay_rawdata <evio_file_or_dir> [more files/dirs...] -o output_dir\n"
                  << "       [-f max_files] [-j threads] [-c daq_config.json] [-d hycal_map.json] "
                  <<         "[-n N] [-z zerosup_override] [-p] [--Ecalib] [--noWaveform]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -f  max files to process (default: all)\n";
        std::cerr << "  -j  number of threads (default: 4)\n";
        std::cerr << "  -c  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -d  HyCal map JSON (default: <db>/hycal_map.json)\n";
        std::cerr << "  -n  max events per file (default: all)\n";
        std::cerr << "  -p  include peak analysis branches (soft + firmware DAQ-mode)\n";
        std::cerr << "  -x  run the X17 reconstruction path\n";
        std::cerr << "  -z  override GEM zero-suppression threshold (sigma)\n";
        std::cerr << "  --Ecalib  enable Ecalib mode, throw out extra branches to reduce output size (including peak analysis)\n";
        std::cerr << "  --noWaveform  disable saving raw waveform samples (including peak analysis)\n";
        return 1;
    }
    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::cout << "Processing " << num_files << " files with "
              << num_threads << " threads\n";

    if(daq_map.empty()) daq_map = db_dir + "/hycal_map.json";

    const std::vector<std::string> inputs(evio_files.begin(), evio_files.begin() + num_files);
    EnsureGainCorr(inputs, db_dir, daq_config, daq_map, num_threads);

    int run_num = get_run_int(evio_files[0]);
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    const std::string recon_config = db_dir + (x17 ? "/reconstruction_config_x17.json"
                                                   : "/reconstruction_config.json");

    std::cerr << "Using HyCal map: " << daq_map << "\n";
    const int errors = RunReplayPool(inputs, num_threads, daq_config, daq_map,
        [&](const std::string &in) { return output_dir + "/" + makeOutputFile(in); },
        [&](Replay &replay, const std::string &in, const std::string &out) {
            return replay.Process(in, out, gRunConfig, db_dir, recon_config, max_events, peaks, daq_config,
                zerosup_override, Ecalib, noWaveform);
        });

    std::cout << "Done: " << num_files << " files"
              << (errors > 0 ? ", " + std::to_string(errors) + " errors" : "")
              << "\n";

    return errors > 0 ? 1 : 0;
}
