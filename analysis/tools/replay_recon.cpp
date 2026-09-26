//=============================================================================
// replay_recon — convert multiple EVIO files to reconstructed ROOT trees (multi-threaded)
//
// Usage: prad2ana_replay_recon <evio(or _raw.root)_file_or_dir> [more files/dirs...]
//                     -o output_dir [-f max_files] [-prad1] [-j num_threads]
//                     [-r recon_config.json] [-c daq_config.json] [-d daq_map.json]
//                     [-g gem_pedestal.json] [-z zerosup_threshold] [-m merge_files]
//                     [-x17] [-x17_full] [-random] [-gem_hit]
//   -o  output directory (REQUIRED)
//   -f  max files to process (default: all)
//   -prad1  read PRad-1 data and do not include GEM
//   -j  number of threads (default: 4)
//   -r  reconstruction config JSON
//   -c  DAQ configuration file
//   -d  HyCal map file (default: <db>/hycal_map.json)
//   -g  GEM pedestal file
//   -z  zero-suppression threshold override (not available for reading raw files)
//   -m  merge this many split ROOT files per hadd output (default: 62; 0 disables)
//   -x17  run the X17 reconstruction path (default without a mode option: PRad2)
//   -x17_full  run the X17 full reconstruction path
//   -random enable reconstruction for random trigger runs
//   -gem_hit  include GEM hit-level branches in the output (default: disabled)
//=============================================================================

#include "Replay.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "GainCorrCompute.h"
#include "ToolUtils.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <getopt.h>
#include <filesystem>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <atomic>
#include <mutex>

using namespace analysis;

static bool isRawReplayFile(const std::string &path)
{
    const auto name = std::filesystem::path(path).filename().string();
    return name.find("_raw") != std::string::npos
           && name.find(".root") != std::string::npos;
}

static std::string makeOutputFile(const std::string &input_path)
{
    std::string out = std::filesystem::path(input_path).filename().string();
    if (isRawReplayFile(input_path)) {
        // Preserve foo_raw.root -> foo_recon.root, while also producing a
        // valid name for variants such as foo_raw_filter.root.
        const auto root_pos = out.rfind(".root");
        out.resize(root_pos);
        if (out.size() >= 4
            && out.compare(out.size() - 4, 4, "_raw") == 0)
            out.resize(out.size() - 4);
        return out + "_recon.root";
    }
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_recon.root";
    return out;
}

int main(int argc, char *argv[])
{
    analysis::InitRootThreading();

    std::string recon_config,daq_config, daq_map, gem_ped_file, output_dir;
    float zerosup_override = 5.f;
    int max_files = -1;
    int num_threads = 4;
    int merge_batch_size = 62;
    bool prad1 = false;
    bool x17 = false;
    bool x17_blind = true;
    bool random = false;
    bool gem_hit = false;

    std::string db_dir = prad2::database_dir();
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    static const option long_options[] = {
        {"prad1", no_argument, nullptr, 1001},
        {"x17", no_argument, nullptr, 1000},
        {"x17_full", no_argument, nullptr, 1002},
        {"random", no_argument, nullptr, 1003},
        {"gem_hit", no_argument, nullptr, 1004},
        {nullptr, 0, nullptr, 0}
    };

    opterr = 0;
    int opt;
    while ((opt = getopt_long_only(argc, argv, "o:f:r:c:d:j:g:z:m:",
                                   long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'f':
                if (!ParseIntOption(optarg, max_files))
                    return InvalidOptionValue("-f", optarg);
                break;
            case 'r': recon_config = optarg; break;
            case 'c': daq_config = optarg; break;
            case 'd': daq_map = optarg; break;
            case 'j':
                if (!ParseIntOption(optarg, num_threads) || num_threads <= 0)
                    return InvalidOptionValue("-j", optarg);
                break;
            case 'g': gem_ped_file = optarg; break;
            case 'z':
                if (!ParseFloatOption(optarg, zerosup_override)
                        || zerosup_override < 0.f)
                    return InvalidOptionValue("-z", optarg);
                break;
            case 'm':
                if (!ParseIntOption(optarg, merge_batch_size)
                        || merge_batch_size < 0)
                    return InvalidOptionValue("-m", optarg);
                break;
            case 1000: x17 = true; break;
            case 1001: prad1 = true; break;
            case 1002: x17_blind = false; x17 = true; break; // -x17_full
            case 1003: random = true; break; // -random
            case 1004: gem_hit = true; break; // -gem_hit
            case '?':
            default:
                return InvalidOption(argv, optind, optopt);
        }
    }

    if (prad1 && x17) {
        std::cerr << "Options -prad1 and -x17 cannot be used together\n";
        return 2;
    }
    if (random && (prad1 || x17)) {
        std::cerr << "Option -random cannot be used with -prad1 or an X17 mode\n";
        return 2;
    }

    // Collect EVIO and replay_raw ROOT inputs from files, directories, or a mix.
    std::vector<std::string> input_files = CollectInputs(argc, argv, optind,
        [](const std::string &name) { return IsEvioName(name) || isRawReplayFile(name); });

    if (input_files.empty() || output_dir.empty()) {
        std::cerr << "Usage: prad2ana_replay_recon <evio_or_raw_file_or_dir> [more files/dirs...] -o output_dir\n"
                  << "       [-f max_files] [-j threads] [-c daq_config.json] [-d daq_map.json]\n"
                  << "       [-g gem_ped.json] [-z threshold] [-m merge_files] [-prad1] [-x17] [-x17_full] [-random] [-gem_hit]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -f  max files to process (default: all)\n";
        std::cerr << "  -j  number of threads (default: 4)\n";
        std::cerr << "  -r  reconstruction config JSON\n";
        std::cerr << "  -c  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -d  HyCal map JSON (default: <db>/hycal_map.json)\n";
        std::cerr << "  -g  GEM pedestal JSON\n";
        std::cerr << "  -z  zero-suppression threshold override (not available for reading raw files)\n";
        std::cerr << "  -m  merge this many split ROOT files per hadd output (default: 62; 0 disables)\n";
        std::cerr << "  default  PRad2 mode\n";
        std::cerr << "  -prad1  PRad-1 mode (no GEM)\n";
        std::cerr << "  -x17  run the X17 reconstruction path\n";
        std::cerr << "  -x17_full  run the X17 full reconstruction path\n";
        std::cerr << "  -random  enable reconstruction for random trigger runs\n";
        std::cerr << "  -gem_hit  include GEM hit-level branches in the output (default: disabled)\n";
        return 1;
    }
    if (prad1 && std::any_of(input_files.begin(), input_files.end(), isRawReplayFile)) {
        std::cerr << "Raw ROOT reconstruction does not support -prad1\n";
        return 1;
    }
    if (x17 && x17_blind) {
        std::cout << "X17 blind reconstruction mode enabled\n";
    }
    if (x17 && !x17_blind) {
        std::cout << "X17 full reconstruction mode enabled\n";
    }
    if (random) {
        std::cout << "Random trigger reconstruction mode enabled\n";
    }
    std::cout << "Replay mode: " << (x17 ? "X17" : prad1 ? "PRad1" : "PRad2") << "\n";
    if (gem_hit) {
        std::cout << "GEM hit-level branches will be included in the output\n";
    } else {
        std::cout << "GEM hit-level branches will NOT be included in the output\n";
    }
    if(recon_config.empty()) {
        if(x17) recon_config = db_dir + "/reconstruction_config_x17.json";
        else recon_config = db_dir + "/reconstruction_config.json";
    }
    int num_files = static_cast<int>(input_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::cout << "Processing " << num_files << " files with "
              << num_threads << " threads\n";
    
    if(daq_map.empty()) daq_map = db_dir + "/hycal_map.json";

    const std::vector<std::string> inputs(input_files.begin(), input_files.begin() + num_files);

    // Only direct EVIO input can be passed to replay_gainCorr.  replay_raw
    // input is assumed to have been produced after that step.
    std::vector<std::string> evio_inputs;
    std::copy_if(inputs.begin(), inputs.end(), std::back_inserter(evio_inputs),
                 [](const std::string &f) { return !isRawReplayFile(f); });
    EnsureGainCorr(evio_inputs, db_dir, daq_config, daq_map, num_threads);

    int run_num = get_run_int(input_files[0]);
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    const auto output_for = [&](const std::string &in) { return output_dir + "/" + makeOutputFile(in); };
    std::vector<char> output_ok;
    std::cerr << "Using HyCal map: " << daq_map << "\n";
    std::atomic<int> errors{RunReplayPool(inputs, num_threads, daq_config, daq_map, output_for,
        [&](Replay &replay, const std::string &in, const std::string &out) {
            return isRawReplayFile(in)
                ? replay.ProcessRaw2Recon(in, out, gRunConfig, db_dir, recon_config,
                                          daq_config, gem_ped_file, x17, x17_blind, random, gem_hit)
                : replay.ProcessWithRecon(in, out, gRunConfig, db_dir, recon_config,
                                          daq_config, gem_ped_file, zerosup_override,
                                          prad1, x17, x17_blind, random, gem_hit);
        }, &output_ok)};

    std::cout << "Done: " << num_files << " files"
              << (errors > 0 ? ", " + std::to_string(errors.load()) + " errors" : "")
              << "\n";

    if (merge_batch_size == 0)
        return errors > 0 ? 1 : 0;

    std::map<int, std::vector<std::string>> outputs_by_run;
    for (int i = 0; i < num_files; ++i)
        if (output_ok[i])
            outputs_by_run[get_run_int(input_files[i])].push_back(output_for(input_files[i]));

    struct MergeJob { std::string output; std::vector<std::string> inputs; };
    std::vector<MergeJob> merge_jobs;
    for (const auto &[rn, files] : outputs_by_run) {
        for (std::size_t begin = 0, batch = 0; begin < files.size(); begin += merge_batch_size, ++batch) {
            std::size_t end = std::min<std::size_t>(begin + merge_batch_size, files.size());
            if (end - begin < 2)
                continue;

            char name[64];
            std::snprintf(name, sizeof(name), "/prad_%06d_recon_%03zu.root", rn, batch);
            merge_jobs.push_back({
                output_dir + name,
                std::vector<std::string>(files.begin() + begin, files.begin() + end)
            });
        }
    }

    if (!merge_jobs.empty()) {
        std::size_t num_merge_threads = std::min<std::size_t>(merge_jobs.size(), num_threads);
        std::cout << "Merging " << merge_jobs.size()
                  << " batch(es) with " << num_merge_threads
                  << " hadd process(es)\n";

        std::mutex io_mtx;
        ParallelFor(merge_jobs.size(), num_threads, [&](std::size_t idx, int) {
            const auto &job = merge_jobs[idx];
            {
                std::lock_guard<std::mutex> lk(io_mtx);
                std::cout << "Merging " << job.inputs.size() << " ROOT files -> " << job.output << "\n";
            }

            std::vector<std::string> args{"hadd", "-f", job.output};
            args.insert(args.end(), job.inputs.begin(), job.inputs.end());
            const int rc = RunCommand(args);
            if (rc != 0) {
                std::lock_guard<std::mutex> lk(io_mtx);
                std::cerr << "  hadd failed with code " << rc << " for " << job.output << "\n";
                errors++;
                return;
            }

            for (const auto &remove_files : job.inputs) {
                std::error_code ec;
                if (!std::filesystem::remove(remove_files, ec) || ec) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "  failed to remove merged input " << remove_files;
                    if (ec)
                        std::cerr << ": " << ec.message();
                    std::cerr << "\n";
                    errors++;
                }
            }
        });
    }

    return errors > 0 ? 1 : 0;
}
