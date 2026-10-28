#pragma once
//=============================================================================
// ToolUtils.h — small helpers shared by the analysis command-line tools
//=============================================================================

#include <TH1.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace analysis {

class Replay;

// Strict option parsers: the whole argument must be a number in range.
// On failure `value` is left untouched and false is returned.
inline bool ParseIntOption(const char *text, int &value)
{
    if (!text || *text == '\0') return false;
    char *end = nullptr;
    errno = 0;
    const long parsed = std::strtol(text, &end, 10);
    if (errno == ERANGE || *end != '\0'
            || parsed < std::numeric_limits<int>::min()
            || parsed > std::numeric_limits<int>::max())
        return false;
    value = static_cast<int>(parsed);
    return true;
}

inline bool ParseFloatOption(const char *text, float &value)
{
    if (!text || *text == '\0') return false;
    char *end = nullptr;
    errno = 0;
    const float parsed = std::strtof(text, &end);
    if (errno == ERANGE || *end != '\0' || !std::isfinite(parsed)) return false;
    value = parsed;
    return true;
}

// getopt error reporters; both return the exit code 2.
inline int InvalidOptionValue(const char *flag, const char *value)
{
    std::cerr << "Invalid value for " << flag << ": "
              << (value ? value : "(missing)") << "\n";
    return 2;
}

// Pass getopt's optind and optopt after it returned '?'.
inline int InvalidOption(char *argv[], int opt_index, int opt_char)
{
    if (opt_char != 0)
        std::cerr << "Invalid option or missing argument: -"
                  << static_cast<char>(opt_char) << "\n";
    else
        std::cerr << "Invalid option: "
                  << (opt_index > 0 ? argv[opt_index - 1] : "(unknown)") << "\n";
    return 2;
}

// ── Input files ─────────────────────────────────────────────────────────────

// File-name filters for directory inputs (they see the name without its
// directory).  EVIO split files are *.evio.NNNNN, hence the substring test.
inline bool IsEvioName(const std::string &name)
{
    return name.find(".evio") != std::string::npos;
}

inline bool IsRootName(const std::string &name)
{
    return name.find(".root") != std::string::npos;
}

inline bool IsRawRootName(const std::string &name)
{
    return name.find("_raw.root") != std::string::npos;
}

// Also matches the merged prad_NNNNNN_recon_NNN.root files of replay_recon.
inline bool IsReconRootName(const std::string &name)
{
    return name.find("_recon") != std::string::npos
        && name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0;
}

inline bool IsLmsRootName(const std::string &name)
{
    return name.find("_lms.root") != std::string::npos;
}

// A directory expands to its regular files whose name passes keep(name),
// sorted by path; any other path is returned as is.
template <class Keep>
std::vector<std::string> ExpandInputPath(const std::string &path, Keep keep)
{
    namespace fs = std::filesystem;
    if (!fs::is_directory(path)) return {path};
    std::vector<std::string> files;
    for (const auto &entry : fs::directory_iterator(path))
        if (entry.is_regular_file() && keep(entry.path().filename().string()))
            files.push_back(entry.path().string());
    std::sort(files.begin(), files.end());
    return files;
}

// Expand argv[first, argc) in order.  With max_files > 0 at most that many
// files are returned and no argument is expanded once the cap is reached.
template <class Keep>
std::vector<std::string> CollectInputs(int argc, char *argv[], int first, Keep keep,
                                       int max_files = -1)
{
    const auto full = [&](const std::vector<std::string> &files) {
        return max_files > 0 && static_cast<int>(files.size()) >= max_files;
    };
    std::vector<std::string> files;
    for (int i = first; i < argc && !full(files); ++i)
        for (auto &f : ExpandInputPath(argv[i], keep)) {
            if (full(files)) break;
            files.push_back(std::move(f));
        }
    return files;
}

// ── Subprocesses ────────────────────────────────────────────────────────────

// Run argv[0] (looked up in PATH unless it contains a '/') with arguments
// argv[1..], without a shell, and wait for it.  Returns its exit status; 127
// when it could not be started, 1 when it was killed by a signal or could not
// be waited for.  Safe to call from several threads at once.
int RunCommand(const std::vector<std::string> &argv);

// ── ROOT and worker threads ─────────────────────────────────────────────────

// Call at the start of main(), before any thread touches ROOT.  Also loads
// the dictionaries of the I/O and histogram classes the workers use, so they
// are not first initialised concurrently.
void InitRootThreading();

// Entries of TTree `tree` in the ROOT file at path; 0 when the file cannot be
// opened or has no such tree.
Long64_t TreeEntries(const std::string &path, const char *tree);

// Per-file entry limits spending a total of max_events over the files in
// order: each file gets min(its `tree` entries, what is left), files after
// the budget is used up get 0 and are not opened (so max_events == 0 gives
// every file 0).  max_events < 0 means no limit: every file gets -1.
std::vector<Long64_t> DistributeEventBudget(const std::vector<std::string> &files,
                                            const char *tree, Long64_t max_events);

// Process the files in rounds of min(nthreads, N) threads, one file each.
// job(file_index, slot) runs on a worker thread, slot being the thread's
// position within its round, and returns success.  after_round(first, last)
// runs on the calling thread once the threads of files [first, last) have
// joined.  Prints the "Processing N file(s) ..." header and one
// "[worker slot] file i / N-1: path -> OK|FAILED" line per file.
void RunFilesInRounds(const std::vector<std::string> &files, int nthreads,
                      const std::function<bool(int file_index, int slot)> &job,
                      const std::function<void(int first, int last)> &after_round = {});

// Call fn(i, worker) for every i in [0, n) on max(1, min(n_threads, n))
// threads; each thread (worker = 0, 1, ...) takes the next unclaimed index
// until none is left.  Returns once every call has finished.
void ParallelFor(size_t n, int n_threads,
                 const std::function<void(size_t i, int worker)> &fn);

// Replay every input on up to n_threads threads, each with its own Replay set
// up by LoadDaqConfig(daq_config) (when non-empty) and LoadHyCalMap(daq_map).
// Calls run(replay, input, output_for(input)) and logs "  [i/N] input ->
// output" to stdout or "  [i/N] FAILED: input" to stderr.  Returns the number
// of failed inputs; ok, when given, receives a success flag per input.
int RunReplayPool(const std::vector<std::string> &inputs, int n_threads,
                  const std::string &daq_config, const std::string &daq_map,
                  const std::function<std::string(const std::string &input)> &output_for,
                  const std::function<bool(Replay &replay, const std::string &input,
                                           const std::string &output)> &run,
                  std::vector<char> *ok = nullptr);

// ── Histogram bundles ───────────────────────────────────────────────────────

// Histograms registered in booking order, so bundles booked by the same code
// (one per thread plus the merged one) can be merged pairwise.
using HistList = std::vector<TH1 *>;

// Create an H(args...) detached from any directory and register it in reg.
template <class H, class... Args>
std::unique_ptr<H> Book(HistList &reg, Args &&...args)
{
    auto h = std::make_unique<H>(std::forward<Args>(args)...);
    h->SetDirectory(nullptr);
    reg.push_back(h.get());
    return h;
}

// dst[i]->Add(src[i]) for every registered histogram.
inline void AddAll(const HistList &dst, const HistList &src)
{
    if (dst.size() != src.size())
        throw std::logic_error("AddAll: histogram bundles were booked differently");
    for (size_t i = 0; i < dst.size(); ++i) dst[i]->Add(src[i]);
}

// Write to the new ROOT file `output` the sum, by name, of the top-level
// histograms of the files `inputs` (unreadable inputs, other objects and
// subdirectories are skipped).  Returns false when output cannot be created.
bool MergeTopLevelHistograms(const std::vector<std::string> &inputs,
                             const std::string &output);

} // namespace analysis
