//=============================================================================
// ToolUtils.cpp — out-of-line helpers of the analysis command-line tools
//=============================================================================

#include "ToolUtils.h"
#include "Replay.h"

#include <TClass.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TTree.h>

#include <atomic>
#include <map>
#include <mutex>
#include <system_error>
#include <thread>

#include <spawn.h>
#include <sys/wait.h>

extern char **environ;

namespace analysis {

int RunCommand(const std::vector<std::string> &argv)
{
    if (argv.empty()) return 127;
    std::vector<char *> args;
    args.reserve(argv.size() + 1);
    for (const auto &a : argv) args.push_back(const_cast<char *>(a.c_str()));
    args.push_back(nullptr);

    // Our buffered output belongs before the child's in a redirected log.
    std::cout.flush();
    // posix_spawnp, unlike fork + execvp, is safe to use from a threaded
    // process and reports a failed exec back to the caller.
    pid_t pid = 0;
    const int err = posix_spawnp(&pid, args[0], nullptr, nullptr, args.data(), environ);
    if (err != 0) {
        std::cerr << argv[0] << ": cannot run: "
                  << std::generic_category().message(err) << "\n";
        return 127;
    }
    int status = 0;
    while (waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        std::cerr << argv[0] << ": waitpid failed: "
                  << std::generic_category().message(errno) << "\n";
        return 1;
    }
    return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}

void InitRootThreading()
{
    ROOT::EnableThreadSafety();
    for (const char *cls : {"TTree", "TFile", "TBranch", "TH1F", "TH2F"})
        TClass::GetClass(cls);
}

Long64_t TreeEntries(const std::string &path, const char *tree)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return 0;
    auto *t = dynamic_cast<TTree *>(f->Get(tree));
    return t ? t->GetEntries() : 0;
}

std::vector<Long64_t> DistributeEventBudget(const std::vector<std::string> &files,
                                            const char *tree, Long64_t max_events)
{
    std::vector<Long64_t> limits(files.size(), -1);
    if (max_events < 0) return limits;
    Long64_t remaining = max_events;
    for (size_t i = 0; i < files.size(); ++i) {
        limits[i] = remaining > 0 ? std::min(TreeEntries(files[i], tree), remaining) : 0;
        remaining -= limits[i];
    }
    return limits;
}

void RunFilesInRounds(const std::vector<std::string> &files, int nthreads,
                      const std::function<bool(int file_index, int slot)> &job,
                      const std::function<void(int first, int last)> &after_round)
{
    const int n_files = static_cast<int>(files.size());
    const int threads_count = std::max(1, std::min(nthreads, n_files));
    const int rounds = (n_files + threads_count - 1) / threads_count;
    std::cout << "Processing " << n_files << " file(s) with "
              << threads_count << " thread(s), " << rounds << " round(s)\n";

    std::mutex io_mutex;
    for (int round = 0; round < rounds; ++round) {
        const int first = round * threads_count;
        const int last = std::min(first + threads_count, n_files);
        std::vector<std::thread> workers;
        workers.reserve(last - first);
        for (int file_index = first; file_index < last; ++file_index) {
            workers.emplace_back([&, file_index, first]() {
                const bool ok = job(file_index, file_index - first);
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cout << "[worker " << (file_index - first) << "] file "
                          << file_index << " / " << (files.size() - 1)
                          << ": " << files[file_index] << " -> "
                          << (ok ? "OK" : "FAILED") << "\n";
            });
        }
        for (auto &worker : workers) worker.join();
        if (after_round) after_round(first, last);
    }
}

void ParallelFor(size_t n, int n_threads,
                 const std::function<void(size_t i, int worker)> &fn)
{
    if (n == 0) return;
    const size_t n_workers = std::min(n, static_cast<size_t>(std::max(n_threads, 1)));
    std::atomic<size_t> next{0};
    std::vector<std::thread> workers;
    workers.reserve(n_workers);
    for (size_t w = 0; w < n_workers; ++w)
        workers.emplace_back([&, w]() {
            for (size_t i = next++; i < n; i = next++) fn(i, static_cast<int>(w));
        });
    for (auto &worker : workers) worker.join();
}

int RunReplayPool(const std::vector<std::string> &inputs, int n_threads,
                  const std::string &daq_config, const std::string &daq_map,
                  const std::function<std::string(const std::string &input)> &output_for,
                  const std::function<bool(Replay &replay, const std::string &input,
                                           const std::string &output)> &run,
                  std::vector<char> *ok)
{
    const size_t n = inputs.size();
    if (ok) ok->assign(n, 0);

    // One Replay per worker thread (own EvChannel, own buffers).
    std::vector<std::unique_ptr<Replay>> replays(std::max(n_threads, 1));
    std::atomic<int> errors{0};
    std::mutex io_mtx;
    ParallelFor(n, n_threads, [&](size_t i, int worker) {
        auto &replay = replays[worker];
        if (!replay) {
            replay = std::make_unique<Replay>();
            if (!daq_config.empty()) replay->LoadDaqConfig(daq_config);
            replay->LoadHyCalMap(daq_map);
        }
        const std::string &input = inputs[i];
        const std::string output = output_for(input);
        const bool good = run(*replay, input, output);
        if (ok) (*ok)[i] = good;

        std::lock_guard<std::mutex> lk(io_mtx);
        if (good) {
            std::cout << "  [" << (i + 1) << "/" << n << "] "
                      << input << " -> " << output << "\n";
        } else {
            std::cerr << "  [" << (i + 1) << "/" << n << "] FAILED: " << input << "\n";
            ++errors;
        }
    });
    return errors.load();
}

bool MergeTopLevelHistograms(const std::vector<std::string> &inputs,
                             const std::string &output)
{
    TFile *merged = TFile::Open(output.c_str(), "RECREATE");
    if (!merged || merged->IsZombie()) {
        std::cerr << "Cannot create merged output " << output << "\n";
        delete merged;
        return false;
    }

    std::map<std::string, std::unique_ptr<TH1>> merged_histograms;
    for (const auto &path : inputs) {
        TFile *input = TFile::Open(path.c_str(), "READ");
        if (!input || input->IsZombie()) {
            delete input;
            continue;
        }
        TIter keys(input->GetListOfKeys());
        while (auto *key = dynamic_cast<TKey *>(keys())) {
            TObject *object = key->ReadObj();
            auto *hist = dynamic_cast<TH1 *>(object);
            if (!hist) {
                delete object;
                continue;
            }

            const std::string name = hist->GetName();
            auto it = merged_histograms.find(name);
            if (it == merged_histograms.end()) {
                std::unique_ptr<TH1> copy(dynamic_cast<TH1 *>(hist->Clone(name.c_str())));
                if (copy) {
                    copy->SetDirectory(nullptr);
                    merged_histograms.emplace(name, std::move(copy));
                }
            } else {
                it->second->Add(hist);
            }
            delete object;
        }
        input->Close();
        delete input;
    }

    merged->cd();
    for (auto &[name, hist] : merged_histograms) {
        hist->SetDirectory(merged);
        hist->Write(name.c_str(), TObject::kOverwrite);
        // Detach before the file goes away, or the histogram's destructor
        // would try to remove itself from a deleted directory.
        hist->SetDirectory(nullptr);
    }
    merged->Close();
    delete merged;
    return true;
}

} // namespace analysis
