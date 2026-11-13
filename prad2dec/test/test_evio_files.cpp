//=============================================================================
// test_evio_files.cpp — run_number_from_path / discover_split_files.
//
// discover_split_files is exercised on a scratch directory created under the
// current working directory (the build tree when run by CTest).
//
// Build target: prad2dec_test_evio_files  (CTest: prad2dec.evio_files).
//=============================================================================

#include "EvioFiles.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {
int g_total    = 0;
int g_failures = 0;

void expect_run(const std::string &path, int want)
{
    ++g_total;
    int got = prad2::run_number_from_path(path);
    if (got != want) {
        ++g_failures;
        std::fprintf(stderr, "  FAIL  run_number_from_path(\"%s\") = %d, want %d\n",
                     path.c_str(), got, want);
    }
}

void expect_files(const std::string &path, const std::vector<std::string> &want,
                  bool expand_siblings = false)
{
    ++g_total;
    auto got = prad2::discover_split_files(path, expand_siblings);
    if (got != want) {
        ++g_failures;
        std::fprintf(stderr, "  FAIL  discover_split_files(\"%s\") returned %zu file(s):\n",
                     path.c_str(), got.size());
        for (auto &g : got) std::fprintf(stderr, "        %s\n", g.c_str());
    }
}
} // namespace

int main()
{
    namespace fs = std::filesystem;

    std::printf("run_number_from_path\n");
    expect_run("prad_023881.evio.00000", 23881);
    expect_run("/data/prad_023881.evio.00003", 23881);
    expect_run("run_1234", 1234);
    expect_run("RUN_0042.dat", 42);
    expect_run("/data/prad_023867/", 23867);
    expect_run("/data/prad_023867//", 23867);
    expect_run("/data/prad_023867", 23867);
    expect_run("/data/prad_023867/some_file.dat", 23867);
    expect_run("/data/prad_023867/prad_023881.evio.00000", 23881);   // name wins
    expect_run("/data/prad_023867/sub/some_file.dat", -1);           // parent only
    expect_run("/data/some_dir/", -1);
    expect_run("some_file.dat", -1);
    expect_run("", -1);
    expect_run("/", -1);

    std::printf("discover_split_files\n");
    const fs::path root = fs::absolute("evio_files_test_tmp");
    std::error_code ec;
    fs::remove_all(root, ec);
    const fs::path dir = root / "prad_023881";
    fs::create_directories(dir);
    std::vector<std::string> want;
    for (const char *name : {"prad_023881.evio.00000", "prad_023881.evio.00001"}) {
        std::ofstream(dir / name) << "x";
        want.push_back((dir / name).string());
    }
    std::ofstream(dir / "prad_023881.00000_raw.root") << "x";

    expect_files(dir.string(), want);
    expect_files(dir.string() + "/", want);
    expect_files(dir.string() + "//", want);
    expect_files((dir / "prad_023881.evio.*").string(), want);
    expect_files((dir / "prad_023881.evio.00001").string(), want, true);
    expect_files((dir / "prad_023881.00000_raw.root").string(),
                 {(dir / "prad_023881.00000_raw.root").string()});

    fs::remove_all(root, ec);

    std::printf("\n%d checks, %d failures\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}
