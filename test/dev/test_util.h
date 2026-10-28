#pragma once
//=============================================================================
// test_util.h — minimal check / summary harness for the ctest unit tests in
// test/dev (each test is a single-TU executable, so this stays header-only).
//=============================================================================

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <unistd.h>

namespace testutil {

inline int failures = 0;

inline void check(bool condition, const std::string &message)
{
    if (condition) return;
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
}

inline bool close_to(float actual, float expected, float tolerance = 1.e-4f)
{
    return std::fabs(actual - expected) <= tolerance;
}

// Creates <system temp>/<prefix><pid> for fixture files; the caller removes it.
inline std::filesystem::path make_temp_dir(const std::string &prefix)
{
    const auto dir = std::filesystem::temp_directory_path() /
        (prefix + std::to_string(static_cast<long long>(::getpid())));
    std::filesystem::create_directories(dir);
    return dir;
}

// Prints "<n> <fail_label> test(s) failed" or "<pass_label> tests passed" and
// returns the process exit code.
inline int finish(const char *fail_label, const char *pass_label)
{
    if (failures != 0) {
        std::cerr << failures << ' ' << fail_label << " test(s) failed\n";
        return 1;
    }
    std::cout << pass_label << " tests passed\n";
    return 0;
}

} // namespace testutil
