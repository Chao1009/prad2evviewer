// cosmic_common.h: file readers shared by the cosmic gain macros.
#pragma once

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

// One module of a cosmic_modules_run<N>.json written by prad2ana_cosmic_test.
// Only the first entry of the module's array is read, which in a file that
// cosmic_test -j appended runs to is the oldest run.
struct CosmicEntry {
    std::string name;                   // "W12", "G156"
    int run;
    double ph_mean, ph_sigma, ph_diff;  // peak height
    double pi_mean, pi_sigma, pi_diff;  // peak integral
    int count;
};

// Appends the W/G module entries of path in file order; false when the file
// cannot be opened.
inline bool read_cosmic_modules(const std::string &path, std::vector<CosmicEntry> &entries)
{
    std::ifstream fin(path);
    if (!fin.is_open()) return false;
    std::string line;
    while (std::getline(fin, line)) {
        CosmicEntry e;
        char key = 0;
        int number = 0;
        if (std::sscanf(line.c_str(),
                        " \"%c%d\": [{\"run\": %d, \"peak_height_mean\": %lf, \"peak_height_sigma\": %lf, \"peak_height_diff\": %lf,"
                        " \"peak_integral_mean\": %lf, \"peak_integral_sigma\": %lf, \"peak_integral_diff\": %lf,"
                        " \"count\": %d",
                        &key, &number, &e.run, &e.ph_mean, &e.ph_sigma, &e.ph_diff,
                        &e.pi_mean, &e.pi_sigma, &e.pi_diff, &e.count) != 10
            || (key != 'W' && key != 'G'))
            continue;
        e.name = key + std::to_string(number);
        entries.push_back(e);
    }
    return true;
}

// Follows a vset_iter<N>.json line by line: name is the channel ("W12",
// "G156") of the last "name" line seen, empty before the first one.
struct VsetLineTracker {
    std::string name;

    // Advances over line; true when it is the current channel's "V0Set" line.
    bool isV0Set(const std::string &line)
    {
        char buf[64];
        if (line.find("\"name\"") != std::string::npos
            && std::sscanf(line.c_str(), " \"name\": \"%63[^\"]\"", buf) == 1)
            name = buf;
        return !name.empty() && line.find("\"V0Set\"") != std::string::npos;
    }
};
