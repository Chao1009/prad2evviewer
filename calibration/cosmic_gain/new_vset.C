#include "cosmic_common.h"

#include <map>

void new_vset(){

    std::string vset_file = "vset_iter1.json";
    std::string cosmic_file = "cosmic_modules_run575.json";
    std::string output_file = "vset_iter2.json";

    // ---- Read vset_file: V0Set per channel ----
    std::map<std::string, double> vset;
    {
        std::ifstream fin(vset_file);
        if (!fin.is_open()) { std::cerr << "Cannot open " << vset_file << std::endl; return; }
        std::string line;
        VsetLineTracker channel;
        while (std::getline(fin, line)) {
            double val;
            if (channel.isV0Set(line) && sscanf(line.c_str(), " \"V0Set\": %lf", &val) == 1)
                vset[channel.name] = val;
        }
    }

    // ---- Read cosmic_file: peak_height_mean ----
    std::vector<CosmicEntry> cosmic;
    if (!read_cosmic_modules(cosmic_file, cosmic)) { std::cerr << "Cannot open " << cosmic_file << std::endl; return; }

    // ---- Compute the new V0Set ----
    std::map<std::string, double> vsetnew;

    int n_valid = 0, n_skip = 0;
    int n_increase[3] = {0}, n_decrease[3] = {0}, n_unchanged = 0;

    for (const auto &c : cosmic) {
        auto it = vset.find(c.name);
        if (it == vset.end() || c.name == "W1019" || c.name == "W1020") { n_skip++; continue; }
        const double v = it->second, mean = c.ph_mean;
        double vn = v;

        if (mean > 45.0) {
            vn = v - 20.0;
            n_decrease[2]++;
        }else if (mean > 40.0) {
            vn = v - 10.0;
            n_decrease[1]++;
        }else if (mean > 37.0) {
            vn = v - 5.0;
            n_decrease[0]++;
        }

        if(mean < 25.0) {
            vn = v + 20.0;
            n_increase[2]++;
        }else if(mean < 30.0) {
            vn = v + 10.0;
            n_increase[1]++;
        }else if(mean < 33.0) {
            vn = v + 5.0;
            n_increase[0]++;
        }

        if(mean >= 33.0 && mean <= 37.0) {
            vn = v;
            n_unchanged++;
        }
        if (vn > 1270.0 && c.name[0] == 'W') vn = 1270.0;
        if (vn > 1800 && c.name[0] == 'G') vn = 1800.0;
        vsetnew[c.name] = vn;
        n_valid++;
    }

    printf("\n=== Summary: %d valid, %d skipped ===\n", n_valid, n_skip);
    printf("Voltage increased (+20V, height<25):  %d\n", n_increase[2]);
    printf("Voltage increased (+10V, height<30):  %d\n", n_increase[1]);
    printf("Voltage increased (+5V, height<33):  %d\n", n_increase[0]);
    printf("Voltage decreased (-20V, height>45):  %d\n", n_decrease[2]);
    printf("Voltage decreased (-10V, height>40):  %d\n", n_decrease[1]);
    printf("Voltage decreased (-5V, height>37):  %d\n", n_decrease[0]);
    printf("Voltage unchanged (33~37):            %d\n\n", n_unchanged);

    printf("%-6s %10s %10s %10s %8s\n", "Ch", "vset", "height", "vsetnew", "action");
    printf("----------------------------------------------------------\n");
    for (const auto &c : cosmic) {
        auto it = vsetnew.find(c.name);
        if (it == vsetnew.end()) continue;
        const char* action = "keep";
        double dv = it->second - vset[c.name];
        if (dv > 15.0) action = "+20V";
        else if (dv < -15.0) action = "-20V";
        else if (dv > 7.0) action = "+10V";
        else if (dv < -7.0) action = "-10V";
        else if (dv > 2.0) action = "+5V";
        else if (dv < -2.0) action = "-5V";
        printf("%-6s %10.1f %10.2f %10.1f %8s\n",
               c.name.c_str(), vset[c.name], c.ph_mean, it->second, action);
    }

    // ---- Write output_file: vset_file with V0Set replaced by vsetnew ----
    {
        std::ifstream fin(vset_file);
        std::ofstream fout(output_file);
        if (!fin.is_open() || !fout.is_open()) {
            std::cerr << "Cannot open files for writing" << std::endl;
            return;
        }
        std::string line;
        VsetLineTracker channel;
        while (std::getline(fin, line)) {
            auto it = channel.isV0Set(line) ? vsetnew.find(channel.name) : vsetnew.end();
            if (it != vsetnew.end()) {
                bool has_comma = (line.find(",") != std::string::npos &&
                                  line.rfind(",") > line.find("V0Set"));
                char buf[256];
                if (has_comma)
                    snprintf(buf, sizeof(buf), "        \"V0Set\": %.1f,", it->second);
                else
                    snprintf(buf, sizeof(buf), "        \"V0Set\": %.1f", it->second);
                fout << buf << "\n";
                continue;
            }
            fout << line << "\n";
        }
        printf("\n%s written.\n", output_file.c_str());
    }
}
