//=============================================================================
// EpicsData.cpp — text → (channel, value) parser for EPICS slow-control banks.
//=============================================================================

#include "EpicsData.h"

#include <cstdlib>
#include <sstream>

namespace epics
{

int ParseEpicsText(const std::string &text, EpicsRecord &out)
{
    out.channel.clear();
    out.value.clear();
    if (text.empty()) return 0;

    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
        // strip trailing \r and surrounding whitespace
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' ||
                                 line.back() == '\t' || line.back() == '\n'))
            line.pop_back();
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;

        // Lines look like "<value>  <channel_name>".  We accept whitespace
        // separators of any length.  Strict format check: must have exactly
        // one numeric leading field.
        size_t mid = line.find_first_of(" \t", start);
        if (mid == std::string::npos) continue;
        std::string val_str  = line.substr(start, mid - start);
        size_t name_start    = line.find_first_not_of(" \t", mid);
        if (name_start == std::string::npos) continue;
        std::string name_str = line.substr(name_start);

        char *endp = nullptr;
        double v = std::strtod(val_str.c_str(), &endp);
        if (endp == val_str.c_str()) continue;          // not a number
        if (name_str.empty())        continue;

        out.channel.push_back(std::move(name_str));
        out.value.push_back(v);
    }
    return static_cast<int>(out.channel.size());
}

} // namespace epics
