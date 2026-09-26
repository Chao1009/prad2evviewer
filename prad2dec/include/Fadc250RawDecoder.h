#pragma once
//=============================================================================
// Fadc250RawDecoder.h — FADC250 hardware-format raw data decoder (0xE109)
//
// Decodes JLab FADC250 self-describing 32-bit word format.
// Used when rol2 composite reformatting is skipped.
// Outputs to the same RocData structures as Fadc250Decoder (composite).
//=============================================================================

#include "Fadc250Data.h"
#include <cstdint>
#include <cstddef>

namespace fdec
{

class Fadc250RawDecoder
{
public:
    // Decode one ROC's FADC250 hardware-format raw data bank (0xE109
    // payload words).  Returns number of slots decoded.
    static int DecodeRoc(const uint32_t *data, size_t nwords, RocData &roc);
};

} // namespace fdec
