#pragma once
//=============================================================================
// Fadc250Decoder.h — decode composite FADC250 data into pre-allocated buffers
//=============================================================================

#include "Fadc250Data.h"
#include <cstddef>

namespace fdec
{

class Fadc250Decoder
{
public:
    // Decode one ROC's composite payload into roc.
    // Format "c,i,l,N(c,Ns)" — packed native-endian.
    // Returns number of slots decoded, or -1 on error.
    static int DecodeRoc(const uint8_t *data, size_t nbytes, RocData &roc);
};

} // namespace fdec
