// hex_index_z3.h
// Z3 Encoding/Decoding for Aperture 3 Hexagons
// Based on DGGRID's DgZ3StringRF.cpp

#pragma once

#include <string>

namespace hexify {
namespace z3 {

// Z3 Encoding: (i,j) → Z3 String
// Based on DGGRID DgQ2DItoZ3StringConverter::convertTypedAddress
std::string encode(long long i, long long j, int resolution);

// Z3 Decoding: Z3 String → (i,j)
// Based on DGGRID DgZ3StringToQ2DIConverter::convertTypedAddress
void decode(const std::string& z3_str, int resolution,
            long long& i, long long& j);

} // namespace z3
} // namespace hexify
