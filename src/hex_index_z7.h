// hex_index_z7.h
// Z7 Encoding/Decoding for Aperture 7 Hexagons
// Exact replication of DGGRID's DgZ7StringRF.cpp logic
// With canonical form support

#pragma once

#include "hex_index_ivec3d.h"
#include <string>

namespace hexify {
namespace z7 {

extern const int adjacentBaseCellTable[12][4];
extern const int inverseAdjacentBaseCellTable[12][2];

std::string encode(int quadNum, long long i, long long j, int resolution);

void decode(const std::string& z7_str, int resolution,
            int& quadNum, long long& i, long long& j);

// Get the canonical form of a Z7 index
// Finds the lexicographically smallest index in the cycle
// max_iterations: maximum number of decode/encode cycles to try (default 128)
std::string canonical_form(const std::string& z7_index, int max_iterations = 128);

} // namespace z7
} // namespace hexify
