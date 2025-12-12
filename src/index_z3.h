// index_z3.h
// Z3 Space-Filling Index for Aperture-3 Hexagonal Grids
//
// Encodes (i,j) coordinates as a hierarchical index for aperture-3 subdivision.
// The bijective mapping preserves locality in the hex grid.
//
// Mathematical basis: Each aperture-3 level subdivides a hex into 3 children
// arranged in a triangular pattern. The encoding maps base-3 digit pairs
// to child positions using the geometric relationship.
//
// References:
// - Sahr et al. (2003) "Geodesic Discrete Global Grid Systems"
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once

#include <string>

namespace hexify {
namespace z3 {

// Z3 Encoding: (i, j) coordinates → Z3 hierarchical string
// Resolution determines number of subdivision levels
std::string encode(long long i, long long j, int resolution);

// Z3 Decoding: Z3 hierarchical string → (i, j) coordinates
void decode(const std::string& z3_str, int resolution,
            long long& i, long long& j);

} // namespace z3
} // namespace hexify
