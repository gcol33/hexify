// index_z7.h
// Z7 Hierarchical Index for Aperture-7 Hexagonal Grids
//
// Implements hierarchical space-filling indexing for aperture-7 hex subdivision.
// Each parent hex subdivides into 7 children (1 center + 6 surrounding).
// The index string encodes the traversal path from the base cell to the target.
//
// Mathematical basis:
// - Aperture-7 scaling rotates ~19.1° (arctan(sqrt(3)/5)) and scales by sqrt(7)
// - Coordinates use cube system (i,j,k) with constraint i+j+k=0
// - Resolution alternates Class II/III orientation
//
// References:
// - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
// - H3 Coordinate Systems (h3geo.org/docs/core-library/coordsystems)
// - Red Blob Games "Hexagonal Grids" (cube coordinates)
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once

#include "ijk_coordinates.h"
#include <string>

namespace hexify {
namespace z7 {

// Bijective aperture-7 hierarchical index (hexify-native). Keeps the quad fixed
// (no DGGRID base-cell reassignment / pentagon skip), so every (quad, i, j)
// round-trips. The leading field is quad + n * seed, n the solid's number of
// quads and seed the unit digit the hierarchy walk arrives at; it is the plain
// two-digit quad DGGRID writes
// for a cell whose whole ancestry lies inside its quad, and about two cells in
// three sit on a quad boundary and carry a nonzero seed instead. (i,j) are
// Class I substrate.
std::string encode_bijective(int quadNum, long long i, long long j, int resolution);

void decode_bijective(const std::string& index, int resolution,
                      int& quadNum, long long& i, long long& j);

// Get the canonical form of a Z7 index
// Finds the lexicographically smallest index in the cycle
// max_iterations: maximum number of decode/encode cycles to try (default 128)
std::string canonical_form(const std::string& z7_index, int max_iterations = 128);

} // namespace z7
} // namespace hexify
