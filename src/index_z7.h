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

// Aperture-7 hierarchical index of a cell given as its quad and Class I
// substrate (i, j). On the icosahedron it is IGEO7's Z7 index, identical to
// DGGRID's: a two-digit base cell followed by one digit per resolution. On
// other solids the leading field is quad + n * seed, n the solid's number of
// quads and seed the unit digit the hierarchy walk arrives at.
std::string encode(int quadNum, long long i, long long j, int resolution);

// The inverse of encode(): the quad and Class I substrate (i, j) of an index.
// The resolution is the number of digits after the leading field.
void decode(const std::string& index, int& quadNum, long long& i, long long& j);

// Whether an IGEO7 index lies in the subsequence a pentagon deletes: its first
// nonzero digit is the direction its base cell lacks. Such a string names no
// cell. Always false on solids without IGEO7 labels.
bool in_deleted_subsequence(const std::string& index);

// Get the canonical form of a Z7 index
// Finds the lexicographically smallest index in the cycle
// max_iterations: maximum number of decode/encode cycles to try (default 128)
std::string canonical_form(const std::string& z7_index, int max_iterations = 128);

} // namespace z7
} // namespace hexify
