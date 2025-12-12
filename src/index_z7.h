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

// Base cell adjacency tables derived from icosahedral topology.
// The icosahedron has 12 vertices; each base cell corresponds to one vertex.
// adjacentBaseCellTable[cell][dir] gives the adjacent cell in direction dir.
// Directions: 0=self, 1=edge1, 2=edge2, 3=edge3 (varies by hemisphere)
extern const int adjacentBaseCellTable[12][4];

// Inverse adjacency for decoding: maps encoded cell back to original quad
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
