// index_zorder.h
// Z-Order (Morton) Curve Indexing for Hexagonal Grids
//
// Implements Morton's space-filling curve (1966) for multi-resolution indexing.
// Z-order interleaves coordinate digits to create a 1D index preserving 2D locality.
//
// For aperture N grids:
// - Convert (i, j) to radix-N representation
// - Interleave digits: i0 j0 i1 j1 i2 j2 ...
// - Result is a string encoding the hierarchical path
//
// References:
// - Morton, G.M. (1966) "A computer oriented geodetic data base"
// - Wikipedia "Z-order curve" for general theory
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once

#include <string>

namespace hexify {
namespace zorder {

// Z-Order for Aperture 3
std::string encode_ap3(long long i, long long j, int resolution);
void decode_ap3(const std::string& z_str, int resolution, 
                long long& i, long long& j);

// Z-Order for Aperture 4
std::string encode_ap4(long long i, long long j, int resolution);
void decode_ap4(const std::string& z_str, long long& i, long long& j);

// Z-Order for Aperture 7
std::string encode_ap7(long long i, long long j, int resolution);
void decode_ap7(const std::string& z_str, int resolution, 
                long long& i, long long& j);

} // namespace zorder
} // namespace hexify
