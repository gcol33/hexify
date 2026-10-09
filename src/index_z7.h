// index_z7.h
// Z7 Hierarchical Index for Aperture-7 Hexagonal Grids
//
// Implements hierarchical space-filling indexing for aperture-7 hex subdivision.
// Each parent hex subdivides into 7 children (1 center + 6 surrounding).
// The index encodes the traversal path from the base cell to the target.
//
// Mathematical basis:
// - Aperture-7 scaling rotates ~19.1° (arctan(sqrt(3)/5)) and scales by sqrt(7)
// - Coordinates use cube system (i,j,k) with constraint i+j+k=0
// - Resolution alternates Class II/III orientation
//
// References:
// - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
// - Sahr (2019) "Central Place Indexing", Cartographica 54(1): 16-29
// - Kmoch, Sahr, Chan, Uuemaa (2025) "IGEO7", AGILE: GIScience Series 6
// - H3 Coordinate Systems (h3geo.org/docs/core-library/coordsystems)
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#pragma once

#include "ijk_coordinates.h"
#include "constants.h"
#include <cstdint>
#include <string>

namespace hexify {
namespace z7 {

typedef IVec3D::Direction Digit;

// IGEO7's packed form holds 20 three-bit digits below a four-bit base cell.
constexpr int kMaxPackedRes = 20;

// The deepest resolution whose monotonic IDs fit in a signed 64-bit integer:
// 12 * 7^21 < 2^63 < 12 * 7^22.
constexpr int kMaxMonotonicRes = 21;

// A cell's aperture-7 label: the leading field and one digit 0-6 per
// resolution. On the icosahedron the leading field is the IGEO7 base cell
// (0-11), the pentagon the cell descends from, and the digits are read in that
// base cell's frame. On other solids it is quad + n * seed, n the solid's
// number of quads and seed the unit digit the hierarchy walk arrives at.
struct Label {
    int lead;
    int res;
    Digit digit[kMaxResolution + 2];  // digit[1..res]
};

// Whether the active solid carries IGEO7 labels (the icosahedron).
bool igeo7_labels();

// The label of a cell given as its quad and Class I substrate (i, j).
Label label_of(int quadNum, long long i, long long j, int resolution);

// The quad and Class I substrate (i, j) a label names.
void cell_of(const Label& label, int& quadNum, long long& i, long long& j);

// Z7 string: the leading field as two decimal digits, then the digits.
std::string to_string(const Label& label);
Label from_string(const std::string& index);

// IGEO7's 64-bit integer: base cell in bits 63-60, digit r in bits
// 59 - 3(r - 1) to 57 - 3(r - 1), digit 7 at every level below the cell's
// resolution. Defined on the icosahedron up to resolution 20.
uint64_t to_packed(const Label& label);

// The label of a packed index; false unless `z` is one: base cell 0-11 and a
// run of digits 0-6 followed only by 7s.
bool from_packed(uint64_t z, Label& out);

// The monotonic ID: base * 7^r plus the digits read as a base-7 number, up to
// kMaxMonotonicRes. It counts the deleted subsequences of the pentagons too,
// so within one resolution it is increasing in the Z7 order and leaves holes.
uint64_t to_monotonic(const Label& label);

// The label of a monotonic ID at a resolution; false if its base passes 11.
bool from_monotonic(uint64_t m, int resolution, Label& out);

// The neighbour one step from an IGEO7 label in direction `dir` (1-6) of its
// base cell's frame, by digit arithmetic: false for a pentagon's deleted
// direction, where it has no neighbour.
bool neighbor(const Label& in, Digit dir, Label& out);

// Aperture-7 hierarchical index string of a cell given as its quad and Class I
// substrate (i, j): to_string(label_of(...)).
std::string encode(int quadNum, long long i, long long j, int resolution);

// The inverse of encode(): the quad and Class I substrate (i, j) of an index.
// The resolution is the number of digits after the leading field.
void decode(const std::string& index, int& quadNum, long long& i, long long& j);

// Whether an IGEO7 label lies in the subsequence a pentagon deletes: its first
// nonzero digit is the direction its base cell lacks. Such a label names no
// cell. Always false on solids without IGEO7 labels.
bool in_deleted_subsequence(const Label& label);
bool in_deleted_subsequence(const std::string& index);

// Get the canonical form of a Z7 index
// Finds the lexicographically smallest index in the cycle
// max_iterations: maximum number of decode/encode cycles to try (default 128)
std::string canonical_form(const std::string& z7_index, int max_iterations = 128);

} // namespace z7
} // namespace hexify
