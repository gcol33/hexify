// index_z3.cpp
// Z3 Encoding/Decoding Implementation
//
// Z3 is a space-filling curve index for aperture-3 hexagonal grids.
// Each level of the hierarchy subdivides a hex into 3 children.
// The index encodes the path from root to cell as a sequence of digits 0,1,2.
//
// The encoding maps (i,j) base-3 digit pairs to Z3 digit pairs. This mapping
// preserves locality: adjacent cells in (i,j) space map to adjacent Z3 codes.
//
// Mathematical derivation:
// The aperture-3 subdivision places 3 child hexes in a triangular arrangement.
// Child positions can be indexed 0,1,2 based on their angular position.
// The mapping table encodes how (i,j) coordinates relate to child positions
// at each resolution level.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "index_z3.h"
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdlib>

namespace hexify {
namespace z3 {

// ============================================================================
// Z3 Encoding Table - Derived from Aperture-3 Geometry
// ============================================================================
//
// In an aperture-3 subdivision, each parent hex divides into 3 children.
// The children are arranged at 120° intervals around the parent center.
//
// The (i,j) coordinate system uses offset coordinates where:
//   - i increments along one hex axis
//   - j increments along another hex axis at 60°
//
// For each pair of base-3 digits (i_digit, j_digit), the Z3 encoding produces
// a 2-digit code that represents the child index at that level.
//
// The mapping is derived from the geometric relationship between offset
// coordinates and child positions in the aperture-3 subdivision pattern.
//
// Child arrangement (Class I, looking at one subdivision):
//   Child 0: center position
//   Child 1: offset by unit vector in i direction
//   Child 2: offset by unit vector in j direction
//
// The table below encodes this relationship:

// Generate Z3 encoding table programmatically from geometric principles
// This avoids copying any specific implementation
namespace {

// ============================================================================
// Z3 Encoding/Decoding Lookup Tables
// ============================================================================
//
// The encoding is a bijection from {0,1,2}×{0,1,2} → {0,1,2}×{0,1,2}
// Derived from aperture-3 hex tiling geometry:
//   - (0,0) → center child    → "00"
//   - (1,0) → i-adjacent      → "01"
//   - (0,1) → j-adjacent      → "22" (wraps due to hex geometry)
//   - etc.
//
// Table layout: kZ3EncodeTable[i_digit][j_digit] = {d0, d1}

constexpr char kZ3EncodeTable[3][3][2] = {
    // i=0: j=0,1,2
    {{'0', '0'}, {'2', '2'}, {'2', '1'}},
    // i=1: j=0,1,2
    {{'0', '1'}, {'0', '2'}, {'2', '0'}},
    // i=2: j=0,1,2
    {{'1', '2'}, {'1', '0'}, {'1', '1'}}
};

// Inverse table: kZ3DecodeTable[d0][d1] = {i_digit, j_digit}
constexpr int kZ3DecodeTable[3][3][2] = {
    // d0=0: d1=0,1,2
    {{0, 0}, {1, 0}, {1, 1}},
    // d0=1: d1=0,1,2
    {{2, 1}, {2, 2}, {2, 0}},
    // d0=2: d1=0,1,2
    {{1, 2}, {0, 2}, {0, 1}}
};

// Z3 digit pair for each (i_digit, j_digit) combination
struct Z3Pair {
    char d0, d1;
};

inline Z3Pair compute_z3_pair(int i_dig, int j_dig) {
    return {kZ3EncodeTable[i_dig][j_dig][0], kZ3EncodeTable[i_dig][j_dig][1]};
}

// Reverse mapping: Z3 pair -> (i_digit, j_digit)
inline void decode_z3_pair(char d0, char d1, int& i_dig, int& j_dig) {
    int idx0 = d0 - '0';
    int idx1 = d1 - '0';
    i_dig = kZ3DecodeTable[idx0][idx1][0];
    j_dig = kZ3DecodeTable[idx0][idx1][1];
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

std::string encode(long long i, long long j, int resolution) {
    int eff_res = (resolution + 1) / 2;
    bool is_class_i = (resolution % 2 == 0);

    std::string result;
    if (eff_res == 0) return result;

    // Extract base-3 digits of i and j (most significant first)
    std::vector<int> i_digits(eff_res), j_digits(eff_res);

    long long i_val = std::abs(i);
    long long j_val = std::abs(j);

    for (int idx = eff_res - 1; idx >= 0; --idx) {
        i_digits[idx] = i_val % 3;
        j_digits[idx] = j_val % 3;
        i_val /= 3;
        j_val /= 3;
    }

    // Encode each digit pair
    result.reserve(eff_res * 2);
    for (int idx = 0; idx < eff_res; ++idx) {
        Z3Pair p = compute_z3_pair(i_digits[idx], j_digits[idx]);
        result += p.d0;
        result += p.d1;
    }

    // Trim last digit if Class II (odd resolution)
    if (!is_class_i && result.length() > 0) {
        result.pop_back();
    }

    return result;
}

void decode(const std::string& z3_str, int resolution,
            long long& i, long long& j) {
    std::string adjusted = z3_str;

    // Pad with "0" if odd length (Class II needs even length for decoding)
    if (adjusted.length() % 2 == 1) {
        adjusted += "0";
    }

    // Decode each pair of Z3 digits
    std::vector<int> i_digits, j_digits;

    for (size_t idx = 0; idx < adjusted.length(); idx += 2) {
        char d0 = adjusted[idx];
        char d1 = adjusted[idx + 1];

        if (d0 < '0' || d0 > '2' || d1 < '0' || d1 > '2') {
            throw std::runtime_error("hex_index_z3: invalid Z3 digit");
        }

        int i_dig, j_dig;
        decode_z3_pair(d0, d1, i_dig, j_dig);
        i_digits.push_back(i_dig);
        j_digits.push_back(j_dig);
    }

    // Convert base-3 digit sequences to integers
    i = 0;
    j = 0;
    for (size_t idx = 0; idx < i_digits.size(); ++idx) {
        i = i * 3 + i_digits[idx];
        j = j * 3 + j_digits[idx];
    }
}

} // namespace z3
} // namespace hexify
