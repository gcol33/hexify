// index_zorder.cpp
// Z-Order (Morton) Curve Indexing Implementation
//
// Encodes (i, j) coordinates using digit interleaving for space-filling indexing.
// The technique dates to Morton (1966) and is widely used in spatial databases.
//
// Algorithm:
// 1. Convert i and j to radix-N strings (N = aperture)
// 2. Interleave digits: result[2k] = i[k], result[2k+1] = j[k]
// 3. For aperture 4, combine into single digit: d = 2*i_bit + j_bit
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include "index_zorder.h"
#include <stdexcept>
#include <cmath>

namespace hexify {
namespace zorder {

// Helper: Convert integer to radix string
static std::string int_to_radix(long long value, int radix, int num_digits) {
  std::string result;
  long long v = std::abs(value);
  
  for (int i = 0; i < num_digits; i++) {
    result = std::to_string(v % radix) + result;
    v /= radix;
  }
  
  return result;
}

// Helper: Convert radix string to integer
static long long radix_to_int(const std::string& radix_str, int radix) {
  long long result = 0;
  for (char c : radix_str) {
    int digit = c - '0';
    if (digit < 0 || digit >= radix) {
      throw std::runtime_error("hex_index_zorder: invalid digit in radix string");
    }
    result = result * radix + digit;
  }
  return result;
}

// ============================================================================
// Aperture 4: Combine i,j bits into single digits
// ============================================================================

std::string encode_ap4(long long i, long long j, int resolution) {
  std::string result;
  
  // Convert to radix-2 strings
  std::string i_str = int_to_radix(i, 2, resolution);
  std::string j_str = int_to_radix(j, 2, resolution);
  
  // Interleave and combine: digit = i_bit * 2 + j_bit
  for (int r = 0; r < resolution; r++) {
    int i_bit = i_str[r] - '0';
    int j_bit = j_str[r] - '0';
    int digit = i_bit * 2 + j_bit;
    result += std::to_string(digit);
  }
  
  return result;
}

void decode_ap4(const std::string& z_str, long long& i, long long& j) {
  i = 0;
  j = 0;
  
  for (char c : z_str) {
    int digit = c - '0';
    int i_bit = digit / 2;
    int j_bit = digit % 2;
    
    i = i * 2 + i_bit;
    j = j * 2 + j_bit;
  }
}

// ============================================================================
// Aperture 3: Radix-3 digit interleaving
// ============================================================================

std::string encode_ap3(long long i, long long j, int resolution) {
  int eff_res = (resolution + 1) / 2;
  bool is_class_i = (resolution % 2 == 0);
  
  std::string result;
  if (eff_res == 0) return result;
  
  std::string i_str = int_to_radix(i, 3, eff_res);
  std::string j_str = int_to_radix(j, 3, eff_res);
  
  // Simple alternation (i then j)
  for (int idx = 0; idx < eff_res; idx++) {
    result += i_str[idx];
    result += j_str[idx];
  }
  
  // Trim last digit if Class II
  if (!is_class_i && result.length() > 0) {
    result.pop_back();
  }
  
  return result;
}

void decode_ap3(const std::string& z_str, int resolution, 
                long long& i, long long& j) {
  bool is_class_i = (resolution % 2 == 0);
  std::string adjusted = z_str;
  
  // For Class II: infer the missing j-digit from last i-digit
  if (!is_class_i && adjusted.length() % 2 == 1) {
    // Get the last i-digit (which is the last character in the string)
    int last_i_digit = adjusted[adjusted.length() - 1] - '0';
    
    // Infer j-digit using the constraint: j%3 == jDigits[i%3]
    // jDigits mapping: 0→0, 1→2, 2→1
    std::string j_digits[] = {"0", "2", "1"};
    adjusted += j_digits[last_i_digit % 3];
  }
  
  // Split alternating digits (i then j)
  std::string i_str, j_str;
  
  for (size_t idx = 0; idx < adjusted.length(); idx += 2) {
    i_str += adjusted[idx];
    if (idx + 1 < adjusted.length()) {
      j_str += adjusted[idx + 1];
    }
  }
  
  i = radix_to_int(i_str, 3);
  j = radix_to_int(j_str, 3);
}

// ============================================================================
// Aperture 7: Simple radix-7 alternation
// ============================================================================

std::string encode_ap7(long long i, long long j, int resolution) {
  if (resolution == 0) return "";
  
  // Convert to radix-7 strings (one digit per resolution level)
  std::string i_str = int_to_radix(i, 7, resolution);
  std::string j_str = int_to_radix(j, 7, resolution);
  
  // Interleave: i0 j0 i1 j1 i2 j2 ...
  std::string result;
  for (int idx = 0; idx < resolution; idx++) {
    result += i_str[idx];
    result += j_str[idx];
  }
  
  return result;
}

void decode_ap7(const std::string& z_str, int resolution, 
                long long& i, long long& j) {
  if (z_str.empty()) {
    i = 0;
    j = 0;
    return;
  }
  
  // De-interleave: i digits are at even positions, j at odd
  std::string i_str, j_str;
  for (size_t idx = 0; idx < z_str.length(); idx += 2) {
    i_str += z_str[idx];
    if (idx + 1 < z_str.length()) {
      j_str += z_str[idx + 1];
    }
  }
  
  // Convert radix-7 strings to integers
  i = radix_to_int(i_str, 7);
  j = radix_to_int(j_str, 7);
}

} // namespace zorder
} // namespace hexify
