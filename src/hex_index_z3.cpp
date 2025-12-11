// hex_index_z3.cpp
// Z3 Encoding/Decoding Implementation
// Based on DGGRID's DgZ3StringRF.cpp

#include "hex_index_z3.h"
#include <stdexcept>
#include <sstream>

namespace hexify {
namespace z3 {

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
      throw std::runtime_error("hex_index_z3: invalid digit in radix string");
    }
    result = result * radix + digit;
  }
  return result;
}

// Z3 lookup table from DGGRID (DgZ3StringRF.cpp lines 135-140)
static const std::string encode_table[3][3] = {
  {"00", "22", "21"}, // i=0, j=0,1,2
  {"01", "02", "20"}, // i=1, j=0,1,2
  {"12", "10", "11"}  // i=2, j=0,1,2
};

std::string encode(long long i, long long j, int resolution) {
  int eff_res = (resolution + 1) / 2;
  bool is_class_i = (resolution % 2 == 0);
  
  std::string result;
  if (eff_res == 0) return result;
  
  std::string i_str = int_to_radix(i, 3, eff_res);
  std::string j_str = int_to_radix(j, 3, eff_res);
  
  // Use lookup table for each (i,j) digit pair
  for (int idx = 0; idx < eff_res; idx++) {
    int i_digit = i_str[idx] - '0';
    int j_digit = j_str[idx] - '0';
    result += encode_table[i_digit][j_digit];
  }
  
  // Trim last digit if Class II
  if (!is_class_i && result.length() > 0) {
    result.pop_back();
  }
  
  return result;
}

void decode(const std::string& z3_str, int resolution,
            long long& i, long long& j) {
  std::string adjusted = z3_str;
  
  // Pad with "0" if Class II (odd length) to make even
  // This matches DGGRID line 220-221
  if (adjusted.length() % 2 == 1) {
    adjusted += "0";
  }
  
  // Reverse lookup table from DGGRID (lines 230-257)
  std::string i_str, j_str;
  
  for (size_t idx = 0; idx < adjusted.length(); idx += 2) {
    std::string z3code = adjusted.substr(idx, 2);
    
    // Decode using DGGRID's exact reverse mapping
    if (z3code == "00") {
      i_str += "0";
      j_str += "0";
    } else if (z3code == "22") {
      i_str += "0";
      j_str += "1";
    } else if (z3code == "21") {
      i_str += "0";
      j_str += "2";
    } else if (z3code == "01") {
      i_str += "1";
      j_str += "0";
    } else if (z3code == "02") {
      i_str += "1";
      j_str += "1";
    } else if (z3code == "20") {
      i_str += "1";
      j_str += "2";
    } else if (z3code == "12") {
      i_str += "2";
      j_str += "0";
    } else if (z3code == "10") {
      i_str += "2";
      j_str += "1";
    } else if (z3code == "11") {
      i_str += "2";
      j_str += "2";
    } else {
      throw std::runtime_error("hex_index_z3: invalid Z3 code: " + z3code);
    }
  }
  
  i = radix_to_int(i_str, 3);
  j = radix_to_int(j_str, 3);
}

} // namespace z3
} // namespace hexify
