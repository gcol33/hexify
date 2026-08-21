// hex_index.cpp
// Unified indexing for ISEA aperture 3, 4, 7

#include "cell_index.h"
#include "index_z3.h"
#include "index_zorder.h"
#include "index_z7.h"
#include "coordinate_transforms.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace hexify {

// ---------------------------------------------------------------------------
// Aperture-7 Z7 index path
//
// hexify stores each ap7 cell as its exact resolution-r surrogate (quad, i, j)
// -- a bijective, geographically consistent coordinate. The Z7 index is the
// bijective hierarchical encoding of it (z7::encode_bijective/decode_bijective):
// the surrogate is expanded to its Class I substrate IJK (exact integer, one
// aperture-7 level at odd resolutions), encoded, and coarsened back on decode.
// The digits follow DGGRID's DgZ7StringRF except in the pentagon regions where
// DGGRID's own encoder is non-injective (it collides distinct cells); there
// hexify keeps a distinct, round-tripping index instead of the colliding one.
// The leading field is quad + 12 * seed rather than the bare quad, seed naming
// the level-0 point the hierarchy walk arrives at, which is the origin only for
// a cell whose whole ancestry lies inside its quad.
// ---------------------------------------------------------------------------
namespace {
  const int MAX_RES_AP3 = 30;
  const int MAX_RES_AP4 = 30;
  const int MAX_RES_AP7 = 20;
  
  std::string format_quad(int quadNum) {
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << quadNum;
    return oss.str();
  }
  
  int parse_quad(const std::string& index) {
    if (index.length() < 2) {
      throw std::runtime_error("hex_index: invalid index string (too short)");
    }
    std::string qstr = index.substr(0, 2);
    return std::stoi(qstr);  // stoi handles leading zeros correctly
  }
}

// Characters one resolution level spends in an index string, after the two
// leading quad digits. Aperture 7's z-order spells a level as a pair; Z3, Z7
// and the other z-order apertures spell it as a single digit. A cell's
// resolution and its parent's index are both read off this one number, so they
// cannot disagree about how deep a string is.
static size_t index_digits_per_level(int aperture, IndexType index_type) {
  if (index_type == IndexType::ZORDER && aperture == 7) return 2;
  return 1;
}

bool is_valid_index_type(int aperture, IndexType index_type) {
  if (index_type == IndexType::AUTO) return true;
  if (index_type == IndexType::ZORDER) return true;
  if (index_type == IndexType::Z3 && aperture == 3) return true;
  if (index_type == IndexType::Z7 && aperture == 7) return true;
  return false;
}

IndexType get_default_index_type(int aperture) {
  if (aperture == 3) return IndexType::Z3;
  if (aperture == 7) return IndexType::Z7;
  if (aperture == 4) return IndexType::ZORDER;
  return IndexType::ZORDER;
}

std::string cell_to_index(int face, long long i, long long j, 
                          int resolution, int aperture,
                          IndexType index_type) {
  // Face validation depends on aperture and resolution
  if (aperture == 3) {
    if (resolution == 0) {
      if (face < 0 || face > 11) {
        throw std::runtime_error("hex_index: invalid face number for aperture 3 resolution 0 (must be 0-11)");
      }
    } else {
      if (face < 0 || face > 19) {
        throw std::runtime_error("hex_index: invalid face number for aperture 3 (must be 0-19)");
      }
    }
  } else {
    // Aperture 4 and 7 only use faces 0-11
    if (face < 0 || face > 11) {
      throw std::runtime_error("hex_index: invalid face number (must be 0-11)");
    }
  }
  
  if (resolution < 0) {
    throw std::runtime_error("hex_index: invalid resolution");
  }
  
  if (aperture == 7 && resolution > MAX_RES_AP7) {
    throw std::runtime_error("hex_index: resolution exceeds max for aperture 7");
  }
  if (aperture == 3 && resolution > MAX_RES_AP3) {
    throw std::runtime_error("hex_index: resolution exceeds max for aperture 3");
  }
  if (aperture == 4 && resolution > MAX_RES_AP4) {
    throw std::runtime_error("hex_index: resolution exceeds max for aperture 4");
  }

  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }
  
  if (!is_valid_index_type(aperture, index_type)) {
    throw std::runtime_error("hex_index: invalid index_type for aperture");
  }
  
  std::string result = format_quad(face);
  
  if (resolution == 0) {
    return result;
  }
  
  if (index_type == IndexType::ZORDER) {
    if (aperture == 3) {
      result += zorder::encode_ap3(i, j, resolution);
    } else if (aperture == 4) {
      result += zorder::encode_ap4(i, j, resolution);
    } else if (aperture == 7) {
      result += zorder::encode_ap7(i, j, resolution);
    }
  } else if (index_type == IndexType::Z3) {
    result += z3::encode(i, j, resolution);
  } else if (index_type == IndexType::Z7) {
    // Z7 encode includes the base cell in its output. Expand the stored
    // surrogate to its Class I substrate IJK (exact, identity for even res),
    // then encode bijectively (quad fixed).
    long long sub_i, sub_j;
    ap7_surrogate_to_substrate_ijk(i, j, resolution, sub_i, sub_j);
    return z7::encode_bijective(face, sub_i, sub_j, resolution);
  }
  
  return result;
}

void index_to_cell(const std::string& index, int aperture,
                   IndexType index_type,
                   int& face, long long& i, long long& j, int& resolution) {
  if (index.length() < 2) {
    throw std::runtime_error("hex_index: invalid index string");
  }
  
  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }
  
  if (!is_valid_index_type(aperture, index_type)) {
    throw std::runtime_error("hex_index: invalid index_type for aperture");
  }
  
  face = parse_quad(index);
  
  if (index.length() == 2) {
    resolution = 0;
    i = 0;
    j = 0;
    return;
  }
  
  std::string index_str = index.substr(2);
  
  if (index_type == IndexType::ZORDER) {
    if (aperture == 3) {
      resolution = index_str.length();
      zorder::decode_ap3(index_str, resolution, i, j);
    } else if (aperture == 4) {
      resolution = index_str.length();
      zorder::decode_ap4(index_str, i, j);
    } else if (aperture == 7) {
      resolution = index_str.length() / 2;
      zorder::decode_ap7(index_str, resolution, i, j);
    }
  } else if (index_type == IndexType::Z3) {
    resolution = index_str.length();
    z3::decode(index_str, resolution, i, j);
  } else if (index_type == IndexType::Z7) {
    // Z7 decode expects the full index (including base cell)
    // Calculate resolution from index length first
    resolution = index.length() - 2;
    int quadNum = face;
    long long sub_i, sub_j;
    z7::decode_bijective(index, resolution, quadNum, sub_i, sub_j);
    // Coarsen the Class I substrate IJK back to the resolution-r surrogate
    // (identity for even res) that hexify stores.
    ap7_substrate_to_surrogate_ijk(sub_i, sub_j, resolution, i, j);
    face = quadNum;
  }
}

uint64_t index_to_uint64(const std::string& index, int aperture,
                         IndexType index_type) {
  // Pack index string into a 64-bit integer:
  // Bits 60-63: face/quad (4 bits, 0-11)
  // Bits 57-59: aperture-7 Z7 hierarchy seed (3 bits), 0 for every other index
  // Remaining bits: resolution digits packed per aperture
  if (index.length() < 2) {
    throw std::runtime_error("hex_index: index too short for uint64 conversion");
  }

  int lead = parse_quad(index);
  uint64_t result = (static_cast<uint64_t>(lead % 12) << 60) |
                    (static_cast<uint64_t>(lead / 12) << 57);

  if (index.length() == 2) return result;

  std::string digits = index.substr(2);
  int bits_per_digit;
  if (aperture == 3) bits_per_digit = 2;       // 0-2 needs 2 bits
  else if (aperture == 4) bits_per_digit = 2;   // 0-3 needs 2 bits
  else if (aperture == 7) bits_per_digit = 3;   // 0-6 needs 3 bits
  else throw std::runtime_error("hex_index: unsupported aperture for uint64");

  // Aperture 7 spends three bits on the hierarchy seed; the others start at 60.
  const int digit_base = (aperture == 7) ? 57 : 60;

  if (static_cast<int>(digits.length()) * bits_per_digit > digit_base) {
    throw std::runtime_error("hex_index: index too long for uint64 (overflow)");
  }

  for (size_t r = 0; r < digits.length(); r++) {
    int digit = digits[r] - '0';
    int shift = digit_base - static_cast<int>((r + 1) * bits_per_digit);
    result |= (static_cast<uint64_t>(digit) << shift);
  }

  return result;
}

std::string uint64_to_index(uint64_t value, int resolution, int aperture,
                            IndexType index_type) {
  // Unpack uint64 back to index string
  int face = static_cast<int>((value >> 60) & 0xF);
  int seed = static_cast<int>((value >> 57) & 0x7);
  std::string result = format_quad(face + 12 * seed);

  if (resolution == 0) return result;

  int bits_per_digit;
  int max_digit;
  if (aperture == 3) { bits_per_digit = 2; max_digit = 2; }
  else if (aperture == 4) { bits_per_digit = 2; max_digit = 3; }
  else if (aperture == 7) { bits_per_digit = 3; max_digit = 6; }
  else throw std::runtime_error("hex_index: unsupported aperture for uint64");

  const int digit_base = (aperture == 7) ? 57 : 60;
  uint64_t mask = (1ULL << bits_per_digit) - 1;

  for (int r = 0; r < resolution; r++) {
    int shift = digit_base - (r + 1) * bits_per_digit;
    int digit = static_cast<int>((value >> shift) & mask);
    if (digit > max_digit) {
      throw std::runtime_error("hex_index: uint64 value has an out-of-range digit for this aperture");
    }
    result += static_cast<char>('0' + digit);
  }

  return result;
}

std::string get_parent_index(const std::string& index, int aperture,
                             IndexType index_type) {
  if (index.length() <= 2) {
    throw std::runtime_error("hex_index: cannot get parent of resolution 0");
  }

  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }

  size_t digits = index_digits_per_level(aperture, index_type);
  if (index.length() - 2 < digits) {
    throw std::runtime_error("hex_index: index is too short for its index type");
  }

  return index.substr(0, index.length() - digits);
}

std::vector<std::string> get_children_indices(const std::string& index,
                                              int aperture,
                                              IndexType index_type) {
  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }

  if (aperture != 3 && aperture != 4 && aperture != 7) {
    throw std::runtime_error("hex_index: invalid aperture");
  }

  std::vector<std::string> children;

  int face, parent_res;
  long long parent_i, parent_j;
  index_to_cell(index, aperture, index_type, face, parent_i, parent_j, parent_res);

  int child_res = parent_res + 1;

  // Children beyond the aperture's max resolution are desired to be omitted
  // (a query one level below the finest supported resolution just yields no
  // children), so this is checked explicitly up front rather than relying
  // on a catch-all around cell_to_index() below to swallow the resulting
  // "resolution exceeds max" error -- which would also hide unrelated bugs.
  int max_res_for_aperture = (aperture == 7) ? MAX_RES_AP7 :
                             (aperture == 3) ? MAX_RES_AP3 :
                             (aperture == 4) ? MAX_RES_AP4 : -1;
  if (max_res_for_aperture >= 0 && child_res > max_res_for_aperture) {
    return children;
  }

  // Every index but aperture 7's z-order spells one refinement level as one
  // digit naming which child was taken, so the children are the parent index
  // with each digit appended -- the exact strings get_parent_index() strips
  // back to this one, which is what makes the two operations inverse.
  if (!(index_type == IndexType::ZORDER && aperture == 7)) {
    for (int digit = 0; digit < aperture; digit++) {
      children.push_back(index + std::to_string(digit));
    }
    return children;
  }

  // Aperture 7's z-order spells a level as a radix-7 digit of i beside one of
  // j, which scales both coordinates by 7 rather than naming a child. Its
  // children are the seven aperture-7 cells around the scaled parent centre.
  static const long long hex_offsets[7][2] = {
    {0, 0}, {1, 0}, {0, 1}, {-1, 1}, {-1, 0}, {0, -1}, {1, -1}
  };

  for (const auto& offset : hex_offsets) {
    children.push_back(cell_to_index(face, parent_i * 7 + offset[0],
                                     parent_j * 7 + offset[1],
                                     child_res, aperture, index_type));
  }

  return children;
}

int compare_indices(const std::string& idx1, const std::string& idx2) {
  return idx1.compare(idx2);
}

int get_index_resolution(const std::string& index, int aperture,
                         IndexType index_type) {
  if (index.length() <= 2) return 0;

  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }

  if (!is_valid_index_type(aperture, index_type)) return 0;

  return static_cast<int>((index.length() - 2) /
                          index_digits_per_level(aperture, index_type));
}

} // namespace hexify
