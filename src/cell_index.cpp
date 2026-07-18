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
// Aperture-7 surrogate <-> DGGRID Q2DI reconciliation (Z7 index path only)
//
// z7::encode/decode are a faithful port of DGGRID's DgZ7StringRF, which
// quantifies a DgQ2DICoord expressed on the quad's *Class I substrate* at
// numClassI = (res + 1) / 2 aperture-7 levels (scale 7^numClassI). hexify's
// canonical ap7 cell coordinate is the "surrogate": Class I for even res,
// Class II for odd res, at scale sqrt(7)^res. For even res the two frames
// coincide (7^(res/2) == sqrt(7)^res, Class I); for odd res they differ by a
// single aperture-7 refinement (the surrogate is one Class II level coarser
// than the Class I substrate the encoder walks).
//
// We reconcile geometrically through the shared, already-validated quad_xy
// frame rather than duplicating any per-quad lattice algebra here, so this
// adapter never touches the forward/inverse quantization pipeline.
// ---------------------------------------------------------------------------
namespace {

inline int ap7_classI_res(int res) { return 2 * ((res + 1) / 2); }  // res if even, res+1 if odd

void ap7_surrogate_to_q2di(long long sur_i, long long sur_j, int res,
                           long long& q_i, long long& q_j) {
  int cres = ap7_classI_res(res);
  if (cres == res) { q_i = sur_i; q_j = sur_j; return; }
  double qx, qy;
  surrogate_ij_to_quad_xy_ap7(sur_i, sur_j, res, qx, qy);
  quad_xy_to_surrogate_ij_ap7(qx, qy, cres, q_i, q_j);
}

void ap7_q2di_to_surrogate(long long q_i, long long q_j, int res,
                           long long& sur_i, long long& sur_j) {
  int cres = ap7_classI_res(res);
  if (cres == res) { sur_i = q_i; sur_j = q_j; return; }
  double qx, qy;
  surrogate_ij_to_quad_xy_ap7(q_i, q_j, cres, qx, qy);
  quad_xy_to_surrogate_ij_ap7(qx, qy, res, sur_i, sur_j);
}

} // anonymous namespace

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
    // Z7 encode already includes the base cell in its output. hexify stores ap7
    // cells as surrogates; the DGGRID Z7 encoder expects a Q2DI on the Class I
    // substrate, so reconcile first (identity for even res).
    long long q_i, q_j;
    ap7_surrogate_to_q2di(i, j, resolution, q_i, q_j);
    return z7::encode(face, q_i, q_j, resolution);
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
    long long q_i, q_j;
    z7::decode(index, resolution, quadNum, q_i, q_j);
    // Decode returns a Class I substrate Q2DI; convert back to a surrogate
    // (identity for even res) so downstream cell conversion sees hexify coords.
    ap7_q2di_to_surrogate(q_i, q_j, resolution, i, j);
    face = quadNum;
  }
}

uint64_t index_to_uint64(const std::string& index, int aperture,
                         IndexType index_type) {
  // Pack index string into a 64-bit integer:
  // Bits 60-63: face/quad (4 bits, 0-11)
  // Remaining bits: resolution digits packed per aperture
  if (index.length() < 2) {
    throw std::runtime_error("hex_index: index too short for uint64 conversion");
  }

  int face = parse_quad(index);
  uint64_t result = static_cast<uint64_t>(face) << 60;

  if (index.length() == 2) return result;

  std::string digits = index.substr(2);
  int bits_per_digit;
  if (aperture == 3) bits_per_digit = 2;       // 0-2 needs 2 bits
  else if (aperture == 4) bits_per_digit = 2;   // 0-3 needs 2 bits
  else if (aperture == 7) bits_per_digit = 3;   // 0-6 needs 3 bits
  else throw std::runtime_error("hex_index: unsupported aperture for uint64");

  // Check we don't overflow (60 bits available for digits)
  if (static_cast<int>(digits.length()) * bits_per_digit > 60) {
    throw std::runtime_error("hex_index: index too long for uint64 (overflow)");
  }

  for (size_t r = 0; r < digits.length(); r++) {
    int digit = digits[r] - '0';
    int shift = 60 - static_cast<int>((r + 1) * bits_per_digit);
    result |= (static_cast<uint64_t>(digit) << shift);
  }

  return result;
}

std::string uint64_to_index(uint64_t value, int resolution, int aperture,
                            IndexType index_type) {
  // Unpack uint64 back to index string
  int face = static_cast<int>((value >> 60) & 0xF);
  std::string result = format_quad(face);

  if (resolution == 0) return result;

  int bits_per_digit;
  int max_digit;
  if (aperture == 3) { bits_per_digit = 2; max_digit = 2; }
  else if (aperture == 4) { bits_per_digit = 2; max_digit = 3; }
  else if (aperture == 7) { bits_per_digit = 3; max_digit = 6; }
  else throw std::runtime_error("hex_index: unsupported aperture for uint64");

  uint64_t mask = (1ULL << bits_per_digit) - 1;

  for (int r = 0; r < resolution; r++) {
    int shift = 60 - (r + 1) * bits_per_digit;
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
  
  if ((index_type == IndexType::Z3 || 
       (index_type == IndexType::ZORDER && aperture == 3))) {
    std::string idx_str = index.substr(2);
    if (idx_str.length() % 2 == 0) {
      return index.substr(0, index.length() - 2);
    } else {
      return index.substr(0, index.length() - 1);
    }
  }
  
  if ((index_type == IndexType::ZORDER || index_type == IndexType::Z7) 
      && aperture == 7) {
    return index.substr(0, index.length() - 1);
  }
  
  return index.substr(0, index.length() - 1);
}

std::vector<std::string> get_children_indices(const std::string& index, 
                                              int aperture,
                                              IndexType index_type) {
  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }
  
  std::vector<std::string> children;
  
  int num_children = (aperture == 3) ? 3 : 
                     (aperture == 4) ? 4 : 
                     (aperture == 7) ? 7 : 0;
  
  if (num_children == 0) {
    throw std::runtime_error("hex_index: invalid aperture");
  }
  
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

  if (aperture == 7 && index_type == IndexType::Z7) {
    // For Z7, children are simply parent_index + digit (0-6)
    // The Z7 encoding is hierarchical - no coordinate computation needed
    for (int digit = 0; digit < 7; digit++) {
      std::string child = index + std::to_string(digit);
      children.push_back(child);
    }
    return children;
  }
  
  // For other aperture 7 cases or coordinate-based approach:
  if (aperture == 7) {
    long long base_i = parent_i * 7;
    long long base_j = parent_j * 7;
    
    std::vector<std::pair<long long, long long>> hex_offsets = {
      {0, 0},
      {1, 0},
      {0, 1},
      {-1, 1},
      {-1, 0},
      {0, -1},
      {1, -1}
    };
    
    for (const auto& offset : hex_offsets) {
      long long child_i = base_i + offset.first;
      long long child_j = base_j + offset.second;

      std::string child = cell_to_index(face, child_i, child_j,
                                        child_res, aperture, index_type);
      children.push_back(child);
    }
    
    return children;
  }
  
  int radix = (aperture == 4) ? 2 : aperture;
  
  if (aperture == 3) {
    radix = 3;
  }
  
  long long i_min = parent_i * radix;
  long long i_max = (parent_i + 1) * radix - 1;
  long long j_min = parent_j * radix;
  long long j_max = (parent_j + 1) * radix - 1;
  
  for (long long child_i = i_min; child_i <= i_max; child_i++) {
    for (long long child_j = j_min; child_j <= j_max; child_j++) {
      
      if (aperture == 3 && child_res % 2 == 1) {
        std::vector<int> j_required = {0, 2, 1};
        int i_mod = child_i % 3;
        int j_mod = child_j % 3;
        int req_j_mod = j_required[i_mod];
        
        if (j_mod != req_j_mod) {
          continue;
        }
      }
      
      std::string child_index = cell_to_index(face, child_i, child_j,
                                               child_res, aperture, index_type);
      children.push_back(child_index);

      if (children.size() >= static_cast<size_t>(num_children)) {
        return children;
      }
    }
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
  
  std::string idx_str = index.substr(2);
  
  if ((index_type == IndexType::Z3 || 
       (index_type == IndexType::ZORDER && aperture == 3))) {
    return idx_str.length();
  } else if (index_type == IndexType::ZORDER && aperture == 4) {
    return idx_str.length();
  } else if (index_type == IndexType::ZORDER && aperture == 7) {
    // Aperture 7 zorder uses 2 digits per resolution level
    return idx_str.length() / 2;
  } else if (index_type == IndexType::Z7) {
    // Z7 uses 1 digit per resolution level
    return idx_str.length();
  }
  
  return 0;
}

} // namespace hexify
