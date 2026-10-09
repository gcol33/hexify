// hex_index.cpp
// Unified indexing for ISEA aperture 3, 4, 7

#include "cell_index.h"
#include "index_z3.h"
#include "index_zorder.h"
#include "index_z7.h"
#include "hex9.h"
#include "hex9_solid.h"
#include "coordinate_transforms.h"
#include "polyhedron.h"
#include "constants.h"
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace hexify {

// ---------------------------------------------------------------------------
// Aperture-7 Z7 index path
//
// hexify stores each ap7 cell as its exact resolution-r surrogate (quad, i, j).
// The Z7 index (z7::encode/decode) expands it to its Class I substrate IJK
// (exact integer, one aperture-7 level at odd resolutions), encodes it, and
// coarsens it back on decode. On the icosahedron the string is IGEO7's Z7
// index as DGGRID writes it: the leading field is the base cell, which is not
// always the cell's quad.
// ---------------------------------------------------------------------------
namespace {
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
  if (aperture == 9 || index_type == IndexType::HEX9) {
    return aperture == 9 && (index_type == IndexType::HEX9 || index_type == IndexType::AUTO);
  }
  if (index_type == IndexType::AUTO) return true;
  if (index_type == IndexType::ZORDER) return true;
  if (index_type == IndexType::Z3 && aperture == 3) return true;
  if (index_type == IndexType::Z7 && aperture == 7) return true;
  return false;
}

IndexType get_default_index_type(int aperture) {
  if (aperture == 9) return IndexType::HEX9;
  if (aperture == 3) return IndexType::Z3;
  if (aperture == 7) return IndexType::Z7;
  if (aperture == 4) return IndexType::ZORDER;
  return IndexType::ZORDER;
}

namespace {

// Aperture 9 writes libhex9's labels only
void check_hex9_type(IndexType t) {
  if (t != IndexType::AUTO && t != IndexType::HEX9) {
    throw std::runtime_error("hex_index: aperture 9 (Hex9) takes index_type hex9");
  }
}

// A Hex9 label's cell, as its quad and stored (i, j), and its level
void hex9_label_cell(const std::string& index, int& quad, long long& i, long long& j,
                     int& resolution) {
  hex9::OctPoint c;
  if (!hex9::parse_label(index, resolution, c)) {
    throw std::runtime_error("hex_index: \"" + index + "\" is no Hex9 label");
  }
  hex9::cell_quad_ij(c, resolution, quad, i, j);
}

// The centre of the Hex9 cell a label names, and its level
hex9::OctPoint hex9_label_centre(const std::string& index, int& resolution) {
  hex9::OctPoint c;
  if (!hex9::parse_label(index, resolution, c)) {
    throw std::runtime_error("hex_index: \"" + index + "\" is no Hex9 label");
  }
  return c;
}

} // anon

std::string cell_to_index(int face, long long i, long long j,
                          int resolution, int aperture,
                          IndexType index_type) {
  if (aperture == 9) {
    check_hex9_type(index_type);
    hex9::OctPoint c;
    if (!hex9::quad_ij_lattice(face, i, j, resolution, c) || !hex9::is_centre(c)) {
      throw std::runtime_error("hex_index: a Hex9 quad coordinate that is no cell");
    }
    return hex9::label(c, resolution);
  }
  // Face validation depends on aperture and resolution
  const SolidTopology& t = topo();
  const int max_face = (aperture == 3 && resolution > 0) ? t.n_faces - 1 : t.n_quads() - 1;
  if (face < 0 || face > max_face) {
    throw std::runtime_error("hex_index: invalid face number (must be 0-" +
                             std::to_string(max_face) + ")");
  }
  
  if (resolution < 0) {
    throw std::runtime_error("hex_index: invalid resolution");
  }
  if (resolution > kMaxResolution) {
    throw std::runtime_error("hex_index: resolution exceeds the maximum resolution");
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
    // then encode.
    long long sub_i, sub_j;
    ap7_surrogate_to_substrate_ijk(i, j, resolution, sub_i, sub_j);
    return z7::encode(face, sub_i, sub_j, resolution);
  }
  
  return result;
}

void index_to_cell(const std::string& index, int aperture,
                   IndexType index_type,
                   int& face, long long& i, long long& j, int& resolution) {
  if (aperture == 9) {
    check_hex9_type(index_type);
    hex9_label_cell(index, face, i, j, resolution);
    return;
  }
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
    // A Z7 index's leading field is not always a bare quad, so its decoder
    // reads it.
    if (index_type == IndexType::Z7) {
      z7::decode(index, face, i, j);
    }
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
    z7::decode(index, quadNum, sub_i, sub_j);
    // Coarsen the Class I substrate IJK back to the resolution-r surrogate
    // (identity for even res) that hexify stores.
    ap7_substrate_to_surrogate_ijk(sub_i, sub_j, resolution, i, j);
    face = quadNum;
  }
}

std::string get_parent_index(const std::string& index, int aperture,
                             IndexType index_type) {
  if (aperture == 9) {
    check_hex9_type(index_type);
    int level;
    const hex9::OctPoint c = hex9_label_centre(index, level);
    if (level == 0) throw std::runtime_error("hex_index: cannot get parent of resolution 0");
    return hex9::label(hex9::ancestor(c, level, level - 1), level - 1);
  }
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

  std::vector<std::string> children;

  if (aperture == 9) {
    check_hex9_type(index_type);
    int level;
    const hex9::OctPoint c = hex9_label_centre(index, level);
    if (level >= hex9::kMaxLevel) return children;
    hex9::OctPoint kids[9];
    hex9::children(c, level, kids);
    for (const hex9::OctPoint& k : kids) children.push_back(hex9::label(k, level + 1));
    return children;
  }

  if (aperture != 3 && aperture != 4 && aperture != 7) {
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
  if (child_res > kMaxResolution) {
    return children;
  }

  // Every index but aperture 7's z-order spells one refinement level as one
  // digit naming which child was taken, so the children are the parent index
  // with each digit appended -- the exact strings get_parent_index() strips
  // back to this one, which is what makes the two operations inverse.
  // A Z7 string whose first nonzero digit is the direction a pentagon lacks
  // names no cell, so a pentagon's descendants are the six that read back.
  if (!(index_type == IndexType::ZORDER && aperture == 7)) {
    for (int digit = 0; digit < aperture; digit++) {
      std::string child = index + std::to_string(digit);
      if (index_type == IndexType::Z7 && z7::in_deleted_subsequence(child)) continue;
      children.push_back(child);
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
  if (aperture == 9) return static_cast<int>(index.substr(0, index.find('.')).size()) - 1;
  if (index.length() <= 2) return 0;

  if (index_type == IndexType::AUTO) {
    index_type = get_default_index_type(aperture);
  }

  if (!is_valid_index_type(aperture, index_type)) return 0;

  return static_cast<int>((index.length() - 2) /
                          index_digits_per_level(aperture, index_type));
}

} // namespace hexify
