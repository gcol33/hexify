#pragma once
// Hex9 (Griffin 2026): the shifted-aperture-9 hexagonal grid on the octahedron.
//
// Lattice. Level L divides every edge of the octahedron into n = 3^(L+1)
// steps. A point of that lattice is written in the solid's standard frame --
// X at vertex 1, Y at vertex 2, Z at vertex 0 -- as integers (x, y, z) with
// |x| + |y| + |z| = n; on its face its weights on the face's three vertices
// are |x|, |y| and |z|. A step between neighbouring lattice points changes
// |x| - |y| by 1 or 2 modulo 3, so (|x| - |y|) mod 3 colours the lattice in
// three classes, and the colour is the same read from either face of an edge.
// The cells of level L are the hexagons centred on the points of colour 1:
// six triangles of the lattice each, three on either side of the line
// through the centre where one weight is a multiple of 3. The vertices of the
// solid have colour 0, so no cell is centred on one; two cells meet at each
// vertex, one corner each, and share both edges there. 12 * 9^L cells in all.
//
// t_cells and d_cells. The octant faces are the level-0 triangles (t_cells);
// each splits into nine of the next level. A t_cell is mode 0 when it points
// the way a face with an even number of negative axes does, mode 1 otherwise
// (Griffin's two-colouring of faces). The line through a cell's centre cuts
// the cell into two half-hexagons (d_cells), each inside one t_cell of the
// cell's own level, one in a mode-0 t_cell and one in a mode-1 t_cell.
//
// Address. A cell is named by the chain of t_cells holding a point of its
// mode-0 half, one digit per level, as libhex9 names it (Griffin 2026,
// Sec. 10b): digit 0 picks one of the twelve level-0 cells, digits 1..L
// one of nine children. The cell ID is that digit string read as a number,
// digit 0 counting 9^L, plus one, so the IDs of a level are 1 .. 12 * 9^L.
// The canonical parent of a cell is the cell of the level above holding its
// mode-0 half; every cell has nine children under that rule. Its canonical
// ancestor k levels up is the cell of that level holding its mode-0 half, as
// libhex9's hex9_cell_ancestor finds it; for k > 1 that is not always the
// parent's parent, which follows the parent's mode-0 half instead.

#include <cstdint>
#include <vector>

namespace hexify {
namespace hex9 {

// The finest level whose 12 * 9^L cell IDs fit a signed 64-bit integer
constexpr int kMaxLevel = 18;

// A lattice point of level L, (x, y, z) in the solid's standard frame
struct OctPoint {
  int64_t c[3];
};

// n = 3^(L+1), the steps along an edge at level L
int64_t level_steps(int level);

// 12 * 9^level, the cells of a level
int64_t level_cells(int level);

// Whether the lattice point is a cell centre (colour 1)
bool is_centre(const OctPoint& p);

// The cell of `level` holding the point of the octahedron with weights
// (wx, wy, wz) >= 0 on the axes of the octant with signs (sx, sy, sz),
// wx + wy + wz = 1. Its centre is the colour-1 corner of the lattice
// triangle the point falls in.
OctPoint locate(const double w[3], const int s[3], int level);

// The ID of the cell centred on `centre`; 0 when `centre` is no cell centre.
int64_t encode(const OctPoint& centre, int level);

// The centre of the cell with ID `id`; false when `id` names no cell of
// `level`.
bool decode(int64_t id, int level, OctPoint& centre);

// The digit string of a cell ID: digits[0] in 0..11, digits[1..level] in
// 0..8. False when `id` is out of range.
bool id_digits(int64_t id, int level, int digits[]);

// The key tail libhex9 writes after a cell's digits, (c2 << 1) | r_mo, 0..5
int key_tail(const OctPoint& centre, int level);

// Whether no two cells of any one level share a digit string, which makes
// the cell ID of every level a bijection onto 1 .. 12 * 9^L. Decided on the
// digit tables by a search over pairs of t_cell chains reading the same
// digits: true when no two chains naming different cells can end together.
bool digit_strings_unique();

// The address tables, flattened in C order, for a reader of addresses
// outside C++ (the globe widget): child codes [2][9], child modes [2][9],
// child offsets [2][9][2], region codes [12], root c2 [2][9], root digits
// [8][3] and region digits [2][12][12][2] (255 for no child).
struct Tables {
  std::vector<int> child_code, child_mode, child_offset, rid_code, root_c2,
      root_digit, reg_hex;
};
Tables tables();

// The centre of the canonical ancestor of the cell at `centre`, of level
// `to_level` < `level`: the canonical parent when to_level = level - 1.
OctPoint ancestor(const OctPoint& centre, int level, int to_level);

// The centres of the nine canonical children, one level down
void children(const OctPoint& centre, int level, OctPoint out[9]);

} // namespace hex9
} // namespace hexify
