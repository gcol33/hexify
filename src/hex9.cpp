// hex9.cpp - Hex9 lattice, addresses and hierarchy (see hex9.h)
//
// The digit tables below are the Hex9 address definition of libhex9
// (https://github.com/MrBenGriffin/libhex9, core/h9_addressing.h,
// core/h9_math.h and core/h9_uv_lattice.h), Copyright (c) 2025-2026 Ben Griffin,
// Apache License 2.0; see inst/COPYRIGHTS. The lattice, the canonical half
// and the hierarchy are computed here in exact integers.

#include "hex9.h"
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace hexify {
namespace hex9 {

namespace {

// ---- libhex9's chart of an octant ----
//
// libhex9 lays each octant out in a plane frame ("UV" lattice units) where a
// mode-0 octant points down, with its X, Z and Y corners at (3, 1), (0, -2)
// and (-3, 1), and a mode-1 octant points up, X, Y and Z at (3, -1), (-3, -1)
// and (0, 2). A t_cell of the chart points down exactly when it has mode 0.
// At layer M the chart is scaled by 3^M, and a point with weights
// (a_x, a_y, a_z) summing to 3^(M+1) sits at u = a_x - a_y,
// v = +-(3^M - a_z) (+ for a mode-0 octant).

// The nine children of a t_cell, in libhex9's order: their chart cell codes,
// modes, and centroid offsets in child-layer units, for a mode-0 (down) and a
// mode-1 (up) parent.
const uint8_t kChildCode[2][9] = {
  {0x21, 0x25, 0x26, 0x2a, 0x2b, 0x35, 0x39, 0x3a, 0x49},
  {0x16, 0x25, 0x26, 0x2a, 0x34, 0x35, 0x39, 0x3a, 0x3e}
};
const int kChildMode[2][9] = {
  {0, 1, 0, 1, 0, 0, 1, 0, 0},
  {1, 1, 0, 1, 1, 0, 1, 0, 1}
};
const int kChildOffset[2][9][2] = {
  {{-6, 2}, {-3, 1}, {0, 2}, {3, 1}, {6, 2}, {-3, -1}, {0, -2}, {3, -1}, {0, -4}},
  {{0, 4}, {-3, 1}, {0, 2}, {3, 1}, {-6, -2}, {-3, -1}, {0, -2}, {3, -1}, {6, -2}}
};

// The region id of each chart cell code, rid 0..11, and its inverse; a
// region's mode is rid & 1.
const uint8_t kRidCode[12] = {
  0x49, 0x16, 0x2b, 0x34, 0x21, 0x3e, 0x26, 0x39, 0x35, 0x2a, 0x3a, 0x25
};

// The c2 edge label of the level-0 half-hexagon holding each child of a
// mode-0 and a mode-1 octant, in kChildCode order.
const uint8_t kRootC2[2][9] = {
  {2, 2, 0, 0, 0, 2, 1, 1, 1},
  {2, 1, 2, 2, 1, 1, 0, 0, 0}
};

// The level-0 cell (digit 0..11) of octant oid's half-hexagon with label c2
const uint8_t kRootDigit[8][3] = {
  {0, 4, 5}, {1, 5, 7}, {2, 6, 4}, {3, 7, 6},
  {0, 8, 10}, {1, 9, 8}, {2, 10, 11}, {3, 11, 9}
};

// kRegHex[gp_mo][parent rid][child rid] = {digit, c2}: the digit naming the
// cell whose half-hexagon in the parent t_cell holds the child t_cell, and
// that half-hexagon's c2 label, given the grandparent's mode. 255 where the
// child is no child of the parent.
const uint8_t kRegHex[2][12][12][2] = {
  {
    {{2,1},{255,255},{4,0},{255,255},{8,2},{255,255},{4,0},{2,1},{8,2},{4,0},{2,1},{8,2}},
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{8,1},{255,255},{2,0},{255,255},{4,2},{255,255},{2,0},{8,1},{4,2},{2,0},{8,1},{4,2}},
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{4,1},{255,255},{8,0},{255,255},{2,2},{255,255},{8,0},{4,1},{2,2},{8,0},{4,1},{2,2}},
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{0,1},{255,255},{3,0},{255,255},{6,2},{255,255},{3,0},{0,1},{6,2},{3,0},{0,1},{6,2}},
    {{255,255},{0,2},{255,255},{6,1},{255,255},{4,0},{0,2},{4,0},{6,1},{0,2},{4,0},{6,1}},
    {{6,1},{255,255},{0,0},{255,255},{3,2},{255,255},{0,0},{6,1},{3,2},{0,0},{6,1},{3,2}},
    {{255,255},{4,2},{255,255},{0,1},{255,255},{6,0},{4,2},{6,0},{0,1},{4,2},{6,0},{0,1}},
    {{3,1},{255,255},{6,0},{255,255},{0,2},{255,255},{6,0},{3,1},{0,2},{6,0},{3,1},{0,2}},
    {{255,255},{6,2},{255,255},{4,1},{255,255},{0,0},{6,2},{0,0},{4,1},{6,2},{0,0},{4,1}}
  },
  {
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{255,255},{2,2},{255,255},{8,1},{255,255},{5,0},{2,2},{5,0},{8,1},{2,2},{5,0},{8,1}},
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{255,255},{5,2},{255,255},{2,1},{255,255},{8,0},{5,2},{8,0},{2,1},{5,2},{8,0},{2,1}},
    {{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255},{255,255}},
    {{255,255},{8,2},{255,255},{5,1},{255,255},{2,0},{8,2},{2,0},{5,1},{8,2},{2,0},{5,1}},
    {{1,1},{255,255},{5,0},{255,255},{7,2},{255,255},{5,0},{1,1},{7,2},{5,0},{1,1},{7,2}},
    {{255,255},{1,2},{255,255},{7,1},{255,255},{3,0},{1,2},{3,0},{7,1},{1,2},{3,0},{7,1}},
    {{7,1},{255,255},{1,0},{255,255},{5,2},{255,255},{1,0},{7,1},{5,2},{1,0},{7,1},{5,2}},
    {{255,255},{3,2},{255,255},{1,1},{255,255},{7,0},{3,2},{7,0},{1,1},{3,2},{7,0},{1,1}},
    {{5,1},{255,255},{7,0},{255,255},{1,2},{255,255},{7,0},{5,1},{1,2},{7,0},{5,1},{1,2}},
    {{255,255},{7,2},{255,255},{3,1},{255,255},{1,0},{7,2},{1,0},{3,1},{7,2},{1,0},{3,1}}
  }
};

// kHexReg[digit][c_mo][c2] = {rid, c_mo, c2}: the backward pass, from the
// state below a level to that level's region and the state above it.
const uint8_t kHexReg[9][2][3][3] = {
  {{{8,0,2},{6,0,0},{10,0,1}},{{11,0,2},{9,0,0},{7,0,1}}},
  {{{8,1,1},{6,1,2},{10,1,0}},{{11,1,1},{9,1,2},{7,1,0}}},
  {{{2,0,0},{0,0,1},{4,0,2}},{{5,1,0},{3,1,1},{1,1,2}}},
  {{{6,0,0},{10,0,1},{8,0,2}},{{7,1,0},{11,1,1},{9,1,2}}},
  {{{0,0,1},{4,0,2},{2,0,0}},{{7,0,1},{11,0,2},{9,0,0}}},
  {{{6,1,2},{10,1,0},{8,1,1}},{{1,1,2},{5,1,0},{3,1,1}}},
  {{{10,0,1},{8,0,2},{6,0,0}},{{9,0,0},{7,0,1},{11,0,2}}},
  {{{10,1,0},{8,1,1},{6,1,2}},{{9,1,2},{7,1,0},{11,1,1}}},
  {{{4,0,2},{2,0,0},{0,0,1}},{{3,1,1},{1,1,2},{5,1,0}}}
};

// Centroid offsets, in units of a t_cell's own layer, of the three cells
// whose half-hexagons lie in it, by the t_cell's mode and the c2 label of the
// edge each lies on (0 the horizontal edge, 1 the edge rising to the right,
// 2 the edge falling to the right), a third of the way along it.
const int kHalfCentre[2][3][2] = {
  {{1, 1}, {1, -1}, {-2, 0}},
  {{1, -1}, {-2, 0}, {1, 1}}
};

// ---- Lattice helpers ----

inline int64_t pow3(int k) {
  int64_t p = 1;
  for (int i = 0; i < k; i++) p *= 3;
  return p;
}

inline int64_t mod3(int64_t a) { return ((a % 3) + 3) % 3; }

inline int64_t iabs(int64_t a) { return a < 0 ? -a : a; }

// libhex9's octant id of a sign triple: bit k set for a negative axis k
inline int oid_of(const int s[3]) {
  return (s[0] < 0 ? 1 : 0) | (s[1] < 0 ? 2 : 0) | (s[2] < 0 ? 4 : 0);
}

inline int mode_of_oid(int oid) {
  return ((oid >> 0) ^ (oid >> 1) ^ (oid >> 2)) & 1;
}

// The cell-code table index of `code` among a parent's nine children, or -1
int child_index(int parent_mode, uint8_t code) {
  for (int j = 0; j < 9; j++) {
    if (kChildCode[parent_mode][j] == code) return j;
  }
  return -1;
}

int rid_of_code(uint8_t code) {
  for (int r = 0; r < 12; r++) {
    if (kRidCode[r] == code) return r;
  }
  return -1;
}

// libhex9's band classifier of a point (u, v), relative to a parent t_cell
// centroid at parent scale `scale`, into a chart cell code: three families
// of parallel lines (h: v, p: v - u, n: v + u) at multiples of the scale.
uint8_t classify(int64_t u, int64_t v, int64_t scale, int p_mo) {
  const int64_t s1 = scale, s2 = 2 * scale;
  const int64_t ymx = v - u, ypx = v + u;
  int h, p, n;
  if (p_mo == 0) {
    h = (v > s2) ? 0 : (v >= s1) ? 1 : (v >= 0) ? 2 : (v > -s1) ? 3 : (v > -s2) ? 4 : 5;
    p = (ymx >= s2) ? 0 : (ymx >= 0) ? 1 : (ymx >= -s2) ? 2 : 3;
    n = (ypx < -s2) ? 0 : (ypx <= 0) ? 1 : (ypx < s2) ? 2 : 3;
  } else {
    h = (v > s2) ? 0 : (v > s1) ? 1 : (v >= 0) ? 2 : (v > -s1) ? 3 : (v > -s2) ? 4 : 5;
    p = (ymx > s2) ? 0 : (ymx > 0) ? 1 : (ymx > -s2) ? 2 : 3;
    n = (ypx <= -s2) ? 0 : (ypx <= 0) ? 1 : (ypx <= s2) ? 2 : 3;
  }
  return static_cast<uint8_t>((h << 4) | (p << 2) | n);
}

// A point inside one octant: the octant's signs and the point's weights on
// its X, Y and Z axes at some scale.
struct OctantPoint {
  int s[3];
  int64_t a[3];
};

// The octants holding a lattice point (two on an edge), with its weights.
int octants_of(const OctPoint& p, OctantPoint out[2]) {
  int zero = -1;
  OctantPoint o;
  for (int k = 0; k < 3; k++) {
    o.s[k] = p.c[k] < 0 ? -1 : 1;
    o.a[k] = iabs(p.c[k]);
    if (p.c[k] == 0) zero = k;
  }
  out[0] = o;
  if (zero < 0) return 1;
  out[1] = o;
  out[1].s[zero] = -1;
  return 2;
}

// Chart coordinates at layer M of a point with weights `a` summing to
// 3^(M+1) in an octant of mode `mode`.
inline void chart_uv(const int64_t a[3], int M, int mode, int64_t& u, int64_t& v) {
  u = a[0] - a[1];
  v = (mode == 0) ? pow3(M) - a[2] : a[2] - pow3(M);
}

// The inverse of chart_uv(): weights summing to 3^(M+1)
inline void chart_weights(int64_t u, int64_t v, int M, int mode, int64_t a[3]) {
  const int64_t n = pow3(M + 1);
  a[2] = (mode == 0) ? pow3(M) - v : pow3(M) + v;
  a[0] = (n - a[2] + u) / 2;
  a[1] = (n - a[2] - u) / 2;
}

// The six directions from a lattice point to the centroids of the triangles
// around it, in weights at three times the lattice scale.
const int64_t kAround[6][3] = {
  {2, -1, -1}, {1, 1, -2}, {-1, 2, -1}, {-2, 1, 1}, {-1, -1, 2}, {1, -2, 1}
};

// Mode of the t_cell of a layer holding a point strictly inside it: weights
// `c` at scale `scale` times the layer's division 3^layer of the octant.
inline int tcell_mode(const int64_t c[3], int64_t cell, int64_t division, int oct_mode) {
  const int64_t f = c[0] / cell + c[1] / cell + c[2] / cell;
  const bool inverted = (f == division - 2);
  return oct_mode ^ (inverted ? 1 : 0);
}

// The chain of t_cells holding a point, from the octant down `M` layers:
// region ids rid[0..M] (rid[0] the octant's) and the centroid of each t_cell
// in its own layer's chart units, origin[k] for k = 0..M.
struct Chain {
  int oid;
  int mode;
  int M;
  int rid[kMaxLevel + 3];
  int64_t origin[kMaxLevel + 3][2];
};

// Descends the point at chart (u, v) of layer M in octant oid.
Chain descend(int oid, int64_t u, int64_t v, int M) {
  Chain ch;
  ch.oid = oid;
  ch.mode = mode_of_oid(oid);
  ch.M = M;
  ch.rid[0] = ch.mode;
  ch.origin[0][0] = 0;
  ch.origin[0][1] = 0;
  int64_t scale = pow3(M);
  int p_mo = ch.mode;
  int64_t ou = 0, ov = 0;
  for (int k = 0; k < M; k++) {
    const uint8_t code = classify(u - ou, v - ov, scale, p_mo);
    const int j = child_index(p_mo, code);
    if (j < 0) throw std::logic_error("hex9: a point left its t_cell during descent");
    const int64_t child = scale / 3;
    ou += kChildOffset[p_mo][j][0] * child;
    ov += kChildOffset[p_mo][j][1] * child;
    ch.rid[k + 1] = rid_of_code(code);
    p_mo = kChildMode[p_mo][j];
    scale = child;
    // The child's centroid in its own layer's units, k + 1
    const int64_t down = pow3(M - k - 1);
    ch.origin[k + 1][0] = ou / down;
    ch.origin[k + 1][1] = ov / down;
  }
  return ch;
}

// The digit and c2 label at level k of a chain, 0 <= k < M: of the cell
// whose half-hexagon in t_cell k holds t_cell k + 1.
inline void chain_digit(const Chain& ch, int k, int& digit, int& c2) {
  if (k == 0) {
    const int j = child_index(ch.mode, kRidCode[ch.rid[1]]);
    c2 = kRootC2[ch.mode][j];
    digit = kRootDigit[ch.oid][c2];
    return;
  }
  const uint8_t* e = kRegHex[ch.rid[k - 1] & 1][ch.rid[k]][ch.rid[k + 1]];
  if (e[0] == 255) throw std::logic_error("hex9: inconsistent t_cell chain");
  digit = e[0];
  c2 = e[1];
}

// The centre, at level k, of the cell holding the chain's point
OctPoint chain_centre(const Chain& ch, int k) {
  int digit, c2;
  chain_digit(ch, k, digit, c2);
  const int mode = ch.rid[k] & 1;
  const int64_t u = ch.origin[k][0] + kHalfCentre[mode][c2][0];
  const int64_t v = ch.origin[k][1] + kHalfCentre[mode][c2][1];
  int64_t a[3];
  chart_weights(u, v, k, ch.mode, a);
  const int s[3] = {(ch.oid & 1) ? -1 : 1, (ch.oid & 2) ? -1 : 1, (ch.oid & 4) ? -1 : 1};
  OctPoint p;
  for (int i = 0; i < 3; i++) p.c[i] = s[i] * a[i];
  return p;
}

// A point inside the mode-0 half of the cell at `centre`: the centroid of
// one of its triangles there, as an octant point at three times the level's
// lattice scale.
OctantPoint mode0_point(const OctPoint& centre, int level) {
  OctantPoint oct[2];
  const int n_oct = octants_of(centre, oct);
  const int64_t coarse = 9;               // a level t_cell, in these units
  const int64_t division = pow3(level);   // level t_cells per octant edge
  for (int o = 0; o < n_oct; o++) {
    const int oct_mode = mode_of_oid(oid_of(oct[o].s));
    for (int d = 0; d < 6; d++) {
      int64_t c[3];
      bool inside = true;
      for (int k = 0; k < 3; k++) {
        c[k] = 3 * oct[o].a[k] + kAround[d][k];
        if (c[k] <= 0) inside = false;
      }
      if (!inside) continue;
      if (tcell_mode(c, coarse, division, oct_mode) == 0) {
        OctantPoint q;
        for (int k = 0; k < 3; k++) {
          q.s[k] = oct[o].s[k];
          q.a[k] = c[k];
        }
        return q;
      }
    }
  }
  throw std::logic_error("hex9: a cell without a mode-0 half");
}

// The chain of the mode-0 half of the cell at `centre`, down to its
// triangles one level below (M = level + 1).
Chain cell_chain(const OctPoint& centre, int level) {
  const OctantPoint q = mode0_point(centre, level);
  const int oid = oid_of(q.s);
  int64_t u, v;
  chart_uv(q.a, level + 1, mode_of_oid(oid), u, v);
  return descend(oid, u, v, level + 1);
}

void check_level(int level) {
  if (level < 0 || level > kMaxLevel) {
    throw std::invalid_argument("hex9: level must be in 0..18");
  }
}

} // anon

int64_t level_steps(int level) {
  check_level(level);
  return pow3(level + 1);
}

int64_t level_cells(int level) {
  check_level(level);
  int64_t n = 12;
  for (int k = 0; k < level; k++) n *= 9;
  return n;
}

bool is_centre(const OctPoint& p) {
  return mod3(iabs(p.c[0]) - iabs(p.c[1])) == 1;
}

OctPoint locate(const double w[3], const int s[3], int level) {
  check_level(level);
  const int64_t n = pow3(level + 1);
  int64_t f[3];
  int64_t sum = 0;
  double frac[3];
  for (int k = 0; k < 3; k++) {
    const double x = w[k] * static_cast<double>(n);
    double fl = static_cast<double>(static_cast<int64_t>(x));
    if (fl > x) fl -= 1.0;
    f[k] = static_cast<int64_t>(fl);
    if (f[k] < 0) f[k] = 0;
    if (f[k] > n) f[k] = n;
    frac[k] = x - static_cast<double>(f[k]);
    sum += f[k];
  }
  // The triangle's corners: f plus a unit vector when its floors sum to
  // n - 1 (pointing as the face does), f plus two when they sum to n - 2.
  // Rounding leaves the floors summing to n (a lattice point) or below
  // n - 2 (a point on a lattice line); step them onto the nearest triangle.
  while (sum > n - 1) {
    int k = 0;
    for (int i = 1; i < 3; i++) if (frac[i] < frac[k] && f[i] > 0) k = i;
    if (f[k] == 0) for (int i = 0; i < 3; i++) if (f[i] > 0) { k = i; break; }
    f[k]--;
    frac[k] += 1.0;
    sum--;
  }
  while (sum < n - 2) {
    int k = 0;
    for (int i = 1; i < 3; i++) if (frac[i] > frac[k]) k = i;
    f[k]++;
    frac[k] -= 1.0;
    sum++;
  }
  const bool upright = (sum == n - 1);
  for (int v = 0; v < 3; v++) {
    int64_t a[3];
    for (int k = 0; k < 3; k++) {
      a[k] = upright ? f[k] + (k == v ? 1 : 0) : f[k] + (k == v ? 0 : 1);
    }
    OctPoint p;
    for (int k = 0; k < 3; k++) p.c[k] = s[k] * a[k];
    if (is_centre(p)) return p;
  }
  throw std::logic_error("hex9: a lattice triangle without a cell centre");
}

bool id_digits(int64_t id, int level, int digits[]) {
  check_level(level);
  if (id < 1 || id > level_cells(level)) return false;
  int64_t rest = id - 1;
  for (int k = level; k >= 1; k--) {
    digits[k] = static_cast<int>(rest % 9);
    rest /= 9;
  }
  digits[0] = static_cast<int>(rest);
  return digits[0] < 12;
}

int64_t encode(const OctPoint& centre, int level) {
  check_level(level);
  if (!is_centre(centre)) return 0;
  const Chain ch = cell_chain(centre, level);
  int64_t id = 0;
  for (int k = 0; k <= level; k++) {
    int digit, c2;
    chain_digit(ch, k, digit, c2);
    id = id * 9 + digit;
  }
  return id + 1;
}

int key_tail(const OctPoint& centre, int level) {
  check_level(level);
  const Chain ch = cell_chain(centre, level);
  int digit, c2;
  chain_digit(ch, level, digit, c2);
  return (c2 << 1) | ch.mode;
}

bool decode(int64_t id, int level, OctPoint& centre) {
  int digits[kMaxLevel + 1];
  if (!id_digits(id, level, digits)) return false;
  // The digits fix the cell together with the key tail (c2, r_mo), which
  // seeds the backward pass; a canonical cell is read through its mode-0
  // half, so the pass starts in mode 0. The tail that reads back to the
  // same ID is the cell's.
  for (int r_mo = 0; r_mo < 2; r_mo++) {
    int oid = -1, root_c2 = -1;
    for (int o = 0; o < 8 && oid < 0; o++) {
      if (mode_of_oid(o) != r_mo) continue;
      for (int c = 0; c < 3; c++) {
        if (kRootDigit[o][c] == digits[0]) { oid = o; root_c2 = c; }
      }
    }
    if (oid < 0) continue;
    for (int tail_c2 = 0; tail_c2 < 3; tail_c2++) {
      int rid[kMaxLevel + 2];
      int c_mo = 0, c2 = tail_c2;
      for (int l = level; l >= 1; l--) {
        const uint8_t* e = kHexReg[digits[l]][c_mo][c2];
        rid[l] = e[0];
        c_mo = e[1];
        c2 = e[2];
      }
      if (c2 != root_c2) continue;
      rid[0] = r_mo;
      // The t_cell the cell's mode-0 half lies in, and the cell's centre on
      // its edge labelled tail_c2
      int64_t ou = 0, ov = 0;
      int p_mo = r_mo;
      bool valid = true;
      for (int k = 1; k <= level; k++) {
        const int j = child_index(p_mo, kRidCode[rid[k]]);
        if (j < 0) { valid = false; break; }
        const int64_t child = pow3(level - k);
        ou += kChildOffset[p_mo][j][0] * child;
        ov += kChildOffset[p_mo][j][1] * child;
        p_mo = kChildMode[p_mo][j];
      }
      if (!valid || p_mo != 0) continue;
      int64_t a[3];
      chart_weights(ou + kHalfCentre[0][tail_c2][0], ov + kHalfCentre[0][tail_c2][1],
                    level, r_mo, a);
      if (a[0] < 0 || a[1] < 0 || a[2] < 0) continue;
      OctPoint p;
      for (int k = 0; k < 3; k++) p.c[k] = ((oid >> k) & 1) ? -a[k] : a[k];
      if (encode(p, level) != id) continue;
      centre = p;
      return true;
    }
  }
  return false;
}

Tables tables() {
  Tables t;
  for (int m = 0; m < 2; m++) {
    for (int j = 0; j < 9; j++) {
      t.child_code.push_back(kChildCode[m][j]);
      t.child_mode.push_back(kChildMode[m][j]);
      t.child_offset.push_back(kChildOffset[m][j][0]);
      t.child_offset.push_back(kChildOffset[m][j][1]);
      t.root_c2.push_back(kRootC2[m][j]);
    }
  }
  for (int r = 0; r < 12; r++) t.rid_code.push_back(kRidCode[r]);
  for (int o = 0; o < 8; o++) {
    for (int c = 0; c < 3; c++) t.root_digit.push_back(kRootDigit[o][c]);
  }
  for (int g = 0; g < 2; g++) {
    for (int p = 0; p < 12; p++) {
      for (int c = 0; c < 12; c++) {
        t.reg_hex.push_back(kRegHex[g][p][c][0]);
        t.reg_hex.push_back(kRegHex[g][p][c][1]);
      }
    }
  }
  return t;
}

bool digit_strings_unique() {
  // A chain is in state (m, r) after naming t_cell r, whose parent has mode
  // m; the next t_cell r' of a valid child emits the digit of the cell whose
  // half-hexagon in r holds r'. A pair of chains reading the same digits
  // tracks whether they have named different t_cells (or octants) yet. Two
  // canonical cells of one level share their digits exactly when such a pair,
  // diverged and both at a mode-0 t_cell, can emit the same last digit.
  struct Start { int oid, state, digit; };
  std::vector<Start> starts;
  for (int o = 0; o < 8; o++) {
    const int m = mode_of_oid(o);
    for (int j = 0; j < 9; j++) {
      const int r = rid_of_code(kChildCode[m][j]);
      starts.push_back({o, m * 12 + r, kRootDigit[o][kRootC2[m][j]]});
    }
  }
  auto next = [](int state, int j, int& to, int& digit) {
    const int m = state / 12, r = state % 12;
    const int r2 = rid_of_code(kChildCode[r & 1][j]);
    const uint8_t* e = kRegHex[m][r][r2];
    if (e[0] == 255) throw std::logic_error("hex9: digit table misses a child");
    to = (r & 1) * 12 + r2;
    digit = e[0];
  };
  std::vector<char> seen(24 * 24 * 2, 0);
  std::vector<int> queue;
  auto push = [&](int a, int b, int div) {
    const int key = (a * 24 + b) * 2 + div;
    if (!seen[key]) { seen[key] = 1; queue.push_back(key); }
  };
  for (const Start& s1 : starts) {
    for (const Start& s2 : starts) {
      if (s1.digit != s2.digit) continue;
      push(s1.state, s2.state, (s1.oid != s2.oid || s1.state != s2.state) ? 1 : 0);
    }
  }
  for (size_t q = 0; q < queue.size(); q++) {
    const int key = queue[q];
    const int div = key % 2, a = (key / 2) / 24, b = (key / 2) % 24;
    int ta[9], da[9], tb[9], db[9];
    for (int j = 0; j < 9; j++) {
      next(a, j, ta[j], da[j]);
      next(b, j, tb[j], db[j]);
    }
    const bool mode0 = (a % 12) % 2 == 0 && (b % 12) % 2 == 0;
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        if (da[i] != db[j]) continue;
        if (div && mode0) return false;
        push(ta[i], tb[j], (div || ta[i] % 12 != tb[j] % 12) ? 1 : 0);
      }
    }
  }
  return true;
}

OctPoint ancestor(const OctPoint& centre, int level, int to_level) {
  check_level(level);
  if (to_level < 0 || to_level >= level) {
    throw std::invalid_argument("hex9: an ancestor is of a coarser level, 0 or finer");
  }
  return chain_centre(cell_chain(centre, level), to_level);
}

void children(const OctPoint& centre, int level, OctPoint out[9]) {
  check_level(level + 1);
  // The cell's triangles one level down whose t_cell has mode 0 hold the
  // mode-0 halves of its nine children, three each.
  const int64_t n = pow3(level + 1);
  OctantPoint oct[2];
  const int n_oct = octants_of(centre, oct);
  int found = 0;
  for (int o = 0; o < n_oct; o++) {
    const int oid = oid_of(oct[o].s);
    const int oct_mode = mode_of_oid(oid);
    for (int d = 0; d < 6; d++) {
      int64_t c[3];
      bool inside = true;
      for (int k = 0; k < 3; k++) {
        c[k] = 3 * oct[o].a[k] + kAround[d][k];
        if (c[k] <= 0) inside = false;
      }
      if (!inside || tcell_mode(c, 3, n, oct_mode) != 0) continue;
      int64_t u, v;
      chart_uv(c, level + 1, oct_mode, u, v);
      for (int c2 = 0; c2 < 3; c2++) {
        if (found == 9) throw std::logic_error("hex9: a cell with more than nine children");
        int64_t a[3];
        chart_weights(u + kHalfCentre[0][c2][0], v + kHalfCentre[0][c2][1], level + 1,
                      oct_mode, a);
        OctPoint p;
        for (int k = 0; k < 3; k++) p.c[k] = oct[o].s[k] * a[k];
        out[found++] = p;
      }
    }
  }
  if (found != 9) throw std::logic_error("hex9: a cell without nine children");
}

} // namespace hex9
} // namespace hexify
