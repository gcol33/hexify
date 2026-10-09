// globe_table.cpp
// The values of a grid's cells as the globe shader reads them
//
// Every resolution of the grid, from its own down to 0, lays its cells out
// as one block of the table per diamond quad, in that quad's substrate
// coordinates (Jubair et al. 2016 store icosahedral data as ten such
// diamonds). Row u of a block holds the cells of substrate column u, the
// lattice points (u, v) with v = c * u (mod index), in the order of v; the
// rows and columns run past the quad's box on every side, and a slot there
// holds the cell that owns its point across the quad's edge. The shader then
// reads the cell nearest a point, and the corners of the lattice triangle
// around it, at the place it finds them in the quad of the point's face,
// without moving them into the quad that owns them.
//
// The table is a list of 32-bit words. When the cells given fill much of the
// blocks, a block slot's word sits at the slot's own place in the list.
// Otherwise the slots of the given cells are packed by a perfect hash
// (Lefebvre and Hoppe 2006): slot key k goes to (k' mod M + offset[h(k)])
// mod M, with k' the low 32 bits of k, where h hashes k into a much smaller
// table of offsets chosen so that no two keys meet, and the key is stored
// beside its word.
//
// A coarser level's cell takes the finer cells whose centres lie in it. Laid
// out slot by slot, it reads them from the finer level's block around its
// own centre through a fixed stencil (gather_level()); packed by the hash,
// each given finer cell adds itself to the coarser cells nearest it
// (push_level()). The two agree cell by cell (tests/testthat/test-globe.R).
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include "cell_id.h"
#include "coordinate_transforms.h"
#include "globe_mesh.h"
#include "globe_table.h"
#include "grid_math.h"
#include "polyhedron.h"
#include "rcpp_icosa.h"

using namespace Rcpp;

namespace {

using hexify::GlobeFrame;

constexpr uint32_t kAbsent = 0xFFFFFFFFu;   // no cell given at this slot
constexpr uint32_t kNA = 0x7FC00000u;       // a cell whose value is NA
constexpr uint64_t kSide = 8192;            // a texture side every WebGPU device takes,
                                            // 2^13: a slot s is texel (s mod 2^13,
                                            // s / 2^13 mod 2^13) of layer s / 2^26
constexpr uint64_t kLayers = 256;           // texture array layers likewise

// A coarse cell is kept while the cells given cover at least this share of
// it; below, its fill would be at most this opaque, and a lone cell would
// cost an entry at every coarser level.
constexpr double kLeastCover = 1.0 / 4;

// A resolution's part of the table. Each diamond quad's block has hp rows
// of wp slots: row u + pr holds substrate column u, slot (v - residue(u)) /
// index + pc of it the lattice point (u, v). The block starts `base` slots
// into the table, quads in order.
struct Level {
  GlobeFrame f;
  long long pr, pc, hp, wp;
  uint64_t base;
};

// A cell: its quad and substrate (u, v); a vertex quad's single cell is
// (0, 0) of it.
struct Cell {
  int quad;
  long long u, v;
};

inline long long residue(const GlobeFrame& f, long long u) {
  if (f.index == 1) return 0;
  return ((f.c * u) % f.index + f.index) % f.index;
}

// x / d and x % d for x < 2^53, through the reciprocal of d, corrected by
// one where the product's rounding crossed a multiple of d
struct Divider {
  uint64_t d = 1;
  double inv = 1.0;
  explicit Divider(uint64_t den = 1) : d(den), inv(1.0 / static_cast<double>(den)) {}
  uint64_t div(uint64_t x, uint64_t& rem) const {
    uint64_t q = static_cast<uint64_t>(static_cast<double>(x) * inv);
    if (q * d > x) q--;
    else if ((q + 1) * d <= x) q++;
    rem = x - q * d;
    return q;
  }
};

// floor(v / index): a lattice point's place along its row of the quad,
// (v - residue(u)) / index, for any lattice point (u, v)
inline long long row_place(const GlobeFrame& f, long long v) {
  return f.index == 1 ? v : hexify::floor_div(v, f.index);
}

// A lattice point of a quad's box and its padding sits this far past the
// box, along either axis, from the centre of any cell whose triangle of
// neighbouring centres it is a corner of: a point of the box lies within a
// generator's length |g| = sqrt(index) of those corners, 2 |g| / sqrt(3)
// steps of either skew axis; one more row covers the far side's first row
// past the box, one more a point a rounding error outside.
long long triangle_reach(const GlobeFrame& f) {
  return static_cast<long long>(
      std::floor(2.0 * std::sqrt(static_cast<double>(f.index)) / std::sqrt(3.0))) + 2;
}

// A level's frame with padding `reach` past its box, kept within one quad of
// it, which a single crossing of an edge reaches.
Level make_level(const GlobeFrame& f, long long reach, uint64_t base) {
  Level L;
  L.f = f;
  L.pr = std::min(f.dim, reach);
  L.pc = (L.pr + 2 * (f.index - 1)) / f.index;
  L.hp = f.dim + 2 * L.pr;
  L.wp = f.dim / f.index + 2 * L.pc;
  L.base = base;
  return L;
}

uint64_t level_size(const Level& L) {
  return static_cast<uint64_t>(hexify::topo().n_diamonds()) * L.hp * L.wp;
}

// The 0-based index of a cell, its ID - 1, in the shader's numbering
uint64_t cell_index(const GlobeFrame& f, const Cell& x) {
  const hexify::SolidTopology& t = hexify::topo();
  if (x.quad == 0) return 0;
  if (x.quad == t.south_pole()) return f.n_cells - 1;
  return 1 + static_cast<uint64_t>(x.quad - 1) * f.per_quad +
         static_cast<uint64_t>(x.u) * (f.dim / f.index) +
         static_cast<uint64_t>(row_place(f, x.v));
}

// The cells of a frame by index, each read once: an index one past the last
// is a step along its row, wrapping into the next row and quad; any other
// index is decoded.
class Walk {
 public:
  explicit Walk(const GlobeFrame& f)
      : f_(f), per_row_(f.dim / f.index), south_(hexify::topo().south_pole()),
        quad_div_(f.per_quad), row_div_(static_cast<uint64_t>(f.dim / f.index)) {}

  const Cell& at(uint64_t i) {
    if (i == idx_ + 1 && x_.quad != 0 && x_.quad != south_ && i != f_.n_cells - 1) {
      col_++;
      x_.v += f_.index;
      if (col_ == per_row_) {
        col_ = 0;
        if (++x_.u == f_.dim) {
          x_.u = 0;
          x_.quad++;
        }
        x_.v = residue(f_, x_.u);
      }
    } else {
      decode(i);
    }
    idx_ = i;
    return x_;
  }
  // The place of the last cell along its row of the quad
  long long col() const { return col_; }

 private:
  void decode(uint64_t i) {
    col_ = 0;
    if (i == 0) {
      x_ = {0, 0, 0};
      return;
    }
    if (i == f_.n_cells - 1) {
      x_ = {south_, 0, 0};
      return;
    }
    uint64_t within, col;
    x_.quad = static_cast<int>(quad_div_.div(i - 1, within)) + 1;
    x_.u = static_cast<long long>(row_div_.div(within, col));
    col_ = static_cast<long long>(col);
    x_.v = col_ * f_.index + residue(f_, x_.u);
  }

  const GlobeFrame& f_;
  long long per_row_;
  int south_;
  Divider quad_div_, row_div_;
  uint64_t idx_ = ~uint64_t(0);
  Cell x_ = {0, 0, 0};
  long long col_ = 0;
};

// The slot of point (u, v) of diamond quad q in level L's part; false off
// the padded block or off the lattice.
inline bool slot_of(const Level& L, int q, long long u, long long v, uint64_t& key) {
  const long long row = u + L.pr;
  if (row < 0 || row >= L.hp) return false;
  if (L.f.index != 1 && ((v - residue(L.f, u)) % L.f.index) != 0) return false;
  const long long col = row_place(L.f, v) + L.pc;
  if (col < 0 || col >= L.wp) return false;
  key = L.base + (static_cast<uint64_t>(q - 1) * L.hp + row) * L.wp + col;
  return true;
}

// A point (u, v) of a quad from its place across edge e: `along` the edge
// and d past it along the crossed axis.
inline void across(int e, long long top, long long along, long long d,
                   long long& u, long long& v) {
  switch (e) {
    case hexify::kEdgeLeft:  u = d;       v = along;   break;
    case hexify::kEdgeDown:  u = along;   v = d;       break;
    case hexify::kEdgeRight: u = top + d; v = along;   break;
    default:                 u = along;   v = top + d; break;
  }
}

// How far past an edge a level's padding reaches along either axis
inline long long strip_of(const Level& L) {
  return std::max(L.pr, L.pc * L.f.index + L.f.index);
}

// Every padding slot of level L that holds cell x: each padding point the
// quad edge maps (canonicalize_q2d(), as the shader's canonicalize() reads
// them) carry into it. A point past one edge maps across it, or to the
// vertex quad where a far edge starts when along = 0; a point past both far
// edges is the far corner's cell; past any other two edges no cell owns it.
// The edge maps are lattice isometries, so each preimage is solved for
// exactly. Only a cell within strip_of() of an edge of its box has any.
template <typename Emit>
void padding_slots(const Level& L, const Cell& x, Emit emit) {
  const hexify::SolidTopology& t = hexify::topo();
  const long long top = L.f.dim;
  const bool vertex_quad = t.is_pole(x.quad);
  const long long strip = strip_of(L);
  uint64_t key;
  for (int a = 1; a <= t.n_diamonds(); a++) {
    if (t.corner[a][hexify::kCornerFar] == x.quad && x.u == 0 && x.v == 0) {
      for (long long u = top; u < top + L.pr; u++) {
        for (long long v = top; v < top + strip; v++) {
          if (slot_of(L, a, u, v, key)) emit(key);
        }
      }
    }
    for (int e = 0; e < 4; e++) {
      const hexify::QuadEdgeMap& m = t.edge[a][e];
      const bool far = e == hexify::kEdgeRight || e == hexify::kEdgeUp;
      if (vertex_quad) {
        if (m.pole != x.quad) continue;
        for (long long d = far ? 0 : -strip; d < (far ? strip : 0); d++) {
          long long u, v;
          across(e, top, 0, d, u, v);
          if (slot_of(L, a, u, v, key)) emit(key);
        }
        continue;
      }
      if (m.quad != x.quad) continue;
      const long long det = static_cast<long long>(m.k[0][1]) * m.k[1][2] -
                            static_cast<long long>(m.k[0][2]) * m.k[1][1];
      if (det != 1 && det != -1) Rcpp::stop("hexify internal error: a quad edge map is not unimodular");
      const long long di = x.u - m.k[0][0] * top, dj = x.v - m.k[1][0] * top;
      const long long along = (m.k[1][2] * di - m.k[0][2] * dj) / det;
      const long long d = (m.k[0][1] * dj - m.k[1][1] * di) / det;
      if (along < 0 || along >= top || (m.pole >= 0 && along == 0)) continue;
      if (far ? d < 0 : d >= 0) continue;
      long long u, v;
      across(e, top, along, d, u, v);
      if (slot_of(L, a, u, v, key)) emit(key);
    }
  }
}

inline bool near_edge(const Level& L, const Cell& x) {
  const long long top = L.f.dim;
  return std::min(std::min(x.u, x.v), std::min(top - 1 - x.u, top - 1 - x.v)) <= strip_of(L);
}

// The own slot of the cell of diamond quad q at row u, `col` along it
inline uint64_t own_slot(const Level& L, int q, long long u, long long col) {
  return L.base + (static_cast<uint64_t>(q - 1) * L.hp + u + L.pr) * L.wp + col + L.pc;
}

// A cell's own slot, `col` its place along its row, and near an edge its
// padding slots
template <typename Emit>
void cell_slots(const Level& L, const Cell& x, long long col, Emit emit) {
  if (hexify::topo().is_pole(x.quad)) {
    padding_slots(L, x, emit);
    return;
  }
  emit(own_slot(L, x.quad, x.u, col));
  if (near_edge(L, x)) padding_slots(L, x, emit);
}

// The word a cell is stored as. At the grid's resolution, its value as a
// 32-bit float, NA as one NaN. At a coarser one, where the shader reads the
// fill only, its place on the colour ramp in the high 16 bits (0xffff for
// NA) and its cover in the low 16, both in steps of 1 / 65534, so that no
// word is ABSENT.
uint32_t value_word(double v) {
  if (std::isnan(v)) return kNA;
  const float f = static_cast<float>(v);
  uint32_t w;
  std::memcpy(&w, &f, 4);
  return w;
}

uint32_t coarse_word(double pos, double cover) {
  const uint32_t p = std::isnan(pos) ? 0xffffU : static_cast<uint32_t>(pos * 65534.0 + 0.5);
  const uint32_t c = std::max(1U, static_cast<uint32_t>(cover * 65534.0 + 0.5));
  return (p << 16) | c;
}

// A stored cell as a coarser level reads it: its place on the ramp (NaN for
// NA) and its cover
struct Sample {
  double pos, cover;
};

struct Ramp {
  double at[3];
  double place(double v) const {
    return std::min(std::max((v - at[0]) * at[1] + at[2], 0.0), 1.0);
  }
  Sample read(uint32_t word, bool finest) const {
    if (finest) {
      if ((word & 0x7fffffffU) > 0x7f800000U) return {NAN, 1.0};
      float f;
      std::memcpy(&f, &word, 4);
      return {place(f), 1.0};
    }
    const uint32_t p = word >> 16;
    return {p == 0xffffU ? NAN : p / 65534.0, (word & 0xffffU) / 65534.0};
  }
};

// The share of a coarse cell's area a fine cell's worth is: a step of
// aperture a gives a cell with n neighbours 1 + n (a - 1) / 6 children's
// worth, a for a hexagon, less at a vertex of the solid
inline double whole(const hexify::SolidTopology& t, double step, const Cell& y) {
  const bool vertex = t.is_pole(y.quad) || (y.u == 0 && y.v == 0);
  return vertex ? 1.0 + t.valence[y.quad] * (step - 1.0) / 6.0 : step;
}

// The aperture of the step from `fine` to `coarse`
inline double step_of(const GlobeFrame& fine, const GlobeFrame& coarse) {
  return static_cast<double>(fine.per_quad / coarse.per_quad);
}

// A coarse cell from the sums of its children: kept when they cover at least
// kLeastCover of it, its place on the ramp the mean of theirs weighted by the
// area they cover
struct Acc {
  double cover = 0, weight = 0, sum = 0;
  void add(const Sample& s, double share) {
    cover += share * s.cover;
    if (!std::isnan(s.pos)) {
      weight += share * s.cover;
      sum += share * s.cover * s.pos;
    }
  }
  bool kept(double whole_share, uint32_t& word) const {
    const double c = cover / whole_share;
    if (c < kLeastCover) return false;
    word = coarse_word(weight > 0 ? sum / weight : NAN, std::min(c, 1.0));
    return true;
  }
};

// The fine cells whose centres lie in a coarse cell, as substrate offsets
// from the coarse centre in the fine frame, each with its share: a centre on
// the coarse cell's edge or corner splits evenly among the cells meeting
// there (nearest_eisenstein_points()). Coarse centres are fine lattice points
// and both lattices are translation invariant, so one stencil serves every
// coarse cell. `reach` bounds the offsets along either axis.
struct Stencil {
  std::vector<long long> du, dv;
  std::vector<double> share;
  long long reach = 0;
};

// The coarse cells nearest the fine point (u, v) of a quad, in the coarse
// substrate of that quad, with the share of each
inline int coarse_nearest(const GlobeFrame& fine, const GlobeFrame& coarse,
                          long long u, long long v, long long out[3][2]) {
  const long long k = fine.dim / coarse.dim;
  const long long a = coarse.ga, b = coarse.gb;
  const long long D = k * hexify::eisenstein_norm(a - b, b);
  // (u + v*omega) * ((a - b) - b*omega), using omega^2 = -1 - omega
  long long near[3][2];
  const int n = hexify::nearest_eisenstein_points(u * (a - b) + v * b, -u * b + v * a, D, near);
  for (int m = 0; m < n; m++) {
    out[m][0] = near[m][0] * a - near[m][1] * b;
    out[m][1] = near[m][0] * b + near[m][1] * a - near[m][1] * b;
  }
  return n;
}

Stencil child_stencil(const GlobeFrame& fine, const GlobeFrame& coarse) {
  Stencil s;
  const long long k = fine.dim / coarse.dim;
  const long long R = static_cast<long long>(
      std::ceil(2.0 * k * std::sqrt(static_cast<double>(coarse.index)) / 3.0)) + 1;
  for (long long du = -R; du <= R; du++) {
    for (long long dv = -R; dv <= R; dv++) {
      if (residue(fine, du) != ((dv % fine.index) + fine.index) % fine.index) continue;
      long long near[3][2];
      const int n = coarse_nearest(fine, coarse, du, dv, near);
      for (int m = 0; m < n; m++) {
        if (near[m][0] != 0 || near[m][1] != 0) continue;
        s.du.push_back(du);
        s.dv.push_back(dv);
        s.share.push_back(1.0 / n);
        s.reach = std::max(s.reach, std::max(std::llabs(du), std::llabs(dv)));
      }
    }
  }
  return s;
}

// Every cell of `coarse` from the words of `fine` laid out slot by slot:
// each cell sums its stencil's slots around its centre, in its own quad. A
// cell at a vertex of the solid (a diamond quad's origin, or a vertex quad's
// cell) is seen whole from no single quad, the solid folding there: it sums
// the stencil's slots inside the box of each quad it is a corner of, which
// hold each child once, and the cell at its vertex, which a vertex quad
// holds outside every box. Needs the fine padding to reach the stencil.
template <typename Keep>
void gather_level(const Level& fine, const Level& coarse, const Stencil& st,
                  const std::vector<uint32_t>& words, const Ramp& ramp, bool fine_finest,
                  Keep keep) {
  const hexify::SolidTopology& t = hexify::topo();
  const long long k = fine.f.dim / coarse.f.dim;
  const long long top = fine.f.dim;
  const long long nf = fine.f.index;
  const double step = step_of(fine.f, coarse.f);
  const size_t n = st.du.size();
  auto read = [&](uint64_t key, double share, Acc& acc) {
    const uint32_t w = words[key];
    if (w != kAbsent) acc.add(ramp.read(w, fine_finest), share);
  };
  // A stencil point's slot from the slot of the centre, for each residue
  // r = pv mod index of the centre's v: rows du apart, and places along the
  // row floor((r + dv) / index) apart
  std::vector<long long> delta(static_cast<size_t>(nf) * n);
  for (long long r = 0; r < nf; r++) {
    for (size_t m = 0; m < n; m++) {
      delta[r * n + m] = st.du[m] * fine.wp + hexify::floor_div(r + st.dv[m], nf);
    }
  }
  const long long per_row = coarse.f.dim / coarse.f.index;
  for (int q = 1; q <= t.n_diamonds(); q++) {
    const uint64_t block = fine.base + static_cast<uint64_t>(q - 1) * fine.hp * fine.wp;
    for (long long u = 0; u < coarse.f.dim; u++) {
      const long long r = residue(coarse.f, u);
      const long long pu = k * u;
      const uint64_t row = block + static_cast<uint64_t>(pu + fine.pr) * fine.wp + fine.pc;
      for (long long col = (u == 0 ? 1 : 0); col < per_row; col++) {
        const Cell y = {q, u, col * coarse.f.index + r};
        const long long pv = k * y.v;
        const uint64_t centre = row + static_cast<uint64_t>(pv / nf);
        const long long* d = &delta[(pv % nf) * n];
        Acc acc;
        for (size_t m = 0; m < n; m++) read(centre + d[m], st.share[m], acc);
        keep(y, col, acc, whole(t, step, y));
      }
    }
  }
  for (int vertex = 0; vertex < t.n_quads(); vertex++) {
    Acc acc;
    for (int a = 1; a <= t.n_diamonds(); a++) {
      for (int c = 0; c < 4; c++) {
        if (t.corner[a][c] != vertex) continue;
        const long long pu = (c == hexify::kCornerI || c == hexify::kCornerFar) ? top : 0;
        const long long pv = (c == hexify::kCornerJ || c == hexify::kCornerFar) ? top : 0;
        for (size_t m = 0; m < n; m++) {
          const long long u = pu + st.du[m], v = pv + st.dv[m];
          uint64_t key;
          if (u >= 0 && v >= 0 && u < top && v < top && slot_of(fine, a, u, v, key)) {
            read(key, st.share[m], acc);
          }
        }
      }
    }
    if (t.is_pole(vertex)) {
      const Cell x = {vertex, 0, 0};
      bool done = false;
      padding_slots(fine, x, [&](uint64_t key) {
        if (!done) read(key, 1.0, acc);
        done = true;
      });
    }
    const Cell y = {vertex, 0, 0};
    keep(y, 0, acc, whole(t, step, y));
  }
}

// A cell given at a level: its index and how a coarser level reads it
struct Given {
  uint64_t idx;
  Sample s;
};

// The cells of `coarse` from the cells given at `fine`, each adding itself
// to the coarse cells nearest its centre, a share 1 / ties of it each:
// (coarse cell, fine cell, ties) records, sorted by coarse cell and summed
// run by run.
template <typename Keep>
void push_level(const GlobeFrame& fine, const GlobeFrame& coarse,
                const std::vector<Given>& cells, Keep keep) {
  const hexify::SolidTopology& t = hexify::topo();
  const double step = step_of(fine, coarse);
  struct Share {
    uint64_t idx;
    uint32_t child;
    uint32_t ties;
  };
  std::vector<Share> shares;
  shares.reserve(cells.size() + cells.size() / 2);
  Walk walk(fine);
  for (size_t c = 0; c < cells.size(); c++) {
    const Cell& x = walk.at(cells[c].idx);
    const uint32_t child = static_cast<uint32_t>(c);
    if (t.is_pole(x.quad)) {
      shares.push_back({x.quad == 0 ? 0 : coarse.n_cells - 1, child, 1});
      continue;
    }
    long long near[3][2];
    const int n = coarse_nearest(fine, coarse, x.u, x.v, near);
    for (int m = 0; m < n; m++) {
      Cell y = {x.quad, near[m][0], near[m][1]};
      const bool inside = y.u >= 0 && y.v >= 0 && y.u < coarse.dim && y.v < coarse.dim;
      if (!inside && !hexify::substrate_ij_canonicalize(y.quad, y.u, y.v, coarse.dim)) continue;
      shares.push_back({cell_index(coarse, y), child, static_cast<uint32_t>(n)});
    }
  }
  // Least significant digit first radix sort by coarse cell, 11 bits a pass
  {
    uint64_t largest = 0;
    for (const Share& sh : shares) largest = std::max(largest, sh.idx);
    std::vector<Share> other(shares.size());
    for (int shift = 0; shift < 64 && (largest >> shift) != 0; shift += 11) {
      size_t count[2049] = {0};
      for (const Share& sh : shares) count[((sh.idx >> shift) & 2047) + 1]++;
      for (int b = 0; b < 2048; b++) count[b + 1] += count[b];
      for (const Share& sh : shares) other[count[(sh.idx >> shift) & 2047]++] = sh;
      shares.swap(other);
    }
  }
  // No cell of the solid is worth less than one at a vertex of valence 3
  const double least = kLeastCover * (1.0 + 3.0 * (step - 1.0) / 6.0);
  Walk coarse_walk(coarse);
  for (size_t m = 0; m < shares.size();) {
    Acc acc;
    size_t e = m;
    for (; e < shares.size() && shares[e].idx == shares[m].idx; e++) {
      acc.add(cells[shares[e].child].s, 1.0 / shares[e].ties);
    }
    if (acc.cover >= least) {
      const Cell& y = coarse_walk.at(shares[m].idx);
      keep(y, coarse_walk.col(), acc, whole(t, step, y));
    }
    m = e;
  }
}

uint32_t mix32(uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

// The offset a key reads; globe.wgsl's table_word() is the same
uint32_t bucket_of(uint64_t key, uint32_t seed, uint32_t n_buckets) {
  return mix32(static_cast<uint32_t>(key) ^
               mix32(static_cast<uint32_t>(key >> 32) ^ seed)) % n_buckets;
}

struct Entry {
  uint64_t key;
  uint32_t word;
};

// Offsets for `n_buckets` buckets placing every key at its own one of M
// slots, or false. Buckets are placed largest first, each at an offset that
// lands all its keys on free slots: a random free slot for its first key,
// tried until the others are free too, and for a single key the next free
// slot.
bool place(const std::vector<Entry>& e, uint32_t M, uint32_t n_buckets, uint32_t seed,
           std::vector<uint32_t>& offset, std::vector<uint32_t>& slot) {
  const size_t E = e.size();
  std::vector<uint32_t> bucket(E), start(n_buckets + 1, 0);
  for (size_t k = 0; k < E; k++) {
    bucket[k] = bucket_of(e[k].key, seed, n_buckets);
    start[bucket[k] + 1]++;
  }
  uint32_t largest = 0;
  for (uint32_t b = 0; b < n_buckets; b++) {
    largest = std::max(largest, start[b + 1]);
    start[b + 1] += start[b];
  }
  std::vector<uint32_t> member(E), fill(start.begin(), start.end() - 1);
  for (size_t k = 0; k < E; k++) member[fill[bucket[k]]++] = static_cast<uint32_t>(k);

  std::vector<uint32_t> by_size_start(largest + 2, 0), order(n_buckets);
  for (uint32_t b = 0; b < n_buckets; b++) by_size_start[largest - (start[b + 1] - start[b]) + 1]++;
  for (uint32_t s = 0; s <= largest; s++) by_size_start[s + 1] += by_size_start[s];
  for (uint32_t b = 0; b < n_buckets; b++) {
    order[by_size_start[largest - (start[b + 1] - start[b])]++] = b;
  }

  std::vector<uint8_t> used(M, 0);
  uint64_t state = 0x9E3779B97F4A7C15ULL ^ seed;
  auto next = [&state]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
  };
  offset.assign(n_buckets, 0);
  slot.assign(E, 0);
  uint32_t cursor = 0;
  std::vector<uint32_t> base(largest);
  for (uint32_t b : order) {
    const uint32_t first = start[b], size = start[b + 1] - first;
    if (size == 0) break;
    for (uint32_t m = 0; m < size; m++) {
      base[m] = static_cast<uint32_t>(e[member[first + m]].key) % M;
    }
    uint32_t o = 0;
    if (size == 1) {
      while (used[cursor]) cursor++;
      o = cursor >= base[0] ? cursor - base[0] : cursor + M - base[0];
    } else {
      for (uint32_t m = 1; m < size; m++) {
        for (uint32_t l = 0; l < m; l++) {
          if (base[l] == base[m]) return false;
        }
      }
      bool ok = false;
      for (int tries = 0; tries < (1 << 22) && !ok; tries++) {
        const uint32_t s0 = static_cast<uint32_t>(((next() >> 32) * M) >> 32);
        if (used[s0]) continue;
        o = s0 >= base[0] ? s0 - base[0] : s0 + M - base[0];
        ok = true;
        for (uint32_t m = 1; m < size && ok; m++) {
          const uint32_t s = base[m] + o;
          ok = !used[s >= M ? s - M : s];
        }
      }
      if (!ok) return false;
    }
    offset[b] = o;
    for (uint32_t m = 0; m < size; m++) {
      const uint32_t s = base[m] + o >= M ? base[m] + o - M : base[m] + o;
      used[s] = 1;
      slot[member[first + m]] = s;
    }
  }
  return true;
}

// The side lengths of a texture array holding n words in rows of kSide,
// layers of kSide rows, the shader finding a word by shifts alone
void texture_shape(uint64_t n, uint32_t& w, uint32_t& h, uint32_t& layers) {
  n = std::max<uint64_t>(1, n);
  w = static_cast<uint32_t>(std::min(n, kSide));
  const uint64_t rows = (n + kSide - 1) / kSide;
  const uint64_t l = (rows + kSide - 1) / kSide;
  if (l > kLayers) {
    Rcpp::stop("hex_globe() holds at most %.0f values per grid on a graphics card",
               static_cast<double>(kSide * kSide * kLayers));
  }
  layers = static_cast<uint32_t>(l);
  h = static_cast<uint32_t>(l == 1 ? rows : kSide);
}

} // namespace

// The table of a grid's cells for the globe shader. `levels` are the grid's
// resolutions from its own down, each list(resolution, aperture, ap_seq) as
// the C++ entry points take them. With `all_cells`, `values` holds one value
// per cell of the grid in ID order; otherwise the cells are `cell_id`, with
// `values` one per cell or empty when the cells carry none. `ramp` places a
// value v on the colour ramp at clamp((v - ramp[0]) * ramp[1] + ramp[2], 0, 1),
// as ramp_map() in R; coarser levels hold places, not values. `layout` is
// "auto", or "slots" or "hash" to choose the layout; with `listing`, the
// result also lists each level's cells and words.
//
// Returns the levels' frames and blocks (one row each) and the table:
// `values`, the word of each slot (base64 of little-endian 32-bit words),
// laid out as a texture array of `size` = width, height, layers; with
// `keyed`, the perfect hash's `keys` (two words per slot, high first) and
// `offsets` (a texture of `offsets_size`), read with `seed`, `buckets`
// buckets and `m` slots.
// [[Rcpp::export]]
List cpp_globe_table(NumericVector icosa, List levels, NumericVector cell_id,
                     NumericVector values, bool all_cells, NumericVector ramp,
                     std::string layout = "auto", bool listing = false) {
  activate_grid(icosa);
  const hexify::SolidTopology& t = hexify::topo();
  const int n_levels = levels.size();
  if (n_levels == 0) stop("cpp_globe_table needs a level");
  if (n_levels > 1 && ramp.size() != 3) stop("coarser levels need the colour ramp's map");
  Ramp rmap = {{ramp.size() == 3 ? ramp[0] : 0.0, ramp.size() == 3 ? ramp[1] : 0.0,
                ramp.size() == 3 ? ramp[2] : 0.0}};

  std::vector<GlobeFrame> frames;
  for (int k = 0; k < n_levels; k++) {
    List spec = levels[k];
    IntegerVector seq = spec["ap_seq"];
    frames.push_back(hexify::globe_frame(as<int>(spec["resolution"]), as<int>(spec["aperture"]),
                                         std::vector<int>(seq.begin(), seq.end())));
    if (k > 0 && (frames[k - 1].dim % frames[k].dim != 0 ||
                  frames[k - 1].per_quad % frames[k].per_quad != 0)) {
      stop("each level of the globe table must be the parent grid of the one before");
    }
  }
  const GlobeFrame& finest = frames[0];

  const uint64_t n0 = all_cells ? finest.n_cells : static_cast<uint64_t>(cell_id.size());
  if (all_cells && static_cast<uint64_t>(values.size()) != finest.n_cells) {
    stop("values must hold one value per cell of the grid");
  }
  const bool valued = values.size() > 0;
  if (!all_cells && valued && values.size() != cell_id.size()) {
    stop("values must hold one value per cell");
  }

  // Each level's stencil into the one finer, and the padding it needs
  std::vector<Stencil> stencils(n_levels);
  for (int k = 1; k < n_levels; k++) stencils[k] = child_stencil(frames[k - 1], frames[k]);
  std::vector<Level> lv;
  uint64_t total = 0;
  for (int k = 0; k < n_levels; k++) {
    long long reach = triangle_reach(frames[k]);
    if (k + 1 < n_levels) reach = std::max(reach, stencils[k + 1].reach + 1);
    lv.push_back(make_level(frames[k], reach, total));
    total += level_size(lv.back());
  }
  // Coarser levels add at most half the cells given (aperture 3)
  const bool keyed = layout == "hash" ||
      (layout != "slots" && (total >= (uint64_t(1) << 31) || total > 4.5 * n0 + 4096));
  if (!keyed && total >= (uint64_t(1) << 31)) {
    stop("hex_globe() holds at most 2^31 slots laid out one by one");
  }

  std::vector<std::vector<double>> list_idx;
  std::vector<std::vector<uint32_t>> list_word(n_levels);
  list_idx.resize(n_levels);
  auto note = [&](int k, const GlobeFrame& f, const Cell& x, uint32_t word) {
    if (listing) {
      list_idx[k].push_back(static_cast<double>(cell_index(f, x)));
      list_word[k].push_back(word);
    }
  };
  std::vector<int> given_n(n_levels, 0);

  uint32_t w, h, layers, seed = 0, n_buckets = 1, ow = 1, oh = 1, o_layers = 1;
  uint64_t m;
  std::vector<uint32_t> words, keys, offsets(1, 0);

  if (!keyed) {
    m = total;
    texture_shape(m, w, h, layers);
    words.assign(static_cast<size_t>(w) * h * layers, kAbsent);
    // A vertex quad's cell has no slot of its own, only padding slots
    std::vector<uint32_t> pole_word(2 * n_levels, kAbsent);
    auto pole_at = [&](int k, int quad) -> uint32_t& {
      return pole_word[2 * k + (quad == 0 ? 0 : 1)];
    };
    auto own = [&](int k, const Cell& x, long long col, uint32_t word) {
      uint32_t& slot = t.is_pole(x.quad) ? pole_at(k, x.quad)
                                         : words[own_slot(lv[k], x.quad, x.u, col)];
      if (slot != kAbsent) stop("cells must not repeat");
      slot = word;
    };
    // Level k's padding: each cell near an edge of its box, and each vertex
    // quad's cell, into the slots past the boxes that hold it
    auto pad = [&](int k) {
      const Level& L = lv[k];
      const long long strip = strip_of(L), top = L.f.dim, n = L.f.index;
      const long long per_row = top / n;
      auto spread = [&](const Cell& x, uint32_t word) {
        padding_slots(L, x, [&](uint64_t key) {
          if (words[key] != kAbsent) stop("cells must not repeat");
          words[key] = word;
        });
      };
      auto each = [&](int q, long long u, long long r, long long from, long long to) {
        for (long long col = std::max(from, 0LL); col < std::min(to, per_row); col++) {
          const uint32_t word = words[own_slot(L, q, u, col)];
          if (word != kAbsent) spread({q, u, col * n + r}, word);
        }
      };
      for (int q = 1; q <= t.n_diamonds(); q++) {
        for (long long u = 0; u < top; u++) {
          const long long r = residue(L.f, u);
          if (u <= strip || u >= top - 1 - strip) {
            each(q, u, r, 0, per_row);
            continue;
          }
          // v = col * n + r within strip of either end of the row
          const long long low = strip >= r ? (strip - r) / n + 1 : 0;
          const long long high = hexify::floor_div(top - 1 - strip - r + n - 1, n);
          if (low >= high) {
            each(q, u, r, 0, per_row);
          } else {
            each(q, u, r, 0, low);
            each(q, u, r, high, per_row);
          }
        }
      }
      for (int quad : {0, t.south_pole()}) {
        if (pole_at(k, quad) != kAbsent) spread({quad, 0, 0}, pole_at(k, quad));
      }
    };
    // The grid's resolution: each given cell at its own slot
    Walk walk(finest);
    for (uint64_t k = 0; k < n0; k++) {
      uint64_t idx = k;
      if (!all_cells) {
        const int64_t id = hexify::cell_id_get(cell_id[k]);
        if (id == hexify::kCellIdNA || id < 1 || static_cast<uint64_t>(id) > finest.n_cells) {
          stop("cell IDs must name cells of the grid");
        }
        idx = static_cast<uint64_t>(id - 1);
      }
      const Cell& x = walk.at(idx);
      const uint32_t word = valued ? value_word(values[k]) : 0U;
      own(0, x, walk.col(), word);
      note(0, finest, x, word);
    }
    given_n[0] = static_cast<int>(n0);
    pad(0);
    // Each coarser resolution from the one finer
    for (int k = 1; k < n_levels; k++) {
      auto keep = [&](const Cell& y, long long col, const Acc& acc, double whole_share) {
        uint32_t word;
        if (!acc.kept(whole_share, word)) return;
        own(k, y, col, word);
        note(k, frames[k], y, word);
        given_n[k]++;
      };
      if (lv[k - 1].pr >= stencils[k].reach + 1) {
        gather_level(lv[k - 1], lv[k], stencils[k], words, rmap, k == 1, keep);
      } else {
        // A level too small for its padding to reach the stencil
        std::vector<Given> fine;
        Walk fw(frames[k - 1]);
        for (uint64_t idx = 0; idx < frames[k - 1].n_cells; idx++) {
          const Cell& x = fw.at(idx);
          const uint32_t word = t.is_pole(x.quad) ? pole_at(k - 1, x.quad)
                                                  : words[own_slot(lv[k - 1], x.quad, x.u, fw.col())];
          if (word != kAbsent) fine.push_back({idx, rmap.read(word, k == 1)});
        }
        push_level(frames[k - 1], frames[k], fine, keep);
      }
      pad(k);
    }
  } else {
    std::vector<Entry> entries;
    auto put = [&](const Level& L, const Cell& x, long long col, uint32_t word) {
      cell_slots(L, x, col, [&](uint64_t key) { entries.push_back({key, word}); });
    };
    std::vector<Given> cells;
    cells.reserve(n0);
    Walk walk(finest);
    for (uint64_t k = 0; k < n0; k++) {
      uint64_t idx = k;
      if (!all_cells) {
        const int64_t id = hexify::cell_id_get(cell_id[k]);
        if (id == hexify::kCellIdNA || id < 1 || static_cast<uint64_t>(id) > finest.n_cells) {
          stop("cell IDs must name cells of the grid");
        }
        idx = static_cast<uint64_t>(id - 1);
      }
      const Cell& x = walk.at(idx);
      const uint32_t word = valued ? value_word(values[k]) : 0U;
      put(lv[0], x, walk.col(), word);
      note(0, finest, x, word);
      if (n_levels > 1) cells.push_back({idx, rmap.read(word, true)});
    }
    given_n[0] = static_cast<int>(n0);
    for (int k = 1; k < n_levels; k++) {
      std::vector<Given> next;
      push_level(frames[k - 1], frames[k], cells, [&](const Cell& y, long long col, const Acc& acc,
                                                      double ws) {
        uint32_t word;
        if (!acc.kept(ws, word)) return;
        put(lv[k], y, col, word);
        note(k, frames[k], y, word);
        next.push_back({cell_index(frames[k], y), rmap.read(word, false)});
      });
      given_n[k] = static_cast<int>(next.size());
      cells.swap(next);
    }
    const uint64_t E = entries.size();
    m = std::max<uint64_t>(1, E + E / 16 + 16);
    if (m >= (uint64_t(1) << 31)) stop("hex_globe() holds at most 2^31 cells on a graphics card");
    texture_shape(m, w, h, layers);
    n_buckets = static_cast<uint32_t>(std::max<uint64_t>(1, E / 2));
    std::vector<uint32_t> slot;
    bool placed = false;
    for (int attempt = 0; attempt < 40 && !placed; attempt++) {
      seed = mix32(0x9e3779b9U + static_cast<uint32_t>(attempt));
      placed = place(entries, static_cast<uint32_t>(m), n_buckets, seed, offsets, slot);
      if (!placed) n_buckets += n_buckets / 8 + 1;
    }
    if (!placed) stop("cells must not repeat");
    words.assign(static_cast<size_t>(w) * h * layers, kAbsent);
    keys.assign(2 * words.size(), kAbsent);
    for (size_t k = 0; k < E; k++) {
      words[slot[k]] = entries[k].word;
      keys[2 * slot[k]] = static_cast<uint32_t>(entries[k].key >> 32);
      keys[2 * slot[k] + 1] = static_cast<uint32_t>(entries[k].key);
    }
    texture_shape(n_buckets, ow, oh, o_layers);
    if (o_layers > 1) stop("hex_globe() holds at most %.0f cells on a graphics card",
                           static_cast<double>(3 * kSide * kSide));
    offsets.resize(static_cast<size_t>(ow) * oh, 0);
  }

  NumericMatrix out_levels(n_levels, 11);
  for (int k = 0; k < n_levels; k++) {
    const Level& L = lv[k];
    const double row[11] = {static_cast<double>(L.f.dim), static_cast<double>(L.f.index),
                            static_cast<double>(L.f.c), static_cast<double>(L.f.ga),
                            static_cast<double>(L.f.gb), static_cast<double>(L.pr),
                            static_cast<double>(L.pc), static_cast<double>(L.hp),
                            static_cast<double>(L.wp), static_cast<double>(L.base),
                            static_cast<double>(L.f.n_cells)};
    for (int c = 0; c < 11; c++) out_levels(k, c) = row[c];
  }
  colnames(out_levels) = CharacterVector::create("dim", "index", "c", "ga", "gb", "pr", "pc",
                                                 "hp", "wp", "base", "n_cells");

  List out = List::create(
      _["levels"] = out_levels,
      _["given"] = IntegerVector(given_n.begin(), given_n.end()),
      _["keyed"] = keyed,
      _["m"] = static_cast<double>(m),
      _["size"] = NumericVector::create(w, h, layers),
      _["values"] = hexify::base64_words(words),
      _["keys"] = keyed ? hexify::base64_words(keys) : std::string(),
      _["seed"] = static_cast<double>(seed),
      _["buckets"] = static_cast<double>(n_buckets),
      _["offsets_size"] = NumericVector::create(ow, oh),
      _["offsets"] = hexify::base64_words(offsets));
  if (listing) {
    List cells(n_levels);
    for (int k = 0; k < n_levels; k++) {
      NumericVector word(list_word[k].begin(), list_word[k].end());
      cells[k] = List::create(_["idx"] = NumericVector(list_idx[k].begin(), list_idx[k].end()),
                              _["word"] = word);
    }
    out["cells"] = cells;
  }
  return out;
}
