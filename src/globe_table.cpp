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
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <random>
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
constexpr uint64_t kSide = 8192;            // a texture side every WebGPU device takes
constexpr uint64_t kLayers = 256;           // texture array layers likewise

// A resolution's part of the table. Each diamond quad's block has hp rows
// of wp slots: row u + pr holds substrate column u, slot (v - residue(u)) /
// index + pc of it the lattice point (u, v). The block starts `base` slots
// into the resolution's part, quads in order.
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
  return ((f.c * u) % f.index + f.index) % f.index;
}

// A level's frame with its padding. A point of a quad's box lies within a
// generator's length |g| = sqrt(index) of the corners of the lattice
// triangle around it, which therefore lie within 2 |g| / sqrt(3) substrate
// steps of the box along either axis; one more row covers the far side's
// first row past the box, and one more a point a rounding error outside.
// The padding stays within one quad of the box, which a single crossing of
// an edge reaches.
Level make_level(const GlobeFrame& f, uint64_t base) {
  Level L;
  L.f = f;
  const long long reach =
      static_cast<long long>(std::floor(2.0 * std::sqrt(static_cast<double>(f.index)) /
                                        std::sqrt(3.0))) + 2;
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
         static_cast<uint64_t>((x.v - residue(f, x.u)) / f.index);
}

Cell cell_of(const GlobeFrame& f, uint64_t idx) {
  const hexify::SolidTopology& t = hexify::topo();
  if (idx == 0) return {0, 0, 0};
  if (idx == f.n_cells - 1) return {t.south_pole(), 0, 0};
  const uint64_t k = idx - 1;
  const uint64_t within = k % f.per_quad;
  const long long per_row = f.dim / f.index;
  Cell x;
  x.quad = static_cast<int>(k / f.per_quad) + 1;
  x.u = static_cast<long long>(within / per_row);
  x.v = static_cast<long long>(within % per_row) * f.index + residue(f, x.u);
  return x;
}

// The slot of lattice point (u, v) of diamond quad q in level L's part;
// false off the padded block.
bool slot_of(const Level& L, int q, long long u, long long v, uint64_t& key) {
  const long long row = u + L.pr;
  if (row < 0 || row >= L.hp) return false;
  const long long r = v - residue(L.f, u);
  if (r % L.f.index != 0) return false;
  const long long col = r / L.f.index + L.pc;
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

// Every slot of level L that holds cell x: its own place in its quad's box,
// and each padding point the quad edge maps (canonicalize_q2d(), as the
// shader's canonicalize() reads them) carry into it. A point past one edge
// maps across it, or to the vertex quad where a far edge starts when
// along = 0; a point past both far edges is the far corner's cell; past any
// other two edges no cell owns it. The edge maps are lattice isometries, so
// each preimage is solved for exactly.
template <typename Emit>
void cell_slots(const Level& L, const Cell& x, Emit emit) {
  const hexify::SolidTopology& t = hexify::topo();
  const long long top = L.f.dim;
  const bool vertex_quad = t.is_pole(x.quad);
  // The padding reaches this far past an edge along either axis, so only a
  // cell this near an edge of its box is carried into it.
  const long long strip = std::max(L.pr, L.pc * L.f.index + L.f.index);
  uint64_t key;
  if (!vertex_quad) {
    if (slot_of(L, x.quad, x.u, x.v, key)) emit(key);
    if (std::min(std::min(x.u, x.v), std::min(top - 1 - x.u, top - 1 - x.v)) > strip) return;
  }
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

// The word a cell is stored as. At the grid's resolution, its value as a
// 32-bit float, NA as one NaN. At a coarser one, where the shader reads the
// fill only, its place on the colour ramp in the high 16 bits (0xffff for
// NA) and its cover in the low 16, both in steps of 1 / 65534, so that no
// word is ABSENT.
uint32_t value_word(double v) {
  if (ISNAN(v)) return kNA;
  const float f = static_cast<float>(v);
  uint32_t w;
  std::memcpy(&w, &f, 4);
  return w;
}

uint32_t coarse_word(double v, double cover, const double ramp[3]) {
  uint32_t pos = 0xffffU;
  if (!ISNAN(v)) {
    const double x = std::min(std::max((v - ramp[0]) * ramp[1] + ramp[2], 0.0), 1.0);
    pos = static_cast<uint32_t>(std::lround(x * 65534.0));
  }
  const uint32_t c = static_cast<uint32_t>(std::max(1L, std::lround(cover * 65534.0)));
  return (pos << 16) | c;
}

// A cell given at a level: its index, the share of its area covered by the
// cells given at the grid's resolution, and its value (NaN for NA).
struct Present {
  uint64_t idx;
  double cover;
  double value;
};

struct Acc {
  double cover = 0, weight = 0, sum = 0;
};

// A coarse cell is kept while the cells given cover at least this share of
// it; below, its fill would be at most this opaque, and a lone cell would
// cost an entry at every coarser level.
constexpr double kLeastCover = 1.0 / 4;

// The cells of the coarser level `coarse` from those of `fine`. A coarse
// cell takes the fine cells whose centres lie in it, a centre on its edge or
// corner split evenly among the cells meeting there (nearest_eisenstein_
// points()); its cover is the share of it they cover, and its value the
// mean of theirs, weighted by the area they cover. A step of aperture a
// gives a cell with n neighbours 1 + n (a - 1) / 6 children's worth: a for a
// hexagon, less at a vertex of the solid. A vertex quad's cell takes the cell
// at its vertex alone besides its neighbours.
std::vector<Present> coarser(const GlobeFrame& fine, const GlobeFrame& coarse,
                             const std::vector<Present>& cells) {
  const hexify::SolidTopology& t = hexify::topo();
  // Shares are summed in a table over all coarse cells when the cells given
  // are about as many, else gathered as (cell, share) pairs and summed in
  // order of cell.
  const bool dense = coarse.n_cells <= 4 * static_cast<uint64_t>(cells.size()) + 1024;
  std::vector<Acc> table(dense ? coarse.n_cells : 0);
  struct Share {
    uint64_t idx;
    Acc acc;
  };
  std::vector<Share> shares;
  if (!dense) shares.reserve(2 * cells.size());

  const long long k = fine.dim / coarse.dim;
  const long long a = coarse.ga, b = coarse.gb;
  const long long D = k * hexify::eisenstein_norm(a - b, b);
  for (const Present& p : cells) {
    const Cell x = cell_of(fine, p.idx);
    auto add = [&](uint64_t idx, double share) {
      Acc one;
      one.cover = share * p.cover;
      if (!ISNAN(p.value)) {
        one.weight = share * p.cover;
        one.sum = share * p.cover * p.value;
      }
      if (dense) {
        Acc& s = table[idx];
        s.cover += one.cover;
        s.weight += one.weight;
        s.sum += one.sum;
      } else {
        shares.push_back({idx, one});
      }
    };
    if (t.is_pole(x.quad)) {
      add(x.quad == 0 ? 0 : coarse.n_cells - 1, 1.0);
      continue;
    }
    // (u + v*omega) * ((a - b) - b*omega), using omega^2 = -1 - omega
    const long long U = x.u * (a - b) + x.v * b;
    const long long V = -x.u * b + x.v * a;
    long long nearest[3][2];
    const int n = hexify::nearest_eisenstein_points(U, V, D, nearest);
    for (int m = 0; m < n; m++) {
      Cell y;
      y.quad = x.quad;
      y.u = nearest[m][0] * a - nearest[m][1] * b;
      y.v = nearest[m][0] * b + nearest[m][1] * a - nearest[m][1] * b;
      if (!hexify::substrate_ij_canonicalize(y.quad, y.u, y.v, coarse.dim)) continue;
      add(cell_index(coarse, y), 1.0 / n);
    }
  }

  const double step = static_cast<double>(fine.per_quad / coarse.per_quad);
  std::vector<Present> out;
  auto keep = [&](uint64_t idx, const Acc& s) {
    const Cell y = cell_of(coarse, idx);
    const bool vertex = t.is_pole(y.quad) || (y.u == 0 && y.v == 0);
    const double whole = vertex ? 1.0 + t.valence[y.quad] * (step - 1.0) / 6.0 : step;
    const double cover = s.cover / whole;
    if (cover < kLeastCover) return;
    out.push_back({idx, std::min(cover, 1.0), s.weight > 0 ? s.sum / s.weight : NAN});
  };
  if (dense) {
    for (uint64_t idx = 0; idx < coarse.n_cells; idx++) {
      if (table[idx].cover > 0) keep(idx, table[idx]);
    }
  } else {
    std::sort(shares.begin(), shares.end(),
              [](const Share& p, const Share& q) { return p.idx < q.idx; });
    for (size_t k = 0; k < shares.size();) {
      Acc s;
      size_t m = k;
      for (; m < shares.size() && shares[m].idx == shares[k].idx; m++) {
        s.cover += shares[m].acc.cover;
        s.weight += shares[m].acc.weight;
        s.sum += shares[m].acc.sum;
      }
      keep(shares[k].idx, s);
      k = m;
    }
  }
  return out;
}

uint32_t mix32(uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

// The offset a key reads; globe.wgsl's bucket() is the same
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
  std::mt19937_64 rng(seed);
  offset.assign(n_buckets, 0);
  slot.assign(E, 0);
  uint32_t cursor = 0;
  std::vector<uint32_t> base;
  for (uint32_t b : order) {
    const uint32_t first = start[b], size = start[b + 1] - first;
    if (size == 0) break;
    base.resize(size);
    for (uint32_t m = 0; m < size; m++) {
      base[m] = static_cast<uint32_t>(e[member[first + m]].key) % M;
    }
    uint32_t o = 0;
    if (size == 1) {
      while (used[cursor]) cursor++;
      o = (cursor + M - base[0]) % M;
    } else {
      std::vector<uint32_t> sorted(base);
      std::sort(sorted.begin(), sorted.end());
      if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) return false;
      bool ok = false;
      for (int tries = 0; tries < (1 << 22) && !ok; tries++) {
        const uint32_t s0 = static_cast<uint32_t>(rng() % M);
        if (used[s0]) continue;
        o = (s0 + M - base[0]) % M;
        ok = true;
        for (uint32_t m = 1; m < size && ok; m++) ok = !used[(base[m] + o) % M];
      }
      if (!ok) return false;
    }
    offset[b] = o;
    for (uint32_t m = 0; m < size; m++) {
      const uint32_t s = (base[m] + o) % M;
      used[s] = 1;
      slot[member[first + m]] = s;
    }
  }
  return true;
}

// The side lengths of a texture array holding n words, `width` wide when
// it fits, its rows spread evenly over the fewest layers
void texture_shape(uint64_t n, uint64_t width, uint32_t& w, uint32_t& h, uint32_t& layers) {
  w = static_cast<uint32_t>(std::max<uint64_t>(1, std::min(width, kSide)));
  const uint64_t rows = std::max<uint64_t>(1, (n + w - 1) / w);
  const uint64_t l = (rows + kSide - 1) / kSide;
  if (l > kLayers) {
    Rcpp::stop("hex_globe() holds at most %.0f values per grid on a graphics card",
               static_cast<double>(kSide * kSide * kLayers));
  }
  layers = static_cast<uint32_t>(l);
  h = static_cast<uint32_t>((rows + l - 1) / l);
}

} // namespace

// The table of a grid's cells for the globe shader. `levels` are the grid's
// resolutions from its own down, each list(resolution, aperture, ap_seq) as
// the C++ entry points take them. With `all_cells`, `values` holds one value
// per cell of the grid in ID order; otherwise the cells are `cell_id`, with
// `values` one per cell or empty when the cells carry none. `ramp` places a
// value v on the colour ramp at clamp((v - ramp[0]) * ramp[1] + ramp[2], 0, 1),
// as ramp_map() in R, for the coarser levels' words.
//
// Returns the levels' frames and blocks (one row each) and the table:
// `values`, the word of each slot (base64 of little-endian 32-bit words),
// laid out as a texture array of `size` = width, height, layers; with
// `keyed`, the perfect hash's `keys` (two words per slot, high first) and
// `offsets` (a texture of `offsets_size`), read with `seed` and `m` slots.
// [[Rcpp::export]]
List cpp_globe_table(NumericVector icosa, List levels, NumericVector cell_id,
                     NumericVector values, bool all_cells, NumericVector ramp) {
  activate_grid(icosa);
  const int n_levels = levels.size();
  if (n_levels == 0) stop("cpp_globe_table needs a level");
  if (n_levels > 1 && ramp.size() != 3) stop("coarser levels need the colour ramp's map");
  const double ramp_at[3] = {ramp.size() == 3 ? ramp[0] : 0.0, ramp.size() == 3 ? ramp[1] : 0.0,
                             ramp.size() == 3 ? ramp[2] : 0.0};

  std::vector<Level> lv;
  uint64_t total = 0;
  for (int k = 0; k < n_levels; k++) {
    List spec = levels[k];
    IntegerVector seq = spec["ap_seq"];
    const GlobeFrame f = hexify::globe_frame(as<int>(spec["resolution"]),
                                             as<int>(spec["aperture"]),
                                             std::vector<int>(seq.begin(), seq.end()));
    if (k > 0 && (lv.back().f.dim % f.dim != 0 || lv.back().f.per_quad % f.per_quad != 0)) {
      stop("each level of the globe table must be the parent grid of the one before");
    }
    lv.push_back(make_level(f, total));
    total += level_size(lv.back());
  }

  // The cells given at each level
  std::vector<std::vector<Present>> given(n_levels);
  const GlobeFrame& finest = lv[0].f;
  if (all_cells) {
    if (static_cast<uint64_t>(values.size()) != finest.n_cells) {
      stop("values must hold one value per cell of the grid");
    }
    given[0].resize(finest.n_cells);
    for (uint64_t idx = 0; idx < finest.n_cells; idx++) given[0][idx] = {idx, 1.0, values[idx]};
  } else {
    hexify::require_cell_ids(cell_id);
    const bool valued = values.size() > 0;
    if (valued && values.size() != cell_id.size()) stop("values must hold one value per cell");
    given[0].reserve(cell_id.size());
    for (R_xlen_t k = 0; k < cell_id.size(); k++) {
      const int64_t id = hexify::cell_id_get(cell_id[k]);
      if (id == hexify::kCellIdNA || id < 1 || static_cast<uint64_t>(id) > finest.n_cells) {
        stop("cell IDs must name cells of the grid");
      }
      given[0].push_back({static_cast<uint64_t>(id - 1), 1.0, valued ? values[k] : 0.0});
    }
    std::sort(given[0].begin(), given[0].end(),
              [](const Present& p, const Present& q) { return p.idx < q.idx; });
  }
  for (int k = 1; k < n_levels; k++) given[k] = coarser(lv[k - 1].f, lv[k].f, given[k - 1]);

  uint64_t n_given = 0;
  for (const auto& g : given) n_given += g.size();
  const bool keyed = total >= (uint64_t(1) << 31) || total > 3 * n_given + 4096;

  // Every slot of every given cell, through `emit(key, word)`
  auto each_slot = [&](auto emit) {
    for (int k = 0; k < n_levels; k++) {
      for (const Present& p : given[k]) {
        const uint32_t word = k == 0 ? value_word(p.value) : coarse_word(p.value, p.cover, ramp_at);
        cell_slots(lv[k], cell_of(lv[k].f, p.idx), [&](uint64_t key) { emit(key, word); });
      }
    }
  };

  uint32_t w, h, layers, seed = 0, n_buckets = 1, ow = 1, oh = 1, o_layers = 1;
  uint64_t m;
  std::vector<uint32_t> words, keys, offsets(1, 0);
  if (!keyed) {
    m = total;
    texture_shape(m, static_cast<uint64_t>(lv[0].wp), w, h, layers);
    words.assign(static_cast<size_t>(w) * h * layers, kAbsent);
    each_slot([&](uint64_t key, uint32_t word) {
      if (words[key] != kAbsent) stop("hexify internal error: two cells share a globe table slot");
      words[key] = word;
    });
  } else {
    std::vector<Entry> entries;
    each_slot([&](uint64_t key, uint32_t word) { entries.push_back({key, word}); });
    std::sort(entries.begin(), entries.end(),
              [](const Entry& p, const Entry& q) { return p.key < q.key; });
    for (size_t k = 1; k < entries.size(); k++) {
      if (entries[k].key == entries[k - 1].key) {
        stop("hexify internal error: two cells share a globe table slot");
      }
    }
    const uint64_t E = entries.size();
    m = std::max<uint64_t>(1, E + E / 32 + 16);
    if (m >= (uint64_t(1) << 31)) stop("hex_globe() holds at most 2^31 cells on a graphics card");
    texture_shape(m, kSide, w, h, layers);
    n_buckets = static_cast<uint32_t>(std::max<uint64_t>(1, E / 3));
    std::vector<uint32_t> slot;
    bool placed = false;
    for (int attempt = 0; attempt < 40 && !placed; attempt++) {
      seed = mix32(0x9e3779b9U + static_cast<uint32_t>(attempt));
      placed = place(entries, static_cast<uint32_t>(m), n_buckets, seed, offsets, slot);
      if (!placed) n_buckets += n_buckets / 8 + 1;
    }
    if (!placed) stop("hexify internal error: the globe table's cells could not be hashed");
    words.assign(static_cast<size_t>(w) * h * layers, kAbsent);
    keys.assign(2 * words.size(), kAbsent);
    for (size_t k = 0; k < E; k++) {
      words[slot[k]] = entries[k].word;
      keys[2 * slot[k]] = static_cast<uint32_t>(entries[k].key >> 32);
      keys[2 * slot[k] + 1] = static_cast<uint32_t>(entries[k].key);
    }
    texture_shape(n_buckets, kSide, ow, oh, o_layers);
    if (o_layers > 1) stop("hex_globe() holds at most %.0f cells on a graphics card",
                           static_cast<double>(3 * kSide * kSide));
    offsets.resize(static_cast<size_t>(ow) * oh, 0);
  }

  NumericMatrix frames(n_levels, 11);
  for (int k = 0; k < n_levels; k++) {
    const Level& L = lv[k];
    const double row[11] = {static_cast<double>(L.f.dim), static_cast<double>(L.f.index),
                            static_cast<double>(L.f.c), static_cast<double>(L.f.ga),
                            static_cast<double>(L.f.gb), static_cast<double>(L.pr),
                            static_cast<double>(L.pc), static_cast<double>(L.hp),
                            static_cast<double>(L.wp), static_cast<double>(L.base),
                            static_cast<double>(L.f.n_cells)};
    for (int c = 0; c < 11; c++) frames(k, c) = row[c];
  }
  colnames(frames) = CharacterVector::create("dim", "index", "c", "ga", "gb", "pr", "pc",
                                             "hp", "wp", "base", "n_cells");
  IntegerVector given_n(n_levels);
  for (int k = 0; k < n_levels; k++) given_n[k] = static_cast<int>(given[k].size());

  return List::create(
      _["levels"] = frames,
      _["given"] = given_n,
      _["keyed"] = keyed,
      _["m"] = static_cast<double>(m),
      _["size"] = NumericVector::create(w, h, layers),
      _["values"] = hexify::base64_words(words),
      _["keys"] = keyed ? hexify::base64_words(keys) : std::string(),
      _["seed"] = static_cast<double>(seed),
      _["buckets"] = static_cast<double>(n_buckets),
      _["offsets_size"] = NumericVector::create(ow, oh),
      _["offsets"] = hexify::base64_words(offsets));
}
