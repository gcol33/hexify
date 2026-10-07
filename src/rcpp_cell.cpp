// rcpp_cell.cpp
// Rcpp bindings for cell ID conversion and coordinate transforms
//
// This file provides the R interface for:
// - Lon/lat to cell ID conversion
// - Cell ID to lon/lat conversion
// - Cell ID to cell info conversion
// - Quad IJ coordinate conversion
// - PLANE coordinate conversions
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include "constants.h"
#include "polyhedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "grid_math.h"
#include "coordinate_transforms.h"
#include "rcpp_icosa.h"
#include "cell_walls.h"

using namespace Rcpp;

// ============================================================================
// Triangle to Quad Coordinate Conversion
// ============================================================================

// [[Rcpp::export]]
Rcpp::List cpp_icosa_tri_to_quad_ij(NumericVector icosa, int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                                     int aperture, int resolution) {
    activate_grid(icosa);
    int quad;
    long long i, j;

    hexify::icosa_tri_to_quad_ij(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y, aperture, resolution, quad, i, j);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("i") = (double)i,
        Rcpp::Named("j") = (double)j
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_icosa_tri_to_quad_xy(NumericVector icosa, int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y) {
    activate_grid(icosa);
    int quad;
    double quad_x, quad_y;

    hexify::icosa_tri_to_quad_xy(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y, quad, quad_x, quad_y);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("quad_x") = quad_x,
        Rcpp::Named("quad_y") = quad_y
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_quad_xy_to_icosa_tri(NumericVector icosa, int quad, double quad_x, double quad_y) {
    activate_grid(icosa);
    int icosa_triangle_face;
    double icosa_triangle_x, icosa_triangle_y;

    hexify::quad_xy_to_icosa_tri(quad, quad_x, quad_y, icosa_triangle_face, icosa_triangle_x, icosa_triangle_y);

    return Rcpp::List::create(
        Rcpp::Named("icosa_triangle_face") = icosa_triangle_face,
        Rcpp::Named("icosa_triangle_x") = icosa_triangle_x,
        Rcpp::Named("icosa_triangle_y") = icosa_triangle_y
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_quad_ij_to_xy(int quad, double i, double j,
                              int aperture, int resolution) {
    double quad_x, quad_y;
    hexify::quad_ij_to_xy(quad, static_cast<long long>(i), static_cast<long long>(j),
                          aperture, resolution, quad_x, quad_y);

    return Rcpp::List::create(
        Rcpp::Named("quad_x") = quad_x,
        Rcpp::Named("quad_y") = quad_y
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_lonlat_to_quad_ij(NumericVector icosa, double lon_deg, double lat_deg,
                                  int aperture, int resolution) {
    activate_grid(icosa);
    // Step 1: Forward project to icosa triangle coordinates
    hexify::ProjectionResult fwd = hexify::snyder_forward(lon_deg, lat_deg);

    // Step 2: Convert icosa triangle coords to quad integer coords
    int quad;
    long long i, j;
    hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y, aperture, resolution, quad, i, j);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("i") = (double)i,
        Rcpp::Named("j") = (double)j,
        Rcpp::Named("icosa_triangle_face") = fwd.face,
        Rcpp::Named("icosa_triangle_x") = fwd.icosa_triangle_x,
        Rcpp::Named("icosa_triangle_y") = fwd.icosa_triangle_y
    );
}

// ============================================================================
// Substrate sublattice packing
// ============================================================================
// Cells sit on a sublattice of index N in the substrate. Writing a substrate
// point as i + j*omega with omega = exp(2*pi*i/3), the cells are the multiples
// of a generator of norm N, which is the single congruence
//
//   j = c * i   (mod N)
//
// N is 1 (the grid is the substrate), 3 (the 30-degree Class II lattice), 7 or
// 21 (the lattices an odd number of aperture-7 levels leaves). A quad's
// dimension is divisible by N, so each quad holds exactly dim * dim / N cells
// and the numbering below is dense with no gaps.
struct SubstrateLattice {
    long long index;  // N
    long long c;      // j = c * i (mod N)
};

static const SubstrateLattice kAlignedLattice = {1, 0};

// Residue the cells of column i occupy
static inline long long lattice_residue(long long i, const SubstrateLattice& lat) {
    return ((lat.c * i) % lat.index + lat.index) % lat.index;
}

// 2D cell index within a quad
static uint64_t cell_index_2d(long long i, long long j, long long dim,
                              const SubstrateLattice& lat) {
    if (lat.index == 1) {
        return static_cast<uint64_t>(i) * dim + j;
    }
    return static_cast<uint64_t>(i) * (dim / lat.index) +
           (j - lattice_residue(i, lat)) / lat.index;
}

// Inverse: 2D cell index back to (i, j)
static void ij_from_cell_index(uint64_t idx, long long dim,
                               const SubstrateLattice& lat,
                               long long& i, long long& j) {
    if (lat.index == 1) {
        i = static_cast<long long>(idx / dim);
        j = static_cast<long long>(idx % dim);
        return;
    }
    long long per_column = dim / lat.index;
    i = static_cast<long long>(idx / per_column);
    j = static_cast<long long>(idx % per_column) * lat.index + lattice_residue(i, lat);
}

// Inverse of a mod N, for the N in {3, 7, 21} that occur here
static long long lattice_mod_inverse(long long a, long long N) {
    a = (a % N + N) % N;
    for (long long k = 1; k < N; ++k) {
        if ((a * k) % N == 1) return k;
    }
    Rcpp::stop("substrate lattice: " + std::to_string(a) +
               " has no inverse modulo " + std::to_string(N));
}

// Sublattice of a grid form. Its generator m + n*w (w = exp(pi*i/3)) is
// (m + n) + n*omega, and a substrate point is a multiple of it exactly when
// j = n / (m + n) * i (mod N). Both m + n and n are invertible mod N: a factor
// shared with N would divide the other as well and so square-divide N, which
// is 3, 7 or 21.
static SubstrateLattice sublattice_of(const hexify::HexGridForm& form) {
    long long N = hexify::eisenstein_norm(form.m, form.n);
    if (N == 1) return kAlignedLattice;
    if (N != 3 && N != 7 && N != 21) {
        Rcpp::stop("substrate lattice: unexpected lattice norm " + std::to_string(N));
    }
    long long b = (form.n % N + N) % N;
    long long inv_a = lattice_mod_inverse(form.m + form.n, N);
    SubstrateLattice lat;
    lat.index = N;
    lat.c = (b * inv_a) % N;
    return lat;
}

// ============================================================================
// Aperture 7: Surrogate-based encoding
// ============================================================================
// Surrogates are the canonical cell coordinates for aperture 7.
// Each ap7 cell corresponds to exactly one surrogate (i,j), unlike substrates
// where the substrate grid has ~7x more positions than actual cells.
//
// The 19.1 degree rotation between the surrogate and quad frames puts the
// surrogates of one quad in a skewed patch rather than an axis-aligned box, so
// the within-quad index comes from ap7_surrogate_to_quad_index(), which walks
// the cell centres in the quad's substrate box. It spans [0, 7^res) exactly.
// ============================================================================

// The same sublattice written as its generator a + b*omega, in the (i, j)
// coordinates cpp_cell_to_quad_ij() returns. The cells of a quad are its
// multiples over the Eisenstein integers, so multiplying by the six units
// gives the six neighbours, and dividing a difference of two cells by it
// reads that difference in cell steps.
struct LatticeGenerator {
    long long a;
    long long b;
};

// The six cells one step away, as offsets in (i, j): the generator times each
// unit of the Eisenstein integers, using omega^2 = -1 - omega.
static void lattice_unit_steps(const LatticeGenerator& g, long long steps[6][2]) {
    static const long long units[6][2] = {
        { 1,  0}, { 1,  1}, { 0,  1},
        {-1,  0}, {-1, -1}, { 0, -1}
    };
    for (int k = 0; k < 6; k++) {
        long long p = units[k][0], q = units[k][1];
        steps[k][0] = p * g.a - q * g.b;
        steps[k][1] = p * g.b + q * g.a - q * g.b;
    }
}

// ============================================================================
// Quad frame of a grid
// ============================================================================
// How a grid lays its cells out in a quad: how many cells a quad holds, how
// many substrate steps a quad edge measures, which coordinates a cell is
// stored by, and how a stored coordinate maps to and from the quad plane. A
// pure aperture and a mixed aperture sequence differ only in these, so every
// conversion below reads one frame.
//
// Apertures 3 and 4 and every mixed sequence store a cell by the substrate
// coordinates of its centre, which lie on the sublattice of the grid's form.
// Aperture 7 stores its surrogate, the cell's own Class I coordinate, so its
// stored coordinates pack and step as the aligned lattice.
// ============================================================================

struct QuadFrame {
    std::vector<int> ap_seq;      // empty for a pure aperture
    int aperture;                 // 0 for a mixed sequence
    int resolution;
    hexify::HexGridForm form;
    uint64_t nCells;              // the diamond quads of cells plus the two vertex quads
    uint64_t offsetPerQuad;       // cells per quad, the product of the apertures
    long long dim;                // substrate steps along a quad edge
    SubstrateLattice lattice;     // which stored (i, j) are cells
    LatticeGenerator generator;   // that lattice's generator, in stored (i, j)
};

// The sublattice generator m + n*w (w = exp(pi*i/3)) of a form, written in
// the substrate's (i, j) as (m + n) + n*omega.
static LatticeGenerator generator_of(const hexify::HexGridForm& form) {
    return {form.m + form.n, form.n};
}

// The frame of a pure aperture at a resolution (empty 'ap_seq'), or of a mixed
// sequence, whose resolution is one less than its length.
static QuadFrame quad_frame(int resolution, int aperture, std::vector<int> ap_seq) {
    const bool mixed = !ap_seq.empty();
    if (mixed) {
        resolution = static_cast<int>(ap_seq.size()) - 1;
        aperture = 0;
    } else if (aperture != 3 && aperture != 4 && aperture != 7) {
        Rcpp::stop("aperture must be 3, 4, or 7");
    }
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }

    QuadFrame f;
    f.aperture = aperture;
    f.resolution = resolution;
    f.form = mixed ? hexify::hex_form_sequence(ap_seq)
                   : hexify::hex_form_pure(aperture, resolution);
    f.dim = mixed ? hexify::quad_edge_dim(ap_seq)
                  : hexify::quad_edge_dim(aperture, resolution);
    f.offsetPerQuad = 1;
    for (int k = 1; k <= resolution; k++) {
        f.offsetPerQuad *= static_cast<uint64_t>(mixed ? ap_seq[k] : aperture);
    }
    f.nCells = static_cast<uint64_t>(hexify::topo().n_diamonds()) * f.offsetPerQuad + 2;
    if (aperture == 7) {
        f.lattice = kAlignedLattice;
        f.generator = {1, 0};
    } else {
        f.lattice = sublattice_of(f.form);
        f.generator = generator_of(f.form);
    }
    f.ap_seq = std::move(ap_seq);
    return f;
}

// The frame an entry point's (resolution, aperture, ap_seq) names. An empty
// ap_seq names the pure aperture; a mixed sequence comes with aperture 0 and
// the resolution its length gives.
static QuadFrame grid_frame(int resolution, int aperture, const IntegerVector& ap_seq) {
    if (ap_seq.size() == 0) {
        return quad_frame(resolution, aperture, {});
    }
    if (aperture != 0 || resolution != ap_seq.size() - 1) {
        Rcpp::stop("a mixed aperture sequence takes aperture 0 and resolution "
                   "length(ap_seq) - 1");
    }
    return quad_frame(resolution, 0, std::vector<int>(ap_seq.begin(), ap_seq.end()));
}

// The ID, from 1, of the cell at stored (i, j) of a quad. Quad 0 holds the
// north pole alone, ID 1; quad q > 0 follows it and the q - 1 quads before.
static inline double frame_encode(const QuadFrame& f, int quad, long long i, long long j) {
    uint64_t offset = (quad == 0) ? 0 : 1 + static_cast<uint64_t>(quad - 1) * f.offsetPerQuad;
    uint64_t within_quad = (f.aperture == 7)
        ? hexify::ap7_surrogate_to_quad_index(i, j, f.resolution)
        : cell_index_2d(i, j, f.dim, f.lattice);
    return static_cast<double>(offset + within_quad + 1);
}

// The 0-based index of a cell ID. Stops unless the ID names a cell of the
// grid, so no table is ever indexed with a quad, i or j read from NA or an
// out-of-range ID.
static inline uint64_t frame_cell_index(const QuadFrame& f, double cell_id_raw) {
    if (!std::isfinite(cell_id_raw) || cell_id_raw < 1.0 ||
        cell_id_raw > static_cast<double>(f.nCells)) {
        Rcpp::stop("cell_id must be a finite value in [1, %.0f] for resolution %d",
                   static_cast<double>(f.nCells), f.resolution);
    }
    return static_cast<uint64_t>(cell_id_raw) - 1;
}

// A cell's quad and stored (i, j) from its 0-based index
static inline void frame_decode_index(const QuadFrame& f, uint64_t idx,
                                      int& quad, long long& i, long long& j) {
    if (idx == 0) {
        quad = 0;
        i = 0;
        j = 0;
        return;
    }
    idx--;
    quad = static_cast<int>(idx / f.offsetPerQuad) + 1;
    uint64_t within_quad = idx - static_cast<uint64_t>(quad - 1) * f.offsetPerQuad;
    if (f.aperture == 7) {
        hexify::ap7_quad_index_to_surrogate(within_quad, f.resolution, i, j);
    } else {
        ij_from_cell_index(within_quad, f.dim, f.lattice, i, j);
    }
}

// A cell's quad and stored (i, j) from its ID
static inline void frame_decode(const QuadFrame& f, double cell_id_raw,
                                int& quad, long long& i, long long& j) {
    frame_decode_index(f, frame_cell_index(f, cell_id_raw), quad, i, j);
}

static inline bool frame_in_quad(const QuadFrame& f, long long i, long long j) {
    if (f.aperture == 7) {
        return hexify::ap7_surrogate_in_quad(i, j, f.resolution);
    }
    return i >= 0 && j >= 0 && i < f.dim && j < f.dim;
}

// The quad-plane centre of the cell at stored (i, j)
static inline void frame_ij_to_xy(const QuadFrame& f, int quad, long long i, long long j,
                                  double& x, double& y) {
    if (f.ap_seq.empty()) {
        hexify::quad_ij_to_xy(quad, i, j, f.aperture, f.resolution, x, y);
    } else {
        hexify::center_form(f.form, i, j, x, y);
    }
}

// Re-express a stored coordinate that has stepped outside its quad in the quad
// that owns it. False where no quad owns it, at the solid's fold around
// a vertex.
static inline bool frame_canonicalize(const QuadFrame& f, int& quad,
                                      long long& i, long long& j) {
    if (f.ap_seq.empty()) {
        return hexify::quad_ij_canonicalize(quad, i, j, f.aperture, f.resolution);
    }
    return hexify::substrate_ij_canonicalize(quad, i, j, f.dim);
}

// Floor of a / b for b > 0
static inline long long floor_div(long long a, long long b) {
    long long q = a / b;
    return (a % b != 0 && a < 0) ? q - 1 : q;
}

// The parent, one level up, of the cell at stored (i, j) of a quad of a mixed
// sequence: the parent-lattice point nearest the cell's centre, in that quad's
// plane, re-expressed in the quad that owns it.
//
// Both grids store a cell by the substrate coordinates i + j*omega of its
// centre on the same quad axes, and the child's substrate is k = dim_c / dim_p
// (1, 2, 3 or 7) times finer, so the centre reads c / k in the parent's
// substrate. Parent cells are the multiples of the parent generator g, so in
// parent-cell units the centre is u = c * conj(g) / (k * N(g)), a point of
// the Eisenstein lattice divided by D = k * N(g). The nearest lattice point
// is a corner of the rhombus [x0, x0 + 1] x [y0, y0 + 1] containing u (its
// two halves are equilateral triangles), and the squared distances
// N(D * corner - D * u) are integers, so the comparison is exact.
//
// An aperture-3 step puts child centres on parent corners (three nearest
// points) and an aperture-4 step on parent edge midpoints (two); an
// aperture-7 step has no ties. A tie goes to the candidate with the larger
// 2x + y, which differs between any two candidates a unit apart. The rule
// is a translation in parent-cell units, so every hexagonal parent inside a
// quad receives exactly as many children as the step's aperture.
static void frame_parent(const QuadFrame& child, const QuadFrame& parent,
                         int& quad, long long& i, long long& j) {
    const long long k = child.dim / parent.dim;
    const long long a = parent.generator.a, b = parent.generator.b;
    const long long D = k * hexify::eisenstein_norm(a - b, b);

    // (i + j*omega) * ((a - b) - b*omega), using omega^2 = -1 - omega
    const long long c = a - b, d = -b;
    const long long U = i * c - j * d;
    const long long V = i * d + j * c - j * d;
    const long long x0 = floor_div(U, D), y0 = floor_div(V, D);

    long long best_x = 0, best_y = 0, best_dist = -1, best_key = 0;
    for (int corner = 0; corner < 4; corner++) {
        const long long x = x0 + (corner & 1), y = y0 + (corner >> 1);
        const long long dx = D * x - U, dy = D * y - V;
        const long long dist = hexify::eisenstein_norm(dx - dy, dy);
        const long long key = 2 * x + y;
        if (best_dist < 0 || dist < best_dist ||
            (dist == best_dist && key > best_key)) {
            best_x = x;
            best_y = y;
            best_dist = dist;
            best_key = key;
        }
    }

    i = best_x * a - best_y * b;
    j = best_x * b + best_y * a - best_y * b;
    if (!frame_in_quad(parent, i, j) && !frame_canonicalize(parent, quad, i, j)) {
        Rcpp::stop("hexify internal error: a parent cell centre has no owning quad");
    }
}

// The cell a lon/lat point falls in: its quad and stored (i, j)
static inline void frame_locate(const QuadFrame& f, double lon_deg, double lat_deg,
                                int& quad, long long& i, long long& j) {
    hexify::ProjectionResult fwd = hexify::snyder_forward(lon_deg, lat_deg);

    if (f.ap_seq.empty()) {
        hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x,
                                     fwd.icosa_triangle_y, f.aperture,
                                     f.resolution, quad, i, j);
        return;
    }

    int quad_pre;
    double quad_x, quad_y;
    hexify::icosa_tri_to_quad_xy(fwd.face, fwd.icosa_triangle_x,
                                 fwd.icosa_triangle_y, quad_pre, quad_x, quad_y);
    hexify::quad_xy_to_ij_mixed(quad_pre, quad_x, quad_y, f.form, f.dim,
                                quad, i, j);
}

// Substrate steps along a quad edge of a grid
// [[Rcpp::export]]
double cpp_quad_edge_dim(int resolution, int aperture, IntegerVector ap_seq) {
    return static_cast<double>(grid_frame(resolution, aperture, ap_seq).dim);
}

// The generator a + b*omega of the lattice the stored (i, j) of a grid's
// cells occupy, in the coordinates cpp_cell_to_quad_ij() returns
// [[Rcpp::export]]
NumericVector cpp_cell_lattice_generator(int resolution, int aperture,
                                         IntegerVector ap_seq) {
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    return NumericVector::create(static_cast<double>(f.generator.a),
                                 static_cast<double>(f.generator.b));
}

// [[Rcpp::export]]
NumericVector cpp_quad_ij_to_cell(NumericVector icosa, IntegerVector quad, NumericVector i,
                                  NumericVector j, int resolution, int aperture,
                                  IntegerVector ap_seq) {
    activate_grid(icosa);
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    R_xlen_t n = quad.size();
    NumericVector result(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        // A coordinate that has stepped outside its quad names a cell of a
        // neighbouring quad, so re-express it there before packing. The walk
        // up and down the cell hierarchy reaches these: an ancestor of a cell
        // near a quad edge need not lie in the same quad. Coordinates already
        // inside their quad pass through unchanged.
        if (!frame_canonicalize(f, q, ii, jj)) {
            // Outside every adjacent quad, which is where the solid
            // folds at a vertex. No cell owns the coordinate.
            result[k] = NA_REAL;
            continue;
        }
        result[k] = frame_encode(f, q, ii, jj);
    }

    return result;
}

// The parent of each cell of a mixed sequence `ap_seq` in the grid one level
// coarser, `parent_seq`. A grid depends only on the multiset of its refinement
// steps, so the parent grid is any sequence whose steps are the child's with
// one removed; a family spelling such as "4/3" removes a leading step, not
// the last.
// [[Rcpp::export]]
NumericVector cpp_mixed_parent(NumericVector icosa, NumericVector cell_id,
                               IntegerVector ap_seq, IntegerVector parent_seq) {
    activate_grid(icosa);
    if (ap_seq.size() < 2 || parent_seq.size() != ap_seq.size() - 1) {
        Rcpp::stop("cpp_mixed_parent takes a mixed sequence and the sequence one "
                   "level coarser");
    }
    std::vector<int> steps(ap_seq.begin() + 1, ap_seq.end());
    std::vector<int> parent_steps(parent_seq.begin() + 1, parent_seq.end());
    std::sort(steps.begin(), steps.end());
    std::sort(parent_steps.begin(), parent_steps.end());
    if (!std::includes(steps.begin(), steps.end(),
                       parent_steps.begin(), parent_steps.end())) {
        Rcpp::stop("the parent sequence is not the child sequence with one step removed");
    }

    const QuadFrame child = quad_frame(0, 0, std::vector<int>(ap_seq.begin(), ap_seq.end()));
    const QuadFrame parent = quad_frame(0, 0, std::vector<int>(parent_seq.begin(),
                                                               parent_seq.end()));

    R_xlen_t n = cell_id.size();
    NumericVector result(n);
    for (R_xlen_t k = 0; k < n; k++) {
        if (ISNAN(cell_id[k])) {
            result[k] = NA_REAL;
            continue;
        }
        int quad;
        long long i, j;
        frame_decode(child, cell_id[k], quad, i, j);
        frame_parent(child, parent, quad, i, j);
        result[k] = frame_encode(parent, quad, i, j);
    }
    return result;
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell(NumericVector icosa,
                                 NumericVector lon, NumericVector lat,
                                 int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    R_xlen_t n = lon.size();
    NumericVector result(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_locate(f, lon[k], lat[k], quad, i, j);
        result[k] = frame_encode(f, quad, i, j);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_lonlat(NumericVector icosa, NumericVector cell_id,
                             int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    R_xlen_t n = cell_id.size();
    NumericVector lon(n);
    NumericVector lat(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);

        double quad_x, quad_y;
        frame_ij_to_xy(f, quad, i, j, quad_x, quad_y);

        int icosa_triangle_face;
        double icosa_triangle_x, icosa_triangle_y;
        if (!hexify::try_quad_xy_to_icosa_tri(quad, quad_x, quad_y, icosa_triangle_face,
                                              icosa_triangle_x, icosa_triangle_y)) {
            lon[k] = NA_REAL;
            lat[k] = NA_REAL;
            continue;
        }

        auto ll = hexify::face_xy_to_ll(icosa_triangle_x, icosa_triangle_y, icosa_triangle_face);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }

    return DataFrame::create(
        _["lon_deg"] = lon,
        _["lat_deg"] = lat
    );
}

// Cell IDs to the quad and stored (i, j) each cell is packed from; the inverse
// of cpp_quad_ij_to_cell().
// [[Rcpp::export]]
DataFrame cpp_cell_to_quad_ij(NumericVector icosa, NumericVector cell_id, int resolution, int aperture,
                              IntegerVector ap_seq) {
    activate_grid(icosa);
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    R_xlen_t n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_i(n);
    NumericVector out_j(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);
        out_quad[k] = quad;
        out_i[k] = static_cast<double>(i);
        out_j[k] = static_cast<double>(j);
    }

    return DataFrame::create(
        _["quad"] = out_quad,
        _["i"] = out_i,
        _["j"] = out_j
    );
}

// ============================================================================
// Cell ID to Quad XY Conversion
// ============================================================================
// Converts cell IDs to Quad XY coordinates (continuous).
// Pipeline: Cell ID → Quad IJ → Quad XY
// Produces output compatible with standard ISEA Quad XY representation.
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_cell_to_quad_xy(NumericVector icosa, NumericVector cell_id, int resolution,
                               int aperture) {
    activate_grid(icosa);
    QuadFrame f = quad_frame(resolution, aperture, {});
    R_xlen_t n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_qx(n);
    NumericVector out_qy(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);

        double quad_x, quad_y;
        frame_ij_to_xy(f, quad, i, j, quad_x, quad_y);

        out_quad[k] = quad;
        out_qx[k] = quad_x;
        out_qy[k] = quad_y;
    }

    return DataFrame::create(
        _["quad"] = out_quad,
        _["quad_x"] = out_qx,
        _["quad_y"] = out_qy
    );
}

// ============================================================================
// Quad XY to Cell ID Conversion
// ============================================================================
// Converts Quad XY coordinates (continuous) to cell IDs.
// Pipeline: Quad XY → Quad IJ (quantize) → Cell ID
// Produces cell IDs compatible with standard ISEA numbering.
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_quad_xy_to_cell(NumericVector icosa, IntegerVector quad, NumericVector quad_x,
                                   NumericVector quad_y, int resolution,
                                   int aperture) {
    activate_grid(icosa);
    QuadFrame f = quad_frame(resolution, aperture, {});
    R_xlen_t n = quad.size();
    NumericVector result(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int q = quad[k];
        double qx = quad_x[k];
        double qy = quad_y[k];

        int out_quad;
        long long i, j;
        if (aperture == 7) {
            // AP7: exact-integer quantization straight to the surrogate.
            hexify::quad_xy_to_ij(q, qx, qy, 7, resolution, out_quad, i, j);
        } else {
            // AP3/AP4: through the face the point lies on, which names its quad
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            hexify::quad_xy_to_icosa_tri(q, qx, qy, icosa_triangle_face,
                                         icosa_triangle_x, icosa_triangle_y);
            hexify::icosa_tri_to_quad_ij(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y,
                                         aperture, resolution, out_quad, i, j);
        }
        result[k] = frame_encode(f, out_quad, i, j);
    }

    return result;
}

// ============================================================================
// Cell ID to Icosa Triangle Conversion
// ============================================================================
// Converts cell IDs to icosahedral triangle coordinates (face, x, y).
// Pipeline: Cell ID → Quad IJ → Quad XY → Icosa Triangle
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_cell_to_icosa_tri(NumericVector icosa, NumericVector cell_id, int resolution,
                                 int aperture) {
    activate_grid(icosa);
    QuadFrame f = quad_frame(resolution, aperture, {});
    R_xlen_t n = cell_id.size();
    IntegerVector out_face(n);
    NumericVector out_tx(n);
    NumericVector out_ty(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);

        double quad_x, quad_y;
        frame_ij_to_xy(f, quad, i, j, quad_x, quad_y);

        int icosa_triangle_face;
        double icosa_triangle_x, icosa_triangle_y;
        hexify::quad_xy_to_icosa_tri(quad, quad_x, quad_y, icosa_triangle_face,
                                     icosa_triangle_x, icosa_triangle_y);

        out_face[k] = icosa_triangle_face;
        out_tx[k] = icosa_triangle_x;
        out_ty[k] = icosa_triangle_y;
    }

    return DataFrame::create(
        _["icosa_triangle_face"] = out_face,
        _["icosa_triangle_x"] = out_tx,
        _["icosa_triangle_y"] = out_ty
    );
}

// ============================================================================
// Quad IJ to Icosa Triangle Conversion
// ============================================================================
// Converts Quad IJ coordinates to icosahedral triangle coordinates.
// Pipeline: Quad IJ → Quad XY → Icosa Triangle
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_quad_ij_to_icosa_tri(NumericVector icosa, IntegerVector quad, NumericVector i,
                                    NumericVector j, int resolution,
                                    int aperture) {
    activate_grid(icosa);
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_quad_ij_to_icosa_tri: aperture must be 3, 4, or 7");
    }

    int n = quad.size();
    IntegerVector out_face(n);
    NumericVector out_tx(n);
    NumericVector out_ty(n);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        // Convert Quad IJ → Quad XY
        double quad_x, quad_y;
        hexify::quad_ij_to_xy(q, ii, jj, aperture, resolution, quad_x, quad_y);

        // Convert Quad XY → Icosa Triangle
        int icosa_triangle_face;
        double icosa_triangle_x, icosa_triangle_y;
        hexify::quad_xy_to_icosa_tri(q, quad_x, quad_y, icosa_triangle_face,
                                     icosa_triangle_x, icosa_triangle_y);

        out_face[k] = icosa_triangle_face;
        out_tx[k] = icosa_triangle_x;
        out_ty[k] = icosa_triangle_y;
    }

    return DataFrame::create(
        _["icosa_triangle_face"] = out_face,
        _["icosa_triangle_x"] = out_tx,
        _["icosa_triangle_y"] = out_ty
    );
}

// ============================================================================
// Cell Boundary Generation from Cell ID
// ============================================================================
// A cell's boundary is its centre in quad XY plus six corners one circumradius
// out, turned by the lattice rotation the grid's form carries. Each corner is
// projected on its own, and one that lands outside the quad's valid region is
// carried through the centre's face instead.
// ============================================================================

// Circumradius of a hexagon on the unscaled grid, 1/sqrt(3).
constexpr double kHexCircumradius = 0.57735026918962576451;

// Deepest halving of one cell edge, 2^12 pieces.
constexpr int kMaxEdgeSplits = 12;

// An edge whose corner longitudes are a half turn apart to within this many
// degrees runs over a pole. Corners carry the inverse projection's error, a
// few 1e-5 degrees; a chord this close to a half turn passes within about a
// metre of the pole.
constexpr double kPoleEdgeLonSlack = 1e-3;

// Whether an edge a -> b, whose true midpoint is m, must be split before the
// straight lon/lat chord a -> b may stand for it: m lies further from the
// chord's midpoint than 'tolerance' times the chord's length, both measured
// with longitude scaled by the cosine of the mean latitude.
static inline bool chord_needs_split(double alon, double alat,
                                     double mlon, double mlat,
                                     double blon, double blat,
                                     double tolerance) {
    double dlon = blon - alon;
    if (dlon > 180.0) dlon -= 360.0;
    if (dlon < -180.0) dlon += 360.0;
    double clat = 0.5 * (alat + blat);
    double coslat = std::cos(clat * hexify::kDegToRad);
    double off_lon = mlon - (alon + 0.5 * dlon);
    if (off_lon > 180.0) off_lon -= 360.0;
    if (off_lon < -180.0) off_lon += 360.0;
    double off = std::hypot(off_lon * coslat, mlat - clat);
    double len = std::hypot(dlon * coslat, blat - alat);
    return len > 0.0 && off > tolerance * len;
}

// Deepest halving of one piece for its length alone, 2^30 pieces.
constexpr int kMaxArcSplits = 30;

// How closely a drawn boundary follows the true one: 'tolerance' as for
// chord_needs_split(), and 'max_arc' the longest arc, in radians, one drawn
// piece may span. Zero turns either off.
struct EdgeLimits {
    double tolerance;
    double max_arc;
    bool active() const { return tolerance > 0.0 || max_arc > 0.0; }
};

// The arc between two lon/lat points, in radians.
static inline double lonlat_arc(double alon, double alat,
                                double blon, double blat) {
    const double p1 = alat * hexify::kDegToRad, p2 = blat * hexify::kDegToRad;
    const double s_lat = std::sin(0.5 * (p2 - p1));
    const double s_lon = std::sin(0.5 * (blon - alon) * hexify::kDegToRad);
    const double h = s_lat * s_lat + std::cos(p1) * std::cos(p2) * s_lon * s_lon;
    return 2.0 * std::asin(std::sqrt(std::min(1.0, h)));
}

// Whether a piece a -> b, whose true midpoint is m, at halving 'depth', must
// be halved again: its chord strays from the edge (chord_needs_split), or it
// spans more than the longest arc allowed.
static inline bool piece_needs_split(double alon, double alat,
                                     double mlon, double mlat,
                                     double blon, double blat,
                                     const EdgeLimits& lim, int depth) {
    if (lim.tolerance > 0.0 && depth < kMaxEdgeSplits &&
        chord_needs_split(alon, alat, mlon, mlat, blon, blat, lim.tolerance)) {
        return true;
    }
    return lim.max_arc > 0.0 && depth < kMaxArcSplits &&
           lonlat_arc(alon, alat, blon, blat) > lim.max_arc;
}

// piece_needs_split() over vectors of pieces at one halving depth.
// [[Rcpp::export]]
LogicalVector cpp_pieces_need_split(NumericVector alon, NumericVector alat,
                                    NumericVector mlon, NumericVector mlat,
                                    NumericVector blon, NumericVector blat,
                                    double tolerance, double max_arc, int depth) {
    const EdgeLimits lim{tolerance, max_arc};
    LogicalVector out(alon.size());
    for (R_xlen_t i = 0; i < alon.size(); i++) {
        out[i] = piece_needs_split(alon[i], alat[i], mlon[i], mlat[i],
                                   blon[i], blat[i], lim, depth);
    }
    return out;
}

// A piece of a vertex cell's folded edge shorter than this fraction of the
// cell's circumradius is a corner, not an edge.
constexpr double kEdgePieceMin = 1e-9;

// A point of the quad plane in lon/lat. A point past an edge of the quad is
// first carried into the quad that owns it, so it is read on its own face as
// a point there is. Only a point past the far vertex is read on the face
// holding the cell centre, extended across its edge.
static void quad_point_lonlat(int quad, double qx, double qy,
                              double qx_center, double qy_center,
                              double& lon, double& lat) {
    lon = NA_REAL;
    lat = NA_REAL;

    int face;
    double tx, ty;
    int own_quad = quad;
    double own_x = qx, own_y = qy;
    if (hexify::quad_xy_canonicalize(own_quad, own_x, own_y) &&
        hexify::try_quad_xy_to_icosa_tri(own_quad, own_x, own_y, face, tx, ty)) {
        auto ll = hexify::face_xy_to_ll(tx, ty, face);
        lon = ll.first;
        lat = ll.second;
        return;
    }

    int center_face;
    double center_tx, center_ty;
    if (hexify::try_quad_xy_to_icosa_tri(quad, qx_center, qy_center,
                                         center_face, center_tx, center_ty)) {
        auto ll = hexify::face_xy_to_ll(center_tx + (qx - qx_center),
                                        center_ty + (qy - qy_center),
                                        center_face);
        lon = ll.first;
        lat = ll.second;
    }
}

// Where the plane segment from a point on a face (ax, ay) to the dropped
// corner (bx, by) of a vertex cell leaves the faces, found by bisection on
// whether the quad's vertex table reads a point.
static void sector_exit(int quad, double ax, double ay, double bx, double by,
                        double& out_x, double& out_y) {
    double lo = 0.0, hi = 1.0;
    int face;
    double tx, ty;
    for (int it = 0; it < 60; it++) {
        double mid = 0.5 * (lo + hi);
        if (hexify::try_quad_xy_to_icosa_tri(quad, ax + mid * (bx - ax),
                                             ay + mid * (by - ay),
                                             face, tx, ty)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    out_x = ax + lo * (bx - ax);
    out_y = ay + lo * (by - ay);
}

// Whether a point of the quad plane lies on a face around the quad's origin.
static inline bool on_a_face(int quad, double x, double y) {
    int face;
    double tx, ty;
    return hexify::try_quad_xy_to_icosa_tri(quad, x, y, face, tx, ty);
}

// Which corners a cell at a vertex of the solid drops: a run of consecutive
// corners, from 'first', 'count' long. n faces meet at the vertex and the quad
// plane holds six sectors, so the 6 - n sectors no face reads are the solid's
// angular deficit there, and the corners in them go. The corners are a sixth
// of a turn apart, as wide as a sector, so either 6 - n corners lie inside the
// deficit, or 7 - n lie in it or on its bounding rays with the edges between
// them inside; then the earliest of those stays, and its folded edge reaches
// the same point of the globe as the edge into the latest.
static void dropped_corners(int quad, const double vx[6], const double vy[6],
                            int n_drop, int& first, int& count) {
    bool out[6];
    for (int c = 0; c < 6; c++) out[c] = !on_a_face(quad, vx[c], vy[c]);
    first = -1;
    count = 0;
    for (int c = 0; c < 6; c++) {
        if (out[c]) { first = c; break; }
    }
    if (first >= 0) {
        count = 1;
        if (n_drop > 1) {
            while (count < 6 && out[(first + 5) % 6]) { first = (first + 5) % 6; count++; }
            while (count < 6 && out[(first + count) % 6]) count++;
        }
    } else {
        for (int c = 0; c < 6; c++) {
            int d = (c + 1) % 6;
            if (!on_a_face(quad, 0.5 * (vx[c] + vx[d]), 0.5 * (vy[c] + vy[d]))) {
                first = d;
                count = 1;
                break;
            }
        }
    }
    if (first < 0) Rcpp::stop("vertex cell has no corner in the solid's deficit");
    while (count < n_drop) {
        int last = (first + count - 1) % 6, next = (first + count) % 6;
        int prev = (first + 5) % 6;
        if (!on_a_face(quad, 0.5 * (vx[last] + vx[next]), 0.5 * (vy[last] + vy[next]))) {
            count++;
        } else if (!on_a_face(quad, 0.5 * (vx[prev] + vx[first]), 0.5 * (vy[prev] + vy[first]))) {
            first = prev;
            count++;
        } else {
            break;
        }
    }
}

// One straight piece of a cell's boundary in the quad plane, a -> b. A piece
// that leaves a vertex cell's dropped sector starts at the point of the globe
// where the other piece entered it; when that point is a corner of the cell,
// 'from_corner' names it so the corner keeps its own reading. 'wall' counts
// the cell's walls from 0, and both pieces of a folded edge carry one wall.
struct PlaneEdge {
    double ax, ay, bx, by;
    bool from_corner;
    double cx, cy;
    int wall;
};

// The boundary of one cell as straight pieces of the quad plane,
// counter-clockwise: one corner per face around the vertex for a cell at a
// vertex of the solid (five on the icosahedron, four on the octahedron), six
// for every other cell. With 'fold_vertex', the edge of a vertex cell across
// its dropped corners is folded onto the faces. The quad plane holds six
// triangles around the vertex and the globe fewer, so the sectors holding
// those corners are no face: the two rays bounding them are one face edge, and
// the plane edges c0 -> (first dropped) and (last dropped) -> c1 meet them at
// one point of the globe, the hexagon being symmetric under the deficit's turn.
// The cell's edge runs from c0 to that point on the first ray and on from the
// second to c1. When the corners sit on the rays instead (a Class II vertex
// cell), one of the two plane edges lies wholly in the deficit, so its piece
// has no length and the corner it leaves from is the same point of the globe
// as the other piece's start. Without 'fold_vertex' the edge joins c0 and c1
// directly.
static void cell_plane_edges(int quad, double qx_center, double qy_center,
                             double radius, double rotation_deg,
                             bool at_vertex, bool fold_vertex,
                             std::vector<PlaneEdge>& out) {
    double vx[6], vy[6];
    hexify::generate_hex_corners(qx_center, qy_center, radius, rotation_deg,
                                 vx, vy);

    int drop_first = -1, drop_count = 0;
    if (at_vertex) {
        dropped_corners(quad, vx, vy, 6 - hexify::topo().valence[quad],
                        drop_first, drop_count);
    }
    auto dropped = [&](int c) {
        return drop_count > 0 && ((c - drop_first + 6) % 6) < drop_count;
    };
    const int drop_last = (drop_first + drop_count - 1 + 6) % 6;

    std::vector<int> kept;
    for (int c = 0; c < 6; c++) {
        if (!dropped(c)) kept.push_back(c);
    }
    int n_corner = static_cast<int>(kept.size());

    out.clear();
    for (int k = 0; k < n_corner; k++) {
        int c0 = kept[k];
        int c1 = kept[(k + 1) % n_corner];
        if (c1 != (c0 + 1) % 6 && fold_vertex) {
            double ax, ay, bx, by;
            sector_exit(quad, vx[c0], vy[c0], vx[drop_first], vy[drop_first], ax, ay);
            sector_exit(quad, vx[c1], vy[c1], vx[drop_last], vy[drop_last], bx, by);
            double piece_min = kEdgePieceMin * radius;
            bool from_c0 = std::hypot(ax - vx[c0], ay - vy[c0]) > piece_min;
            if (from_c0) {
                out.push_back({vx[c0], vy[c0], ax, ay, false, 0.0, 0.0, k});
            }
            if (std::hypot(vx[c1] - bx, vy[c1] - by) > piece_min) {
                out.push_back({bx, by, vx[c1], vy[c1], !from_c0, vx[c0], vy[c0], k});
            }
        } else {
            out.push_back({vx[c0], vy[c0], vx[c1], vy[c1], false, 0.0, 0.0, k});
        }
    }
}

// Where a point of the quad plane lies on the solid: its face and
// triangle coordinates. False for a point past the far vertex of every quad,
// which lies in no face.
static bool quad_point_face(int quad, double qx, double qy,
                            int& face, double& tx, double& ty) {
    return hexify::quad_xy_canonicalize(quad, qx, qy) &&
           hexify::try_quad_xy_to_icosa_tri(quad, qx, qy, face, tx, ty);
}

// A straight piece of a cell edge on one face, in its triangle coordinates,
// the cell wall it belongs to, and where it starts and ends along the plane
// segment it was cut from (0 to 1).
struct FacePiece {
    int face;
    double ax, ay, bx, by;
    int wall;
    double s0, s1;
};

// Snyder's projection has a crease along each line from a face's centre to
// its corners: its derivative jumps there, so the image of a straight piece
// bends where it crosses one. The face triangle is the same in every face's
// triangle coordinates, with unit edge and its first corner at the top.
// Returns how many such lines the piece a -> b crosses strictly between its
// ends, with the crossings' places along it, ascending, in 'u'. Fuller's
// projection is smooth inside a face, so it has none.
static int radius_crossings(double ax, double ay, double bx, double by,
                            double u[3]) {
    if (hexify::active_projection() != hexify::FaceProjection::ISEA) return 0;
    constexpr double kCx = 0.5;
    constexpr double kCy = 0.28867513459481288225;   // 1 / (2 sqrt(3))
    static const double corner[3][2] = {
        {0.5, 0.86602540378443864676}, {0.0, 0.0}, {1.0, 0.0}};
    // A crossing this close to an end is that end.
    constexpr double kEndSlack = 1e-12;
    const double dx = bx - ax, dy = by - ay;
    int n = 0;
    for (int k = 0; k < 3; k++) {
        const double ex = corner[k][0] - kCx, ey = corner[k][1] - kCy;
        const double den = dx * ey - dy * ex;
        if (den == 0.0) continue;
        // a + t d = c + r e, solved for t (along the piece) and r (along the
        // line from the centre, 0 to 1)
        const double wx = kCx - ax, wy = kCy - ay;
        const double t = (wx * ey - wy * ex) / den;
        const double r = (wx * dy - wy * dx) / den;
        if (t > kEndSlack && t < 1.0 - kEndSlack && r >= 0.0 && r <= 1.0) u[n++] = t;
    }
    for (int i = 1; i < n; i++) {
        for (int j = i; j > 0 && u[j] < u[j - 1]; j--) std::swap(u[j], u[j - 1]);
    }
    return n;
}

// Each face piece cut where it crosses a crease of the projection
// (radius_crossings), so a vertex lies on every crease it crosses.
static void split_at_creases(std::vector<FacePiece>& pieces) {
    std::vector<FacePiece> out;
    out.reserve(pieces.size());
    double u[3];
    for (const FacePiece& p : pieces) {
        int n = radius_crossings(p.ax, p.ay, p.bx, p.by, u);
        double prev_u = 0.0, px = p.ax, py = p.ay;
        for (int i = 0; i <= n; i++) {
            double t = i < n ? u[i] : 1.0;
            double qx = i < n ? p.ax + t * (p.bx - p.ax) : p.bx;
            double qy = i < n ? p.ay + t * (p.by - p.ay) : p.by;
            out.push_back({p.face, px, py, qx, qy, p.wall,
                           p.s0 + prev_u * (p.s1 - p.s0), p.s0 + t * (p.s1 - p.s0)});
            prev_u = t;
            px = qx;
            py = qy;
        }
    }
    pieces.swap(out);
}

// The plane segment a -> b cut where it crosses from one face to the next.
// Each face is a convex triangle of the plane, so the segment's points on the
// face it is in form one interval: bisection on face membership finds where
// the segment leaves it, and the walk goes on from the face it enters there.
// Returns NULL, or what went wrong when part of the segment lies on no face.
static const char* plane_segment_faces(int quad, double ax, double ay,
                                       double bx, double by, int wall,
                                       std::vector<FacePiece>& out) {
    int face_b;
    double bt_x, bt_y;
    if (!quad_point_face(quad, bx, by, face_b, bt_x, bt_y)) {
        return "cell edge ends outside every face of the solid";
    }
    double s = 0.0;
    int face;
    double tx, ty;
    if (!quad_point_face(quad, ax, ay, face, tx, ty)) {
        return "cell edge starts outside every face of the solid";
    }
    // Two faces meet along an edge and at most five around a vertex, so a segment
    // shorter than a face crosses at most a few.
    for (int guard = 0; guard < 8; guard++) {
        if (face == face_b) {
            out.push_back({face, tx, ty, bt_x, bt_y, wall, s, 1.0});
            return nullptr;
        }
        double lo = s, hi = 1.0;
        int f;
        double fx, fy;
        for (int it = 0; it < 60; it++) {
            double mid = 0.5 * (lo + hi);
            if (quad_point_face(quad, ax + mid * (bx - ax), ay + mid * (by - ay),
                                f, fx, fy) && f == face) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double ex, ey;
        quad_point_face(quad, ax + lo * (bx - ax), ay + lo * (by - ay), f, ex, ey);
        out.push_back({face, tx, ty, ex, ey, wall, s, lo});
        if (!quad_point_face(quad, ax + hi * (bx - ax), ay + hi * (by - ay),
                             face, tx, ty)) {
            return "cell edge passes outside every face of the solid";
        }
        s = hi;
    }
    return "cell edge crosses more faces than a cell edge can";
}

// The boundary of one cell in lon/lat, counter-clockwise and left open. A
// cell edge is straight in the quad plane and curved in lon/lat, so with
// active limits each edge is split in the plane, as DGGRID's densification
// does, until every piece is a straight lon/lat chord to within the
// tolerance's fraction of its length and spans no more than the longest arc
// (piece_needs_split). Inactive limits give the corners alone. A pole falls
// inside a cell or on a cell edge, never on a corner, so every corner keeps
// the position the inverse projection gives it.
static void cell_boundary_lonlat(int quad, double qx_center, double qy_center,
                                 double radius, double rotation_deg,
                                 bool at_vertex, const EdgeLimits& lim,
                                 std::vector<double>& out_lon,
                                 std::vector<double>& out_lat) {
    const bool dense = lim.active();
    std::vector<PlaneEdge> edges;
    cell_plane_edges(quad, qx_center, qy_center, radius, rotation_deg,
                     at_vertex, /*fold_vertex=*/dense, edges);

    out_lon.clear();
    out_lat.clear();

    auto project = [&](double x, double y, double& lon, double& lat) {
        quad_point_lonlat(quad, x, y, qx_center, qy_center, lon, lat);
    };

    // The plane segment a -> b in lon/lat, from a and short of b, halved while
    // piece_needs_split() holds.
    std::function<void(double, double, double, double, double, double,
                       double, double, int)> add_segment;
    add_segment = [&](double ax, double ay, double alon, double alat,
                      double bx, double by, double blon, double blat,
                      int depth) {
        if (dense && R_finite(alon) && R_finite(blon)) {
            double mx = 0.5 * (ax + bx), my = 0.5 * (ay + by);
            double mlon, mlat;
            project(mx, my, mlon, mlat);
            if (R_finite(mlon) &&
                piece_needs_split(alon, alat, mlon, mlat, blon, blat, lim, depth)) {
                add_segment(ax, ay, alon, alat, mx, my, mlon, mlat, depth + 1);
                add_segment(mx, my, mlon, mlat, bx, by, blon, blat, depth + 1);
                return;
            }
        }
        out_lon.push_back(alon);
        out_lat.push_back(alat);
    };
    // Cuts closer than this, as a share of the edge, to each other or to an
    // end are one point.
    constexpr double kCutSlack = 1e-9;
    // Densified, an edge is first cut where it crosses a face edge or a crease
    // of the projection (split_at_creases), where its image bends, so a
    // vertex lies on each bend.
    std::vector<FacePiece> pieces;
    auto add_edge = [&](double ax, double ay, double bx, double by) {
        double alon, alat, blon, blat;
        project(ax, ay, alon, alat);
        project(bx, by, blon, blat);
        // An edge over a pole is two meridians, already exact as drawn, and
        // splitting it lands on the pole, where longitude is undefined.
        if (std::fabs(std::fabs(blon - alon) - 180.0) < kPoleEdgeLonSlack) {
            out_lon.push_back(alon);
            out_lat.push_back(alat);
            return;
        }
        pieces.clear();
        if (dense &&
            plane_segment_faces(quad, ax, ay, bx, by, 0, pieces) == nullptr) {
            split_at_creases(pieces);
        } else {
            pieces.clear();
        }
        double px = ax, py = ay, plon = alon, plat = alat;
        double last = 0.0;
        for (size_t i = 0; i + 1 < pieces.size(); i++) {
            const double t = pieces[i].s1;
            // A corner on a face edge leaves a piece of rounding length
            // there; a cut at it would be a second copy of the corner.
            if (t <= last + kCutSlack || t >= 1.0 - kCutSlack) continue;
            last = t;
            const double qx = ax + t * (bx - ax), qy = ay + t * (by - ay);
            double qlon, qlat;
            project(qx, qy, qlon, qlat);
            add_segment(px, py, plon, plat, qx, qy, qlon, qlat, 0);
            px = qx;
            py = qy;
            plon = qlon;
            plat = qlat;
        }
        add_segment(px, py, plon, plat, bx, by, blon, blat, 0);
    };

    for (const PlaneEdge& e : edges) {
        size_t start = out_lon.size();
        add_edge(e.ax, e.ay, e.bx, e.by);
        if (e.from_corner) {
            project(e.cx, e.cy, out_lon[start], out_lat[start]);
        }
    }
}

// A boundary as an (n + 1) x 2 lon/lat matrix, the first point repeated last.
static NumericMatrix closed_ring(const std::vector<double>& lon,
                                 const std::vector<double>& lat) {
    int n = static_cast<int>(lon.size());
    NumericMatrix coords(n + 1, 2);
    for (int v = 0; v < n; v++) {
        coords(v, 0) = lon[v];
        coords(v, 1) = lat[v];
    }
    coords(n, 0) = lon[0];
    coords(n, 1) = lat[0];
    colnames(coords) = CharacterVector::create("lon", "lat");
    return coords;
}

// ============================================================================
// Cell Boundaries
// ============================================================================

// A cell placed in the quad plane: its quad, centre, and whether it sits at an
// icosahedral vertex.
struct CellPlane {
    int quad;
    double qx, qy;
    bool at_vertex;
};

// The cells of a grid in the quad plane, with the hexagon's circumradius and
// turn there.
struct CellPlanes {
    double radius;
    double rotation_deg;
    std::vector<CellPlane> cells;
};

static CellPlanes cell_planes(const NumericVector& cell_id, const QuadFrame& f) {
    CellPlanes out;
    out.radius = kHexCircumradius / f.form.scale;
    out.rotation_deg = hexify::form_rotation_deg(f.form);
    out.cells.resize(cell_id.size());

    for (R_xlen_t k = 0; k < cell_id.size(); k++) {
        CellPlane& c = out.cells[k];
        long long i, j;
        frame_decode(f, cell_id[k], c.quad, i, j);
        frame_ij_to_xy(f, c.quad, i, j, c.qx, c.qy);
        c.at_vertex = (i == 0 && j == 0);
    }
    return out;
}

static List cell_rings(const NumericVector& cell_id, const QuadFrame& f,
                       const EdgeLimits& lim) {
    CellPlanes g = cell_planes(cell_id, f);
    List result(cell_id.size());
    std::vector<double> lon, lat;
    for (R_xlen_t k = 0; k < cell_id.size(); k++) {
        const CellPlane& c = g.cells[k];
        cell_boundary_lonlat(c.quad, c.qx, c.qy, g.radius, g.rotation_deg,
                             c.at_vertex, lim, lon, lat);
        result[k] = closed_ring(lon, lat);
    }
    return result;
}

// [[Rcpp::export]]
List cpp_cell_to_corners(NumericVector icosa, NumericVector cell_id,
                         int resolution, int aperture, IntegerVector ap_seq,
                         double tolerance = 0.0, double max_arc = 0.0) {
    activate_grid(icosa);
    return cell_rings(cell_id, grid_frame(resolution, aperture, ap_seq),
                      EdgeLimits{tolerance, max_arc});
}

// Closed lon/lat rings whose edges are great-circle arcs between corners, as
// H3 draws its cells, with each arc halved while piece_needs_split() holds,
// as ISEA cell edges are.
// [[Rcpp::export]]
List cpp_densify_great_circle(List rings, double tolerance, double max_arc = 0.0) {
    const EdgeLimits lim{tolerance, max_arc};
    auto unit = [](double lon, double lat, double v[3]) {
        double la = lat * hexify::kDegToRad, lo = lon * hexify::kDegToRad;
        v[0] = std::cos(la) * std::cos(lo);
        v[1] = std::cos(la) * std::sin(lo);
        v[2] = std::sin(la);
    };
    List out(rings.size());
    std::vector<double> lon, lat;
    for (R_xlen_t r = 0; r < rings.size(); r++) {
        NumericMatrix ring = rings[r];
        int n = ring.nrow();
        if (n < 2 || !lim.active()) {
            out[r] = ring;
            continue;
        }
        lon.clear();
        lat.clear();
        std::function<void(const double*, double, double, const double*,
                           double, double, int)> add_arc;
        add_arc = [&](const double* a, double alon, double alat,
                      const double* b, double blon, double blat, int depth) {
            double m[3] = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
            double norm = std::sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
            if (norm > 0.0) {
                for (double& x : m) x /= norm;
                double mlon = std::atan2(m[1], m[0]) * hexify::kRadToDeg;
                double mlat = std::asin(std::max(-1.0, std::min(1.0, m[2]))) *
                              hexify::kRadToDeg;
                if (piece_needs_split(alon, alat, mlon, mlat, blon, blat,
                                      lim, depth)) {
                    add_arc(a, alon, alat, m, mlon, mlat, depth + 1);
                    add_arc(m, mlon, mlat, b, blon, blat, depth + 1);
                    return;
                }
            }
            lon.push_back(alon);
            lat.push_back(alat);
        };
        double a[3], b[3];
        for (int i = 0; i + 1 < n; i++) {
            double alon = ring(i, 0), alat = ring(i, 1);
            double blon = ring(i + 1, 0), blat = ring(i + 1, 1);
            if (std::fabs(std::fabs(blon - alon) - 180.0) < kPoleEdgeLonSlack) {
                lon.push_back(alon);
                lat.push_back(alat);
                continue;
            }
            unit(alon, alat, a);
            unit(blon, blat, b);
            add_arc(a, alon, alat, b, blon, blat, 0);
        }
        lon.push_back(ring(n - 1, 0));
        lat.push_back(ring(n - 1, 1));
        NumericMatrix dense(static_cast<int>(lon.size()), 2);
        for (size_t k = 0; k < lon.size(); k++) {
            dense(k, 0) = lon[k];
            dense(k, 1) = lat[k];
        }
        out[r] = dense;
    }
    return out;
}

// The boundary of one cell as straight pieces on the faces it covers, in
// their triangle coordinates, counter-clockwise in the quad plane.
static void cell_face_pieces(const CellPlanes& g, const CellPlane& c,
                             std::vector<PlaneEdge>& edges,
                             std::vector<FacePiece>& pieces) {
    cell_plane_edges(c.quad, c.qx, c.qy, g.radius, g.rotation_deg,
                     c.at_vertex, /*fold_vertex=*/true, edges);
    pieces.clear();
    for (const PlaneEdge& e : edges) {
        const char* err = plane_segment_faces(c.quad, e.ax, e.ay, e.bx, e.by, e.wall, pieces);
        if (err) Rcpp::stop(err);
    }
    split_at_creases(pieces);
}

// Deepest halving of one face piece when a wall is measured, 2^20 pieces. A
// piece across a cusp of the projection halves at the cusp until here.
constexpr int kMaxMeasureSplits = 20;

// The face piece a -> b on the unit sphere as a wall in the form cell_walls.h
// describes, from a and short of b: a, then the point the piece's midpoint
// projects to, the piece halved while that point lies further than
// 'tolerance' times the chord from the plane of the great circle through a
// and b. A point's place along the curve does not enter the test.
static void face_piece_sphere(int face, double ax, double ay, const hexify::UnitVec& a,
                              double bx, double by, const hexify::UnitVec& b,
                              double tolerance, int depth,
                              std::vector<hexify::UnitVec>& out) {
    const double mx = 0.5 * (ax + bx), my = 0.5 * (ay + by);
    hexify::UnitVec m;
    hexify::face_tri_to_sphere(face, mx, my, m.data());
    if (depth < kMaxMeasureSplits) {
        // n = a x (b - a) = a x b, whose length is the chord times the cosine
        // of half the arc, the chord to within a part in 10^4 for any piece of
        // a cell wall. m . n = (m - a) . n, as a . n = 0: formed from the short
        // vectors m - a and b - a, its rounding error scales with the chord,
        // so the test below stays meaningful however short the piece.
        const double ux = m[0] - a[0], uy = m[1] - a[1], uz = m[2] - a[2];
        const double vx = b[0] - a[0], vy = b[1] - a[1], vz = b[2] - a[2];
        const double nx = a[1] * vz - a[2] * vy;
        const double ny = a[2] * vx - a[0] * vz;
        const double nz = a[0] * vy - a[1] * vx;
        const double chord = std::sqrt(nx * nx + ny * ny + nz * nz);
        const double off = std::fabs(ux * nx + uy * ny + uz * nz);
        if (off > tolerance * chord * chord) {
            face_piece_sphere(face, ax, ay, a, mx, my, m, tolerance, depth + 1, out);
            face_piece_sphere(face, mx, my, m, bx, by, b, tolerance, depth + 1, out);
            return;
        }
    }
    out.push_back(a);
    out.push_back(m);
}

// The walls of one cell on the unit sphere in the form cell_walls.h
// describes, each from one corner to the next, counter-clockwise, densified to
// 'tolerance' (face_piece_sphere).
static void cell_walls_sphere(const CellPlanes& g, const CellPlane& c, double tolerance,
                              std::vector<PlaneEdge>& edges,
                              std::vector<FacePiece>& pieces,
                              std::vector<std::vector<hexify::UnitVec>>& walls) {
    cell_face_pieces(g, c, edges, pieces);
    int n_wall = 0;
    for (const FacePiece& p : pieces) n_wall = std::max(n_wall, p.wall + 1);
    walls.assign(n_wall, {});
    hexify::UnitVec a, b;
    for (size_t i = 0; i < pieces.size(); i++) {
        const FacePiece& p = pieces[i];
        hexify::face_tri_to_sphere(p.face, p.ax, p.ay, a.data());
        hexify::face_tri_to_sphere(p.face, p.bx, p.by, b.data());
        face_piece_sphere(p.face, p.ax, p.ay, a, p.bx, p.by, b, tolerance, 0,
                          walls[p.wall]);
        if (i + 1 == pieces.size() || pieces[i + 1].wall != p.wall) {
            walls[p.wall].push_back(b);
        }
    }
}

// A cell's centre on the unit sphere.
static hexify::UnitVec cell_centre_sphere(const CellPlane& c) {
    int face;
    double tx, ty;
    if (!hexify::try_quad_xy_to_icosa_tri(c.quad, c.qx, c.qy, face, tx, ty)) {
        Rcpp::stop("cell centre lies on no face of the solid");
    }
    hexify::UnitVec v;
    hexify::face_tri_to_sphere(face, tx, ty, v.data());
    return v;
}

// Solid angle of each cell, in steradians: its walls are followed on the
// sphere to 'tolerance' (cell_walls_sphere) and the area they enclose summed
// (enclosed_solid_angle).
static NumericVector cell_solid_angles(const NumericVector& cell_id,
                                       const QuadFrame& f, double tolerance) {
    CellPlanes g = cell_planes(cell_id, f);
    NumericVector out(cell_id.size());
    std::vector<PlaneEdge> edges;
    std::vector<FacePiece> pieces;
    std::vector<std::vector<hexify::UnitVec>> walls;
    for (R_xlen_t k = 0; k < cell_id.size(); k++) {
        cell_walls_sphere(g, g.cells[k], tolerance, edges, pieces, walls);
        out[k] = hexify::enclosed_solid_angle(walls);
    }
    return out;
}

// [[Rcpp::export]]
NumericVector cpp_cell_solid_angle(NumericVector icosa, NumericVector cell_id,
                                   int resolution, int aperture, IntegerVector ap_seq,
                                   double tolerance) {
    activate_grid(icosa);
    return cell_solid_angles(cell_id, grid_frame(resolution, aperture, ap_seq),
                             tolerance);
}

// The boundaries of cells on the flat solid and on the sphere, from the same
// points. Each cell edge is cut where it crosses a face edge, and every piece
// is split into steps no longer than 'step' in triangle coordinates (a face
// edge is about 1). A point carries its position on the flat face and on the
// unit sphere, so the two surfaces show the same boundaries. One closed path
// per cell: a matrix with columns cell (position in 'cell_id', from 1), face
// (from 0), solid x, y, z, sphere x, y, z, and the triangle coordinates tx, ty
// on the face, which a layout of the unfolded solid places in the plane. Where
// a path crosses onto another face the crossing point appears on both faces,
// so the points of one face run unbroken, also where a layout puts faces that
// meet on the solid apart.
// [[Rcpp::export]]
NumericMatrix cpp_cell_surface_paths(NumericVector icosa,
                                     NumericVector cell_id, int resolution,
                                     int aperture, IntegerVector ap_seq,
                                     double step) {
    activate_grid(icosa);
    if (!(step > 0.0)) stop("step must be positive");
    CellPlanes g = cell_planes(cell_id, grid_frame(resolution, aperture, ap_seq));

    constexpr int n_col = 10;
    std::vector<double> rows;
    std::vector<PlaneEdge> edges;
    std::vector<FacePiece> pieces;
    double solid[3], sphere[3];
    auto emit = [&](R_xlen_t cell, int face, double tx, double ty) {
        hexify::face_tri_to_solid(face, tx, ty, solid);
        hexify::face_tri_to_sphere(face, tx, ty, sphere);
        rows.insert(rows.end(), {static_cast<double>(cell + 1),
                                 static_cast<double>(face),
                                 solid[0], solid[1], solid[2],
                                 sphere[0], sphere[1], sphere[2], tx, ty});
    };

    for (R_xlen_t k = 0; k < cell_id.size(); k++) {
        cell_face_pieces(g, g.cells[k], edges, pieces);
        for (size_t i = 0; i < pieces.size(); i++) {
            const FacePiece& p = pieces[i];
            double len = std::hypot(p.bx - p.ax, p.by - p.ay);
            int n = std::max(1, static_cast<int>(std::ceil(len / step)));
            for (int s = 0; s < n; s++) {
                double f = static_cast<double>(s) / n;
                emit(k, p.face, p.ax + f * (p.bx - p.ax), p.ay + f * (p.by - p.ay));
            }
            // Where the path crosses onto another face, its end on this face
            // closes the run of points this face holds.
            const FacePiece& next = pieces[(i + 1) % pieces.size()];
            if (next.face != p.face) emit(k, p.face, p.bx, p.by);
        }
        const FacePiece& first = pieces.front();
        emit(k, first.face, first.ax, first.ay);
    }

    R_xlen_t n_row = static_cast<R_xlen_t>(rows.size() / n_col);
    NumericMatrix out(n_row, n_col);
    for (R_xlen_t r = 0; r < n_row; r++) {
        for (int col = 0; col < n_col; col++) out(r, col) = rows[r * n_col + col];
    }
    colnames(out) = CharacterVector::create("cell", "face", "solid_x", "solid_y",
                                            "solid_z", "sphere_x", "sphere_y",
                                            "sphere_z", "tx", "ty");
    return out;
}

// The solid: its vertices on the unit sphere and the vertex indices (from 1)
// of each face.
// [[Rcpp::export]]
List cpp_icosa_solid(NumericVector icosa) {
    activate_icosa(icosa);
    const hexify::PolyData& S = hexify::poly();
    NumericMatrix verts(S.n_verts(), 3);
    for (int v = 0; v < S.n_verts(); v++) {
        double cl = std::cos(S.verts[v].lat);
        verts(v, 0) = cl * std::cos(S.verts[v].lon);
        verts(v, 1) = cl * std::sin(S.verts[v].lon);
        verts(v, 2) = std::sin(S.verts[v].lat);
    }
    IntegerMatrix faces(S.n_faces(), 3);
    for (int f = 0; f < S.n_faces(); f++) {
        for (int k = 0; k < 3; k++) faces(f, k) = S.topo->faces[f][k] + 1;
    }
    return List::create(_["vertices"] = verts, _["faces"] = faces);
}

// Points of one face, given in lon/lat, as triangle coordinates and as their
// place on the flat face. A point off the face lands on the face's plane
// extended, where the projection still reads it.
// [[Rcpp::export]]
NumericMatrix cpp_lonlat_to_face_solid(NumericVector icosa,
                                       int face, NumericVector lon,
                                       NumericVector lat) {
    activate_icosa(icosa);
    if (face < 0 || face >= hexify::poly().n_faces()) stop("face out of range for the solid");
    R_xlen_t n = lon.size();
    NumericMatrix out(n, 5);
    double p[3];
    for (R_xlen_t k = 0; k < n; k++) {
        auto t = hexify::snyder_forward_to_face(face, lon[k], lat[k]);
        hexify::face_tri_to_solid(face, t.first, t.second, p);
        out(k, 0) = t.first;
        out(k, 1) = t.second;
        out(k, 2) = p[0];
        out(k, 3) = p[1];
        out(k, 4) = p[2];
    }
    colnames(out) = CharacterVector::create("tx", "ty", "x", "y", "z");
    return out;
}

// Triangle coordinates of a face as their place on the flat face.
// [[Rcpp::export]]
NumericMatrix cpp_face_tri_to_solid(NumericVector icosa,
                                    int face, NumericVector tx, NumericVector ty) {
    activate_icosa(icosa);
    if (face < 0 || face >= hexify::poly().n_faces()) stop("face out of range for the solid");
    R_xlen_t n = tx.size();
    NumericMatrix out(n, 3);
    double p[3];
    for (R_xlen_t k = 0; k < n; k++) {
        hexify::face_tri_to_solid(face, tx[k], ty[k], p);
        for (int r = 0; r < 3; r++) out(k, r) = p[r];
    }
    colnames(out) = CharacterVector::create("x", "y", "z");
    return out;
}

// ============================================================================
// Neighbor Finding (v0.7.0)
// ============================================================================

// What a renderer needs to find the cell of a quad-plane point by itself.
// Scaled by 'dim', the quad's side in substrate steps, a point's nearest cell
// centre is its nearest multiple of the generator a + b*omega (omega =
// exp(2*pi*i/3)) in the substrate's (i, j); the edge table moves that centre
// into the quad that owns it, and the cell ID counts 'per_quad' cells per
// quad, numbered within a quad as cell_index_2d() numbers them on the
// sublattice j = c * i (mod index). Aperture 7 stores surrogates but numbers
// its cells by their substrate centres, which is the same count.
// [[Rcpp::export]]
List cpp_globe_frame(NumericVector icosa, int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    SubstrateLattice lattice = sublattice_of(f.form);
    LatticeGenerator generator = generator_of(f.form);
    return List::create(
        _["dim"] = static_cast<double>(f.dim),
        _["index"] = static_cast<double>(lattice.index),
        _["c"] = static_cast<double>(lattice.c),
        _["generator"] = NumericVector::create(static_cast<double>(generator.a),
                                               static_cast<double>(generator.b)),
        _["per_quad"] = static_cast<double>(f.offsetPerQuad),
        _["n_cells"] = static_cast<double>(f.nCells));
}

// The cells adjacent to a vertex quad's cell. The vertex is a corner of every
// diamond quad around it, and each of its neighbours lies one lattice step
// from that corner inside one of those quads' boxes, so the neighbours are the
// steps from the corner that land in the box, over every quad the vertex is a
// corner of. No step leaves a box, so none is read across a face edge.
static void pole_neighbors(const QuadFrame& f, int pole, const long long offsets[6][2],
                           std::vector<double>& out) {
    const hexify::SolidTopology& t = hexify::topo();
    for (int q = 1; q <= t.n_diamonds(); q++) {
        const auto& c = t.corner[q];
        long long ci, cj;
        if (c[hexify::kCornerJ] == pole) { ci = 0; cj = f.dim; }
        else if (c[hexify::kCornerI] == pole) { ci = f.dim; cj = 0; }
        else if (c[hexify::kCornerFar] == pole) { ci = f.dim; cj = f.dim; }
        else continue;
        long long i = ci, j = cj;
        if (f.aperture == 7) {
            hexify::ap7_substrate_to_surrogate_ijk(ci, cj, f.resolution, i, j);
        }
        for (int d = 0; d < 6; d++) {
            long long ni = i + offsets[d][0];
            long long nj = j + offsets[d][1];
            if (frame_in_quad(f, ni, nj)) out.push_back(frame_encode(f, q, ni, nj));
        }
    }
}

// The six cells adjacent to each of `cell_id`, in the frame's own grid.
//
// The six neighbours are the generator times the six units of the Eisenstein
// integers, so one step table serves every lattice. A step leaving the quad is
// sent back through the forward pipeline, which names the quad that owns it.
static Rcpp::List neighbors_in_frame(const Rcpp::NumericVector& cell_id,
                                      const QuadFrame& f) {
    int n = cell_id.size();
    Rcpp::List out(n);

    long long offsets[6][2];
    lattice_unit_steps(f.generator, offsets);

    for (int k = 0; k < n; k++) {
        uint64_t idx = frame_cell_index(f, cell_id[k]);

        std::vector<double> neighbor_ids;
        neighbor_ids.reserve(6);

        // Resolution 0 is one base cell per vertex of the solid, a vertex cell
        // each, and each quad holds a single cell, so adjacency there is the
        // solid's vertex graph rather than a step through a quad frame.
        if (f.resolution == 0) {
            const std::vector<int>& nb = hexify::topo().neighbors[idx];
            Rcpp::NumericVector base(nb.size());
            for (size_t d = 0; d < nb.size(); d++) {
                base[d] = nb[d] + 1;
            }
            out[k] = base;
            continue;
        }

        int quad;
        long long i, j;

        if (idx == 0 || idx == f.nCells - 1) {
            // The two vertex quads hold a single cell each -- a vertex of the
            // solid where several quads meet -- so their own frame carries no
            // offsets to step through.
            pole_neighbors(f, idx == 0 ? 0 : hexify::topo().south_pole(), offsets,
                           neighbor_ids);
            std::sort(neighbor_ids.begin(), neighbor_ids.end());
            neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()),
                               neighbor_ids.end());
            out[k] = Rcpp::NumericVector(neighbor_ids.begin(), neighbor_ids.end());
            continue;
        }
        frame_decode_index(f, idx, quad, i, j);

        for (int d = 0; d < 6; d++) {
            long long ni = i + offsets[d][0];
            long long nj = j + offsets[d][1];

            if (frame_in_quad(f, ni, nj)) {
                neighbor_ids.push_back(frame_encode(f, quad, ni, nj));
                continue;
            }

            // The neighbour sits in another quad: send its centre back through
            // the forward pipeline, which names the quad that owns it.
            double nbr_qx, nbr_qy;
            frame_ij_to_xy(f, quad, ni, nj, nbr_qx, nbr_qy);

            int tri_face;
            double tri_x, tri_y;
            if (!hexify::try_quad_xy_to_icosa_tri(quad, nbr_qx, nbr_qy,
                                                   tri_face, tri_x, tri_y)) {
                // Under-runs of a quad's own frame have no image in it, so the
                // projection has nowhere to send them. Step to the owning quad
                // through the edge table instead.
                int alt_quad = quad;
                long long alt_i = ni, alt_j = nj;
                if (frame_canonicalize(f, alt_quad, alt_i, alt_j)) {
                    neighbor_ids.push_back(frame_encode(f, alt_quad, alt_i, alt_j));
                }
                continue;
            }

            auto ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
            int final_quad;
            long long final_i, final_j;
            frame_locate(f, ll.first, ll.second, final_quad, final_i, final_j);

            neighbor_ids.push_back(frame_encode(f, final_quad, final_i, final_j));
        }

        // Remove duplicates (can happen at pentagons / boundary)
        std::sort(neighbor_ids.begin(), neighbor_ids.end());
        neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()),
                           neighbor_ids.end());
        // Remove self
        double self_id = cell_id[k];
        neighbor_ids.erase(
            std::remove(neighbor_ids.begin(), neighbor_ids.end(), self_id),
            neighbor_ids.end());

        out[k] = Rcpp::NumericVector(neighbor_ids.begin(), neighbor_ids.end());
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_get_neighbors_isea(NumericVector icosa, Rcpp::NumericVector cell_id,
                                  int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    return neighbors_in_frame(cell_id, grid_frame(resolution, aperture, ap_seq));
}

// Perimeter of each cell, on the unit sphere, its walls followed to within
// 'tolerance' (face_piece_sphere). With 'walls', also one row per wall: the
// cell (position in 'cell_id', from 1), the neighbour across it, and its
// WallMeasure.
// [[Rcpp::export]]
List cpp_cell_walls(NumericVector icosa, NumericVector cell_id, int resolution,
                    int aperture, IntegerVector ap_seq, double tolerance, bool walls) {
    activate_grid(icosa);
    if (!(tolerance > 0.0)) stop("tolerance must be positive");
    const QuadFrame f = grid_frame(resolution, aperture, ap_seq);
    const CellPlanes g = cell_planes(cell_id, f);
    const R_xlen_t n = cell_id.size();

    List nbrs = walls ? neighbors_in_frame(cell_id, f) : List(0);
    NumericVector perimeter(n);
    hexify::WallRows rows;
    std::vector<double> row_nbr;
    std::vector<PlaneEdge> edges;
    std::vector<FacePiece> pieces;
    std::vector<std::vector<hexify::UnitVec>> lines;
    std::vector<hexify::WallShape> shapes;
    std::vector<hexify::UnitVec> nbr_centres;
    std::vector<int> wall_of;

    for (R_xlen_t k = 0; k < n; k++) {
        cell_walls_sphere(g, g.cells[k], tolerance, edges, pieces, lines);
        shapes.clear();
        double p = 0.0;
        for (const auto& line : lines) {
            shapes.push_back(hexify::wall_shape(line));
            p += shapes.back().length;
        }
        perimeter[k] = p;
        if (!walls) continue;

        NumericVector nb = nbrs[k];
        const CellPlanes ng = cell_planes(nb, f);
        nbr_centres.clear();
        for (const CellPlane& c : ng.cells) nbr_centres.push_back(cell_centre_sphere(c));
        if (!hexify::walls_of_neighbours(shapes, nbr_centres, wall_of)) {
            stop("cell %.0f: its %d walls do not pair one to one with its %d neighbours",
                 cell_id[k], static_cast<int>(shapes.size()), static_cast<int>(nb.size()));
        }
        const hexify::UnitVec centre = cell_centre_sphere(g.cells[k]);
        for (R_xlen_t j = 0; j < nb.size(); j++) {
            rows.add(k, hexify::measure_wall(shapes[wall_of[j]], centre, nbr_centres[j]));
            row_nbr.push_back(nb[j]);
        }
    }

    List out = List::create(_["perimeter"] = perimeter);
    if (walls) {
        out["walls"] = rows.frame("neighbor_id",
                                  NumericVector(row_nbr.begin(), row_nbr.end()));
    }
    return out;
}

// The same measures for cells given as closed lon/lat rings whose walls are
// great-circle arcs between consecutive corners, such as another program's
// corner output: one row per neighbour, the wall chosen as walls_of_neighbours()
// chooses it. 'centres' is a two-column lon/lat matrix, one row per ring, and
// 'neighbour_centres' a list of such matrices; the 'neighbor' column is the
// row of the neighbour's centre there.
// [[Rcpp::export]]
DataFrame cpp_ring_walls(List rings, NumericMatrix centres, List neighbour_centres) {
    hexify::WallRows rows;
    std::vector<double> row_nbr;
    std::vector<hexify::WallShape> shapes;
    std::vector<hexify::UnitVec> nbr_centres;
    std::vector<int> wall_of;
    for (R_xlen_t k = 0; k < rings.size(); k++) {
        NumericMatrix ring = rings[k];
        NumericMatrix nc = neighbour_centres[k];
        shapes.clear();
        for (int i = 0; i + 1 < ring.nrow(); i++) {
            shapes.push_back(hexify::wall_shape(hexify::great_circle_wall(
                {hexify::unit_from_lonlat(ring(i, 0), ring(i, 1)),
                 hexify::unit_from_lonlat(ring(i + 1, 0), ring(i + 1, 1))})));
        }
        nbr_centres.clear();
        for (int j = 0; j < nc.nrow(); j++) {
            nbr_centres.push_back(hexify::unit_from_lonlat(nc(j, 0), nc(j, 1)));
        }
        if (!hexify::walls_of_neighbours(shapes, nbr_centres, wall_of)) {
            stop("ring %d: its %d walls do not pair one to one with its %d neighbours",
                 static_cast<int>(k + 1), static_cast<int>(shapes.size()),
                 static_cast<int>(nbr_centres.size()));
        }
        const hexify::UnitVec centre = hexify::unit_from_lonlat(centres(k, 0), centres(k, 1));
        for (size_t j = 0; j < nbr_centres.size(); j++) {
            rows.add(k, hexify::measure_wall(shapes[wall_of[j]], centre, nbr_centres[j]));
            row_nbr.push_back(static_cast<double>(j + 1));
        }
    }
    return rows.frame("neighbor", NumericVector(row_nbr.begin(), row_nbr.end()));
}

// ============================================================================
// PLANE Coordinate Conversions
// ============================================================================
// PLANE coordinates represent the unfolded solid in 2D.
// The transformation from Icosa Triangle to PLANE involves:
// 1. Rotate the point by rot60 * 60 degrees
// 2. Translate by the triangle's offset position
//
// On the icosahedron this is a flat map layout ~5.5 x 1.73 units containing
// all 20 triangles.
// Produces standard ISEA PLANE coordinates for visualization.
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_icosa_tri_to_plane(NumericVector icosa, IntegerVector icosa_triangle_face,
                                  NumericVector icosa_triangle_x,
                                  NumericVector icosa_triangle_y) {
    activate_icosa(icosa);
    int n = icosa_triangle_face.size();
    NumericVector out_px(n);
    NumericVector out_py(n);

    for (int k = 0; k < n; k++) {
        int face = icosa_triangle_face[k];
        if (face < 0 || face >= hexify::poly().n_faces()) {
            out_px[k] = NA_REAL;
            out_py[k] = NA_REAL;
            continue;
        }

        hexify::face_tri_to_plane(face, icosa_triangle_x[k], icosa_triangle_y[k],
                                  out_px[k], out_py[k]);
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_plane(NumericVector icosa, NumericVector cell_id, int resolution, int aperture) {
    activate_grid(icosa);
    QuadFrame f = quad_frame(resolution, aperture, {});
    R_xlen_t n = cell_id.size();
    NumericVector out_px(n);
    NumericVector out_py(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);

        double quad_x, quad_y;
        frame_ij_to_xy(f, quad, i, j, quad_x, quad_y);

        int tri_face;
        double tri_x, tri_y;
        hexify::quad_xy_to_icosa_tri(quad, quad_x, quad_y, tri_face,
                                     tri_x, tri_y);

        hexify::face_tri_to_plane(tri_face, tri_x, tri_y, out_px[k], out_py[k]);
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}

// [[Rcpp::export]]
DataFrame cpp_lonlat_to_plane(NumericVector icosa,
                              NumericVector lon, NumericVector lat) {
    activate_icosa(icosa);
    int n = lon.size();
    if (lat.size() != n) {
        stop("cpp_lonlat_to_plane: lon and lat must have same length");
    }

    NumericVector out_px(n);
    NumericVector out_py(n);

    for (int k = 0; k < n; k++) {
        // Project to the solid
        auto fwd = hexify::snyder_forward(lon[k], lat[k]);
        hexify::face_tri_to_plane(fwd.face, fwd.icosa_triangle_x,
                                  fwd.icosa_triangle_y, out_px[k], out_py[k]);
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}

