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
#include "cell_id.h"
#include "index_z7.h"
#include "globe_table.h"
#include "plane_clip.h"
#include "hex9.h"
#include "hex9_solid.h"

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
//
// Aperture 9 is Hex9 on the octahedron (hex9.h). Its level-L cells are stored
// as the aperture-3 Class II lattice of resolution 2L + 1 stores its cells,
// on a coset of that lattice off the solid's vertices, so there are no vertex
// quads of one cell; a cell ID is the cell's Hex9 address.
// ============================================================================

struct QuadFrame {
    std::vector<int> ap_seq;      // empty for a pure aperture
    int aperture;                 // 0 for a mixed sequence, 9 for Hex9
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

// Cells per diamond quad and in all, for the apertures of levels 1..resolution
// (`step(k)`), on the active solid. False where the count passes the largest
// int64, so that the grid's cell IDs cannot all be written.
template <typename Step>
static bool grid_cell_count(int resolution, Step step,
                            int64_t& per_quad, int64_t& n_cells) {
    const int64_t n_diamonds = hexify::topo().n_diamonds();
    per_quad = 1;
    for (int k = 1; k <= resolution; k++) {
        const int64_t a = step(k);
        if (per_quad > hexify::kCellIdMax / a) return false;
        per_quad *= a;
    }
    if (per_quad > (hexify::kCellIdMax - 2) / n_diamonds) return false;
    n_cells = n_diamonds * per_quad + 2;
    return true;
}

// Whether a frame is Hex9's
static inline bool frame_hex9(const QuadFrame& f) { return f.aperture == 9; }

// The finest resolution of a pure aperture whose cell IDs fit in an int64
static int pure_max_resolution(int aperture) {
    if (aperture == 9) return hexify::hex9::kMaxLevel;
    int64_t per_quad, n_cells;
    int r = hexify::kMinResolution;
    while (r < hexify::kMaxResolution &&
           grid_cell_count(r + 1, [&](int) { return int64_t(aperture); },
                           per_quad, n_cells)) {
        r++;
    }
    return r;
}

// The frame of a pure aperture at a resolution (empty 'ap_seq'), or of a mixed
// sequence, whose resolution is one less than its length.
static QuadFrame quad_frame(int resolution, int aperture, std::vector<int> ap_seq) {
    const bool mixed = !ap_seq.empty();
    if (mixed) {
        resolution = static_cast<int>(ap_seq.size()) - 1;
        aperture = 0;
    } else if (aperture != 3 && aperture != 4 && aperture != 7 && aperture != 9) {
        Rcpp::stop("aperture must be 3, 4, 7 or 9");
    }
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }

    QuadFrame f;
    if (aperture == 9) {
        if (hexify::topo().solid != hexify::Solid::Octahedron) {
            Rcpp::stop("aperture 9 (Hex9) is defined on the octahedron only");
        }
        if (resolution > hexify::hex9::kMaxLevel) {
            Rcpp::stop("aperture 9 at resolution %d has more cells than 64-bit cell "
                       "IDs can number (2^63 - 1); its finest resolution is %d",
                       resolution, hexify::hex9::kMaxLevel);
        }
        const int ap3_res = hexify::hex9::frame_resolution(resolution);
        f.aperture = 9;
        f.resolution = resolution;
        f.form = hexify::hex_form_pure(3, ap3_res);
        f.dim = hexify::quad_edge_dim(3, ap3_res);
        f.offsetPerQuad = static_cast<uint64_t>(f.dim * f.dim / 3);
        f.nCells = static_cast<uint64_t>(hexify::hex9::level_cells(resolution));
        f.lattice = sublattice_of(f.form);
        f.generator = generator_of(f.form);
        return f;
    }
    int64_t per_quad, n_cells;
    if (!grid_cell_count(resolution,
                         [&](int k) { return int64_t(mixed ? ap_seq[k] : aperture); },
                         per_quad, n_cells)) {
        if (mixed) {
            Rcpp::stop("this aperture sequence has more cells than 64-bit cell IDs "
                       "can number (2^63 - 1)");
        }
        Rcpp::stop("aperture %d at resolution %d has more cells than 64-bit cell "
                   "IDs can number (2^63 - 1); its finest resolution on this solid "
                   "is %d", aperture, resolution, pure_max_resolution(aperture));
    }
    f.aperture = aperture;
    f.resolution = resolution;
    f.form = mixed ? hexify::hex_form_sequence(ap_seq)
                   : hexify::hex_form_pure(aperture, resolution);
    f.dim = mixed ? hexify::quad_edge_dim(ap_seq)
                  : hexify::quad_edge_dim(aperture, resolution);
    f.offsetPerQuad = static_cast<uint64_t>(per_quad);
    f.nCells = static_cast<uint64_t>(n_cells);
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

// The Hex9 lattice point at stored (i, j) of a quad; stops when (i, j) is
// none, which a stored cell coordinate never is.
static inline hexify::hex9::OctPoint frame_hex9_point(const QuadFrame& f, int quad,
                                                      long long i, long long j) {
    hexify::hex9::OctPoint p;
    if (!hexify::hex9::quad_ij_lattice(quad, i, j, f.resolution, p)) {
        Rcpp::stop("hexify internal error: a Hex9 quad coordinate off the lattice");
    }
    return p;
}

// The ID, from 1, of the cell at stored (i, j) of a quad. Quad 0 holds the
// north pole alone, ID 1; quad q > 0 follows it and the q - 1 quads before.
// A Hex9 cell's ID is its address.
static inline int64_t frame_encode(const QuadFrame& f, int quad, long long i, long long j) {
    if (frame_hex9(f)) {
        const int64_t id = hexify::hex9::encode(frame_hex9_point(f, quad, i, j), f.resolution);
        if (id == 0) Rcpp::stop("hexify internal error: a Hex9 coordinate that is no cell");
        return id;
    }
    uint64_t offset = (quad == 0) ? 0 : 1 + static_cast<uint64_t>(quad - 1) * f.offsetPerQuad;
    uint64_t within_quad = (f.aperture == 7)
        ? hexify::ap7_surrogate_to_quad_index(i, j, f.resolution)
        : cell_index_2d(i, j, f.dim, f.lattice);
    return static_cast<int64_t>(offset + within_quad + 1);
}

// The 0-based index of a cell ID. Stops unless the ID names a cell of the
// grid, so no table is ever indexed with a quad, i or j read from NA or an
// out-of-range ID.
static inline uint64_t frame_cell_index(const QuadFrame& f, int64_t id) {
    if (id == hexify::kCellIdNA || id < 1 || static_cast<uint64_t>(id) > f.nCells) {
        Rcpp::stop("cell_id must be a whole number in [1, %lld] for resolution %d",
                   static_cast<long long>(f.nCells), f.resolution);
    }
    return static_cast<uint64_t>(id) - 1;
}

// The 0-based index of the ID in integer64 slot `slot`
static inline uint64_t frame_cell_index(const QuadFrame& f, double slot) {
    return frame_cell_index(f, hexify::cell_id_get(slot));
}

// A cell's quad and stored (i, j) from its 0-based index
static inline void frame_decode_index(const QuadFrame& f, uint64_t idx,
                                      int& quad, long long& i, long long& j) {
    if (frame_hex9(f)) {
        hexify::hex9::OctPoint c;
        if (!hexify::hex9::decode(static_cast<int64_t>(idx) + 1, f.resolution, c)) {
            Rcpp::stop("hexify internal error: a Hex9 ID that names no cell");
        }
        hexify::hex9::cell_quad_ij(c, f.resolution, quad, i, j);
        return;
    }
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
    if (frame_hex9(f)) {
        hexify::quad_ij_to_xy(quad, i, j, 3, hexify::hex9::frame_resolution(f.resolution), x, y);
    } else if (f.ap_seq.empty()) {
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
    if (frame_hex9(f)) return hexify::substrate_ij_canonicalize(quad, i, j, f.dim);
    if (f.ap_seq.empty()) {
        return hexify::quad_ij_canonicalize(quad, i, j, f.aperture, f.resolution);
    }
    return hexify::substrate_ij_canonicalize(quad, i, j, f.dim);
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
// the Eisenstein lattice divided by D = k * N(g), and its nearest lattice
// points are found exactly (nearest_eisenstein_points()).
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
    long long nearest[3][2];
    hexify::nearest_eisenstein_points(U, V, D, nearest);
    const long long best_x = nearest[0][0], best_y = nearest[0][1];

    i = best_x * a - best_y * b;
    j = best_x * b + best_y * a - best_y * b;
    if (!frame_in_quad(parent, i, j) && !frame_canonicalize(parent, quad, i, j)) {
        Rcpp::stop("hexify internal error: a parent cell centre has no owning quad");
    }
}

// The cell a projected point falls in: its quad and stored (i, j)
static inline void frame_locate_projected(const QuadFrame& f,
                                          const hexify::ProjectionResult& fwd,
                                          int& quad, long long& i, long long& j) {
    if (frame_hex9(f)) {
        const hexify::hex9::OctPoint c = hexify::hex9::face_point_cell(
            fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y, f.resolution);
        hexify::hex9::cell_quad_ij(c, f.resolution, quad, i, j);
        return;
    }
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

// The cell a lon/lat point falls in, latitude geodetic on the active ellipsoid
static inline void frame_locate(const QuadFrame& f, double lon_deg, double lat_deg,
                                int& quad, long long& i, long long& j) {
    frame_locate_projected(f, hexify::snyder_forward(lon_deg, lat_deg), quad, i, j);
}

// Substrate steps along a quad edge of a grid
// [[Rcpp::export]]
double cpp_quad_edge_dim(int resolution, int aperture, IntegerVector ap_seq) {
    return static_cast<double>(grid_frame(resolution, aperture, ap_seq).dim);
}

// The number of cells of a grid, as integer64; NA where it passes the largest
// int64, so that its cell IDs cannot all be written
// [[Rcpp::export]]
NumericVector cpp_grid_n_cells(NumericVector icosa, int resolution, int aperture,
                               IntegerVector ap_seq) {
    activate_grid(icosa);
    const bool mixed = ap_seq.size() > 0;
    if (mixed) resolution = static_cast<int>(ap_seq.size()) - 1;
    if (!mixed && aperture == 9) {
        if (resolution < 0 || resolution > hexify::hex9::kMaxLevel) return hexify::cell_id_na(1);
        return hexify::cell_id_vector(
            std::vector<int64_t>{hexify::hex9::level_cells(resolution)});
    }
    int64_t per_quad, n_cells;
    if (!grid_cell_count(resolution,
                         [&](int k) { return int64_t(mixed ? ap_seq[k] : aperture); },
                         per_quad, n_cells)) {
        return hexify::cell_id_na(1);
    }
    return hexify::cell_id_vector(std::vector<int64_t>{n_cells});
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
    NumericVector result = hexify::cell_id_na(n);

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
            continue;
        }
        result[k] = hexify::cell_id_slot(frame_encode(f, q, ii, jj));
    }

    return result;
}

// ============================================================================
// IGEO7 integer forms
// ============================================================================
// An aperture-7 cell on the icosahedron as IGEO7 writes it beside the Z7
// string: the 64-bit packed index, as 16 hexadecimal digits (DGGRID's INT64
// output), and the monotonic ID. The packed index's 64 bits travel through R
// in an integer64 slot unchanged. integer64 is signed, so the indices of base
// cells 8-11, whose top bit is set, read as negative numbers there, and the
// one index whose bits are 2^63, base cell 8's pentagon at resolution 20, is
// the bit pattern bit64 reserves for NA.

enum class Z7Form { Packed, Hex, Monotonic };

static Z7Form z7_form(const std::string& form) {
    if (form == "int") return Z7Form::Packed;
    if (form == "hex") return Z7Form::Hex;
    if (form == "monotonic") return Z7Form::Monotonic;
    Rcpp::stop("form must be \"int\", \"hex\" or \"monotonic\"");
}

static void require_igeo7_frame(const QuadFrame& f) {
    if (f.aperture != 7 || !hexify::z7::igeo7_labels()) {
        Rcpp::stop("IGEO7's integer forms are defined for aperture-7 grids on the "
                   "icosahedron");
    }
}

// The IGEO7 label of the aperture-7 cell at a frame's stored (i, j)
static inline hexify::z7::Label frame_z7_label(const QuadFrame& f, int quad,
                                               long long i, long long j) {
    long long si, sj;
    hexify::ap7_surrogate_to_substrate_ijk(i, j, f.resolution, si, sj);
    return hexify::z7::label_of(quad, si, sj, f.resolution);
}

// The ID of the aperture-7 cell a label names
static inline int64_t frame_z7_cell(const QuadFrame& f, const hexify::z7::Label& label) {
    int quad;
    long long si, sj, i, j;
    hexify::z7::cell_of(label, quad, si, sj);
    hexify::ap7_substrate_to_surrogate_ijk(si, sj, f.resolution, i, j);
    return frame_encode(f, quad, i, j);
}

static std::string packed_hex(uint64_t z) {
    static const char digits[] = "0123456789abcdef";
    std::string out(16, '0');
    for (int k = 15; k >= 0; k--) {
        out[k] = digits[z & 15u];
        z >>= 4;
    }
    return out;
}

// The packed index 1 to 16 hexadecimal digits spell, with or without "0x"
static bool parse_packed_hex(const std::string& s, uint64_t& z) {
    size_t start = (s.size() > 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) ? 2 : 0;
    if (s.size() == start || s.size() - start > 16) return false;
    z = 0;
    for (size_t k = start; k < s.size(); k++) {
        const char c = s[k];
        int v;
        if (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'a' && c <= 'f') v = 10 + c - 'a';
        else if (c >= 'A' && c <= 'F') v = 10 + c - 'A';
        else return false;
        z = (z << 4) | static_cast<uint64_t>(v);
    }
    return true;
}

// IGEO7 integer form of each cell: integer64 for "int" and "monotonic",
// character for "hex"
// [[Rcpp::export]]
SEXP cpp_cell_to_z7(NumericVector icosa, NumericVector cell_id, int resolution,
                    std::string form) {
    activate_grid(icosa);
    const QuadFrame f = grid_frame(resolution, 7, IntegerVector(0));
    require_igeo7_frame(f);
    const Z7Form fm = z7_form(form);
    if (fm != Z7Form::Monotonic && resolution > hexify::z7::kMaxPackedRes) {
        Rcpp::stop("IGEO7's packed index holds resolutions 0 to %d",
                   hexify::z7::kMaxPackedRes);
    }
    hexify::require_cell_ids(cell_id);
    const R_xlen_t n = cell_id.size();
    NumericVector ids = hexify::cell_id_na(n);
    CharacterVector hex(fm == Z7Form::Hex ? n : 0);

    for (R_xlen_t k = 0; k < n; k++) {
        if (hexify::cell_id_get(cell_id[k]) == hexify::kCellIdNA) {
            if (fm == Z7Form::Hex) hex[k] = NA_STRING;
            continue;
        }
        int quad;
        long long i, j;
        frame_decode(f, cell_id[k], quad, i, j);
        const hexify::z7::Label label = frame_z7_label(f, quad, i, j);
        if (fm == Z7Form::Monotonic) {
            ids[k] = hexify::cell_id_slot(static_cast<int64_t>(hexify::z7::to_monotonic(label)));
        } else if (fm == Z7Form::Packed) {
            ids[k] = hexify::cell_id_slot(static_cast<int64_t>(hexify::z7::to_packed(label)));
        } else {
            hex[k] = packed_hex(hexify::z7::to_packed(label));
        }
    }
    if (fm == Z7Form::Hex) return hex;
    return ids;
}

// The cells of IGEO7 integer forms at a grid's resolution: `index` is
// integer64 for "int" and "monotonic" and character for "hex". An index of
// another resolution, or bits that spell no index, stop with its position.
// [[Rcpp::export]]
NumericVector cpp_z7_to_cell(NumericVector icosa, SEXP index, int resolution,
                             std::string form) {
    activate_grid(icosa);
    const QuadFrame f = grid_frame(resolution, 7, IntegerVector(0));
    require_igeo7_frame(f);
    const Z7Form fm = z7_form(form);
    const bool text = (fm == Z7Form::Hex);
    if (text != (TYPEOF(index) == STRSXP)) {
        Rcpp::stop(text ? "hexadecimal indices must be character"
                        : "integer indices must be integer64");
    }
    const R_xlen_t n = Rf_xlength(index);
    NumericVector out = hexify::cell_id_na(n);
    NumericVector values = text ? NumericVector(0) : NumericVector(index);
    if (!text) hexify::require_cell_ids(values);
    CharacterVector strings = text ? CharacterVector(index) : CharacterVector(0);

    for (R_xlen_t k = 0; k < n; k++) {
        uint64_t z;
        if (text) {
            if (CharacterVector::is_na(strings[k])) continue;
            if (!parse_packed_hex(Rcpp::as<std::string>(strings[k]), z)) {
                Rcpp::stop("index %lld is not 1 to 16 hexadecimal digits",
                           static_cast<long long>(k + 1));
            }
        } else {
            const int64_t v = hexify::cell_id_get(values[k]);
            if (v == hexify::kCellIdNA) continue;
            z = static_cast<uint64_t>(v);
        }
        hexify::z7::Label label;
        const bool ok = (fm == Z7Form::Monotonic)
            ? hexify::z7::from_monotonic(z, resolution, label)
            : hexify::z7::from_packed(z, label);
        if (!ok) {
            Rcpp::stop("index %lld is not an IGEO7 %s", static_cast<long long>(k + 1),
                       fm == Z7Form::Monotonic ? "monotonic ID at this resolution"
                                               : "packed index");
        }
        if (label.res != resolution) {
            Rcpp::stop("index %lld is a resolution-%d cell, not one of this "
                       "resolution-%d grid", static_cast<long long>(k + 1),
                       label.res, resolution);
        }
        out[k] = hexify::cell_id_slot(frame_z7_cell(f, label));
    }
    return out;
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

    hexify::require_cell_ids(cell_id);
    R_xlen_t n = cell_id.size();
    NumericVector result = hexify::cell_id_na(n);
    for (R_xlen_t k = 0; k < n; k++) {
        if (hexify::cell_id_get(cell_id[k]) == hexify::kCellIdNA) continue;
        int quad;
        long long i, j;
        frame_decode(child, cell_id[k], quad, i, j);
        frame_parent(child, parent, quad, i, j);
        result[k] = hexify::cell_id_slot(frame_encode(parent, quad, i, j));
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
    NumericVector result = hexify::cell_id_na(n);

    for (R_xlen_t k = 0; k < n; k++) {
        if (ISNAN(lon[k]) || ISNAN(lat[k])) continue;
        int quad;
        long long i, j;
        frame_locate(f, lon[k], lat[k], quad, i, j);
        result[k] = hexify::cell_id_slot(frame_encode(f, quad, i, j));
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_lonlat(NumericVector icosa, NumericVector cell_id,
                             int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    hexify::require_cell_ids(cell_id);
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
    hexify::require_cell_ids(cell_id);
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
    hexify::require_cell_ids(cell_id);
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
    NumericVector result = hexify::cell_id_na(n);

    for (R_xlen_t k = 0; k < n; k++) {
        int q = quad[k];
        double qx = quad_x[k];
        double qy = quad_y[k];

        int out_quad;
        long long i, j;
        if (aperture == 7) {
            // AP7: exact-integer quantization straight to the surrogate.
            hexify::quad_xy_to_ij(q, qx, qy, 7, resolution, out_quad, i, j);
        } else if (aperture == 9) {
            int face;
            double tx, ty;
            hexify::quad_xy_to_icosa_tri(q, qx, qy, face, tx, ty);
            hexify::hex9::cell_quad_ij(hexify::hex9::face_point_cell(face, tx, ty, resolution),
                                       resolution, out_quad, i, j);
        } else {
            // AP3/AP4: through the face the point lies on, which names its quad
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            hexify::quad_xy_to_icosa_tri(q, qx, qy, icosa_triangle_face,
                                         icosa_triangle_x, icosa_triangle_y);
            hexify::icosa_tri_to_quad_ij(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y,
                                         aperture, resolution, out_quad, i, j);
        }
        result[k] = hexify::cell_id_slot(frame_encode(f, out_quad, i, j));
    }

    return result;
}

// Points of the quad planes in lon/lat, latitude geodetic on the active
// ellipsoid; NA where a point has no image (past the quad's fold at a
// vertex).
// [[Rcpp::export]]
DataFrame cpp_quad_xy_to_lonlat(NumericVector icosa, IntegerVector quad,
                                NumericVector quad_x, NumericVector quad_y) {
    activate_grid(icosa);
    const R_xlen_t n = quad.size();
    if (quad_x.size() != n || quad_y.size() != n) {
        stop("quad, quad_x and quad_y must have the same length");
    }
    NumericVector lon(n, NA_REAL), lat(n, NA_REAL);
    for (R_xlen_t k = 0; k < n; k++) {
        int face;
        double tx, ty;
        if (quad[k] == NA_INTEGER ||
            !hexify::try_quad_xy_to_icosa_tri(quad[k], quad_x[k], quad_y[k], face, tx, ty)) {
            continue;
        }
        const auto ll = hexify::face_xy_to_ll(tx, ty, face);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }
    return DataFrame::create(_["lon_deg"] = lon, _["lat_deg"] = lat);
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
    hexify::require_cell_ids(cell_id);
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

// Whether an edge a -> b, whose true midpoint is m, must be split before the
// great-circle arc a -> b may stand for it: m lies further off the arc's
// great circle than 'tolerance' times the arc. A lon/lat polygon is read with
// great-circle edges on the sphere (s2, sf's default for longlat data), and
// such an arc departs from a straight lon/lat chord most where the edge runs
// east-west, so a piece can pass chord_needs_split() and not this.
static inline bool arc_needs_split(double alon, double alat,
                                   double mlon, double mlat,
                                   double blon, double blat,
                                   double tolerance) {
    auto unit = [](double lon, double lat, double v[3]) {
        const double la = lat * hexify::kDegToRad, lo = lon * hexify::kDegToRad;
        v[0] = std::cos(la) * std::cos(lo);
        v[1] = std::cos(la) * std::sin(lo);
        v[2] = std::sin(la);
    };
    double a[3], m[3], b[3];
    unit(alon, alat, a);
    unit(mlon, mlat, m);
    unit(blon, blat, b);
    const double n[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                         a[0] * b[1] - a[1] * b[0]};
    const double sin_arc = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (sin_arc == 0.0) return false;
    const double arc = std::atan2(sin_arc, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    // sin of m's angular distance from the great circle through a and b
    const double off = std::fabs(m[0] * n[0] + m[1] * n[1] + m[2] * n[2]) / sin_arc;
    return off > std::sin(tolerance * arc);
}

// Whether a piece a -> b, whose true midpoint is m, at halving 'depth', must
// be halved again: it strays from the edge read as a lon/lat chord
// (chord_needs_split) or as a great-circle arc (arc_needs_split), or it spans
// more than the longest arc allowed.
static inline bool piece_needs_split(double alon, double alat,
                                     double mlon, double mlat,
                                     double blon, double blat,
                                     const EdgeLimits& lim, int depth) {
    if (lim.tolerance > 0.0 && depth < kMaxEdgeSplits &&
        (chord_needs_split(alon, alat, mlon, mlat, blon, blat, lim.tolerance) ||
         arc_needs_split(alon, alat, mlon, mlat, blon, blat, lim.tolerance))) {
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
static bool quad_point_face(int quad, double qx, double qy,
                            int& face, double& tx, double& ty);

static void quad_point_lonlat(int quad, double qx, double qy,
                              double qx_center, double qy_center,
                              double& lon, double& lat) {
    lon = NA_REAL;
    lat = NA_REAL;

    int face;
    double tx, ty;
    if (quad_point_face(quad, qx, qy, face, tx, ty)) {
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
// which lies in no face. A point below a near edge is read in the fan of
// faces around the quad's origin, or, on a ray bounding the solid's deficit
// there (where the corner of a Hex9 cell beside a vertex lies), through the
// quad across that edge.
static bool quad_point_face(int quad, double qx, double qy,
                            int& face, double& tx, double& ty) {
    int q = quad;
    double x = qx, y = qy;
    if (hexify::quad_xy_canonicalize(q, x, y) &&
        hexify::try_quad_xy_to_icosa_tri(q, x, y, face, tx, ty)) {
        return true;
    }
    q = quad;
    x = qx;
    y = qy;
    return hexify::quad_xy_canonicalize(q, x, y, /*across_near_edges=*/true) &&
           hexify::try_quad_xy_to_icosa_tri(q, x, y, face, tx, ty);
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

// The face triangle, the same in every face's triangle coordinates (unit
// edge, first corner at the top): its corners, then its edge midpoints.
static const double kFaceTrianglePoints[6][2] = {
    {0.5, 0.86602540378443864676}, {0.0, 0.0}, {1.0, 0.0},
    {0.25, 0.43301270189221932338}, {0.5, 0.0}, {0.75, 0.43301270189221932338}};

// The most creases a piece inside one face can cross: the lines from the
// face centre point in six directions 60 degrees apart, and a line misses
// the centre, so it sees them over less than 180 degrees.
constexpr int kMaxCreaseCrossings = 3;

// The face projection has a crease along lines from a face's centre: its
// derivative jumps there, so the image of a straight piece bends where it
// crosses one. Snyder's has one to each corner; the vertex-oriented
// projection (IVEA) one to each corner and one to each edge midpoint;
// Fuller's is smooth inside a face and has none. The face triangle is the
// same in every face's triangle coordinates, with unit edge and its first
// corner at the top. Returns how many creases the piece a -> b crosses
// strictly between its ends, with the crossings' places along it, ascending,
// in 'u'.
static int radius_crossings(double ax, double ay, double bx, double by,
                            double u[kMaxCreaseCrossings]) {
    constexpr double kCx = 0.5;
    constexpr double kCy = 0.28867513459481288225;   // 1 / (2 sqrt(3))
    // The far end of each crease: the corners, then the edge midpoints
    const double (*ends)[2] = kFaceTrianglePoints;
    int n_lines = 0;
    switch (hexify::active_projection()) {
        case hexify::FaceProjection::ISEA: n_lines = 3; break;
        case hexify::FaceProjection::IVEA: n_lines = 6; break;
        case hexify::FaceProjection::Fuller: n_lines = 0; break;
    }
    // A crossing this close to an end is that end.
    constexpr double kEndSlack = 1e-12;
    const double dx = bx - ax, dy = by - ay;
    int n = 0;
    for (int k = 0; k < n_lines && n < kMaxCreaseCrossings; k++) {
        const double ex = ends[k][0] - kCx, ey = ends[k][1] - kCy;
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
    double u[kMaxCreaseCrossings];
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
// inside a cell, on a cell edge, or -- in a Hex9 grid, whose cells meet at the
// solid's vertices -- on a corner, which pole_corners() splits.

// A point within this many degrees of latitude of a pole is on it.
constexpr double kPoleCornerSlack = 1e-9;

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
            // An end on a pole has no longitude of its own; the piece leaves
            // it along the meridian of its other end.
            const double alon_m = (std::fabs(std::fabs(alat) - 90.0) < kPoleCornerSlack) ? mlon : alon;
            const double blon_m = (std::fabs(std::fabs(blat) - 90.0) < kPoleCornerSlack) ? mlon : blon;
            if (R_finite(mlon) &&
                piece_needs_split(alon_m, alat, mlon, mlat, blon_m, blat, lim, depth)) {
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

// A ring point on a pole has no longitude of its own. The ring's edges reach
// it along two meridians, so it is written as two points on the pole at the
// longitudes of the points before and after it, and the lon/lat ring runs
// along both.
static void pole_corners(std::vector<double>& lon, std::vector<double>& lat) {
    const size_t n = lon.size();
    if (n < 3) return;
    bool any = false;
    for (size_t k = 0; k < n; k++) {
        if (std::fabs(std::fabs(lat[k]) - 90.0) < kPoleCornerSlack) any = true;
    }
    if (!any) return;
    std::vector<double> out_lon, out_lat;
    out_lon.reserve(n + 2);
    out_lat.reserve(n + 2);
    for (size_t k = 0; k < n; k++) {
        if (std::fabs(std::fabs(lat[k]) - 90.0) >= kPoleCornerSlack) {
            out_lon.push_back(lon[k]);
            out_lat.push_back(lat[k]);
            continue;
        }
        const double pole = lat[k] > 0.0 ? 90.0 : -90.0;
        out_lon.push_back(lon[(k + n - 1) % n]);
        out_lat.push_back(pole);
        out_lon.push_back(lon[(k + 1) % n]);
        out_lat.push_back(pole);
    }
    lon.swap(out_lon);
    lat.swap(out_lat);
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
    hexify::require_cell_ids(cell_id);
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
        pole_corners(lon, lat);
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

// A cell's part on one face of the solid, a convex polygon in that face's
// triangle coordinates.
struct FacePolygon {
    int face;
    hexify::PlanePolygon poly;
};

// How close, in radians, a face corner lies to a vertex cell's centre when it
// is the vertex the cell is centred on.
constexpr double kVertexCornerSlack = 1e-9;

// The parts of a cell on the faces it covers. The cell is convex on the
// unfolded solid and each face a triangle of it, so each part is convex: the
// hull of the cell's boundary pieces on that face, with the face's corner at
// the vertex for a cell centred on a vertex of the solid.
static void cell_face_polygons(const CellPlanes& g, const CellPlane& c,
                               std::vector<PlaneEdge>& edges,
                               std::vector<FacePiece>& pieces,
                               std::vector<FacePolygon>& out) {
    cell_face_pieces(g, c, edges, pieces);
    out.clear();
    for (const FacePiece& p : pieces) {
        auto it = std::find_if(out.begin(), out.end(),
                               [&](const FacePolygon& fp) { return fp.face == p.face; });
        if (it == out.end()) {
            out.push_back({p.face, {}});
            it = out.end() - 1;
        }
        it->poly.push_back({p.ax, p.ay});
        it->poly.push_back({p.bx, p.by});
    }
    if (c.at_vertex) {
        const hexify::UnitVec centre = cell_centre_sphere(c);
        hexify::UnitVec v;
        for (FacePolygon& fp : out) {
            for (int k = 0; k < 3; k++) {
                hexify::face_tri_to_sphere(fp.face, kFaceTrianglePoints[k][0],
                                           kFaceTrianglePoints[k][1], v.data());
                if (hexify::arc_angle(v, centre) < kVertexCornerSlack) {
                    fp.poly.push_back({kFaceTrianglePoints[k][0], kFaceTrianglePoints[k][1]});
                }
            }
        }
    }
    for (FacePolygon& fp : out) fp.poly = hexify::convex_hull(fp.poly);
    out.erase(std::remove_if(out.begin(), out.end(),
                             [](const FacePolygon& fp) { return fp.poly.size() < 3; }),
              out.end());
}

// Subdivisions of each smooth edge piece at the coarsest level of
// face_polygon_solid_angle(), and how many times it doubles.
constexpr int kAreaBaseSteps = 8;
constexpr int kAreaLevels = 4;

// Solid angle of a convex polygon of one face. Its edges are cut where they
// cross a crease of the projection, so each piece maps to a smooth curve. The
// polygon through n equally spaced points of every piece, fanned from an
// inner point, falls short of the region by a2 / n^2 + a4 / n^4 + ..., the
// error of the midpoint rule over each piece, so Romberg's extrapolation over
// n = 8, 16, 32, 64 removes the terms through n^-6.
static double face_polygon_solid_angle(int face, const hexify::PlanePolygon& poly,
                                       std::vector<FacePiece>& pieces,
                                       std::vector<hexify::UnitVec>& pts) {
    pieces.clear();
    for (size_t i = 0; i < poly.size(); i++) {
        const hexify::PlanePoint& a = poly[i];
        const hexify::PlanePoint& b = poly[(i + 1) % poly.size()];
        pieces.push_back({face, a.x, a.y, b.x, b.y, static_cast<int>(i), 0.0, 1.0});
    }
    split_at_creases(pieces);

    // Every piece at the finest level, its last point left to the next piece
    const int n_fine = kAreaBaseSteps << (kAreaLevels - 1);
    pts.clear();
    hexify::UnitVec v;
    for (const FacePiece& p : pieces) {
        for (int k = 0; k < n_fine; k++) {
            const double t = static_cast<double>(k) / n_fine;
            hexify::face_tri_to_sphere(face, p.ax + t * (p.bx - p.ax),
                                       p.ay + t * (p.by - p.ay), v.data());
            pts.push_back(v);
        }
    }
    hexify::UnitVec o = {0.0, 0.0, 0.0};
    for (size_t k = 0; k < pts.size(); k += n_fine) {
        for (int d = 0; d < 3; d++) o[d] += pts[k][d];
    }
    o = hexify::normalized(o);

    double r[kAreaLevels];
    for (int level = 0; level < kAreaLevels; level++) {
        const size_t stride = static_cast<size_t>(n_fine / (kAreaBaseSteps << level));
        double omega = 0.0;
        for (size_t k = 0; k < pts.size(); k += stride) {
            omega += hexify::triangle_solid_angle(o, pts[k], pts[(k + stride) % pts.size()]);
        }
        r[level] = omega;
    }
    for (int m = 1; m < kAreaLevels; m++) {
        const double f = std::pow(4.0, m);
        for (int level = kAreaLevels - 1; level >= m; level--) {
            r[level] = (f * r[level] - r[level - 1]) / (f - 1.0);
        }
    }
    return std::fabs(r[kAreaLevels - 1]);
}

// Solid angles of cells and of their parts inside coarser cells. Both grids'
// cells are straight on each face, so a cell and a coarser cell are clipped
// against each other face by face (cell_face_polygons, clip_convex) and each
// piece measured on the sphere (face_polygon_solid_angle). 'pair_cell' and
// 'pair_parent' (from 1, grouped by cell) name the pairs; returns 'piece', one
// solid angle per pair, and 'whole', one per cell that appears in a pair.
// [[Rcpp::export]]
List cpp_cell_overlap_solid_angles(NumericVector icosa, NumericVector cell_id,
                                   int resolution, int aperture, IntegerVector ap_seq,
                                   NumericVector parent_id, int parent_resolution,
                                   int parent_aperture, IntegerVector parent_ap_seq,
                                   IntegerVector pair_cell, IntegerVector pair_parent) {
    activate_grid(icosa);
    if (pair_cell.size() != pair_parent.size()) stop("pair_cell and pair_parent differ in length");
    const CellPlanes gc = cell_planes(cell_id, grid_frame(resolution, aperture, ap_seq));
    const CellPlanes gp = cell_planes(parent_id, grid_frame(parent_resolution, parent_aperture,
                                                            parent_ap_seq));
    std::vector<PlaneEdge> edges;
    std::vector<FacePiece> pieces;
    std::vector<hexify::UnitVec> pts;

    std::vector<std::vector<FacePolygon>> parents(parent_id.size());
    for (R_xlen_t k = 0; k < parent_id.size(); k++) {
        cell_face_polygons(gp, gp.cells[k], edges, pieces, parents[k]);
    }

    NumericVector piece(pair_cell.size()), whole(cell_id.size());
    std::vector<FacePolygon> cell;
    R_xlen_t current = -1;
    for (R_xlen_t m = 0; m < pair_cell.size(); m++) {
        const R_xlen_t k = pair_cell[m] - 1, p = pair_parent[m] - 1;
        if (k < 0 || k >= cell_id.size() || p < 0 || p >= parent_id.size()) {
            stop("pair_cell or pair_parent out of range");
        }
        if (k != current) {
            current = k;
            cell_face_polygons(gc, gc.cells[k], edges, pieces, cell);
            double w = 0.0;
            for (const FacePolygon& fp : cell) {
                w += face_polygon_solid_angle(fp.face, fp.poly, pieces, pts);
            }
            whole[k] = w;
        }
        double s = 0.0;
        for (const FacePolygon& cf : cell) {
            for (const FacePolygon& pf : parents[p]) {
                if (cf.face != pf.face) continue;
                const hexify::PlanePolygon cut = hexify::clip_convex(cf.poly, pf.poly);
                if (!cut.empty()) {
                    s += face_polygon_solid_angle(cf.face, cut, pieces, pts);
                }
            }
        }
        piece[m] = s;
    }
    return List::create(_["piece"] = piece, _["whole"] = whole);
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

hexify::GlobeFrame hexify::globe_frame(int resolution, int aperture,
                                       const std::vector<int>& ap_seq) {
    const QuadFrame f = quad_frame(resolution, aperture, ap_seq);
    const SubstrateLattice lattice = sublattice_of(f.form);
    const LatticeGenerator generator = generator_of(f.form);
    return {f.dim, lattice.index, lattice.c, generator.a, generator.b,
            f.offsetPerQuad, f.nCells};
}

// What a renderer needs to find the cell of a quad-plane point by itself:
// the grid's GlobeFrame.
// [[Rcpp::export]]
List cpp_globe_frame(NumericVector icosa, int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    const QuadFrame q = grid_frame(resolution, aperture, ap_seq);
    const hexify::GlobeFrame f = hexify::globe_frame(q.resolution, q.aperture, q.ap_seq);
    return List::create(
        _["dim"] = static_cast<double>(f.dim),
        _["index"] = static_cast<double>(f.index),
        _["c"] = static_cast<double>(f.c),
        _["generator"] = NumericVector::create(static_cast<double>(f.ga),
                                               static_cast<double>(f.gb)),
        _["per_quad"] = hexify::cell_id_vector(
            std::vector<int64_t>{static_cast<int64_t>(f.per_quad)}),
        _["n_cells"] = hexify::cell_id_vector(
            std::vector<int64_t>{static_cast<int64_t>(f.n_cells)}));
}

// The cells adjacent to a vertex quad's cell. The vertex is a corner of every
// diamond quad around it, and each of its neighbours lies one lattice step
// from that corner inside one of those quads' boxes, so the neighbours are the
// steps from the corner that land in the box, over every quad the vertex is a
// corner of. No step leaves a box, so none is read across a face edge.
static void pole_neighbors(const QuadFrame& f, int pole, const long long offsets[6][2],
                           std::vector<int64_t>& out) {
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

// How a neighbour step that leaves its quad finds the cell it reaches.
enum class Crossing {
    Z7Digits,  // aperture 7 on the icosahedron: IGEO7 digit arithmetic
    LonLat     // the centre sent back through the forward pipeline
};

// The neighbours of an aperture-7 cell on the icosahedron by IGEO7 digit
// arithmetic: the cell's label stepped in the six directions of its base
// cell's frame, a pentagon's deleted direction giving none.
static void z7_neighbors(const QuadFrame& f, int quad, long long i, long long j,
                         std::vector<int64_t>& out) {
    const hexify::z7::Label label = frame_z7_label(f, quad, i, j);
    hexify::z7::Label nbr;
    for (int d = 1; d <= 6; d++) {
        if (hexify::z7::neighbor(label, static_cast<hexify::z7::Digit>(d), nbr)) {
            out.push_back(frame_z7_cell(f, nbr));
        }
    }
}

// The six cells adjacent to each of `cell_id`, in the frame's own grid.
//
// The six neighbours are the generator times the six units of the Eisenstein
// integers, so one step table serves every lattice. A cell with a step leaving
// its quad takes all six from IGEO7 digit arithmetic on an aperture-7 grid of
// the icosahedron, which is exact across quad and face edges; on other grids
// such a step's centre is sent back through the forward pipeline, which names
// the quad that owns it.
static Rcpp::List neighbors_in_frame(const Rcpp::NumericVector& cell_id,
                                      const QuadFrame& f,
                                      Crossing crossing = Crossing::Z7Digits) {
    hexify::require_cell_ids(cell_id);
    int n = cell_id.size();
    Rcpp::List out(n);

    long long offsets[6][2];
    lattice_unit_steps(f.generator, offsets);
    const bool z7_crossing = crossing == Crossing::Z7Digits && f.aperture == 7 &&
                             hexify::z7::igeo7_labels();

    for (int k = 0; k < n; k++) {
        uint64_t idx = frame_cell_index(f, cell_id[k]);

        std::vector<int64_t> neighbor_ids;
        neighbor_ids.reserve(6);

        // Resolution 0 is one base cell per vertex of the solid, a vertex cell
        // each, and each quad holds a single cell, so adjacency there is the
        // solid's vertex graph rather than a step through a quad frame.
        if (f.resolution == 0 && !frame_hex9(f)) {
            for (int v : hexify::topo().neighbors[idx]) neighbor_ids.push_back(v + 1);
            out[k] = hexify::cell_id_vector(neighbor_ids);
            continue;
        }

        int quad;
        long long i, j;

        if (!frame_hex9(f) && (idx == 0 || idx == f.nCells - 1)) {
            // The two vertex quads hold a single cell each -- a vertex of the
            // solid where several quads meet -- so their own frame carries no
            // offsets to step through.
            pole_neighbors(f, idx == 0 ? 0 : hexify::topo().south_pole(), offsets,
                           neighbor_ids);
            std::sort(neighbor_ids.begin(), neighbor_ids.end());
            neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()),
                               neighbor_ids.end());
            out[k] = hexify::cell_id_vector(neighbor_ids);
            continue;
        }
        frame_decode_index(f, idx, quad, i, j);

        long long step_i[6], step_j[6];
        bool inside[6];
        bool all_inside = true;
        for (int d = 0; d < 6; d++) {
            step_i[d] = i + offsets[d][0];
            step_j[d] = j + offsets[d][1];
            inside[d] = frame_in_quad(f, step_i[d], step_j[d]);
            all_inside = all_inside && inside[d];
        }

        if (!all_inside && z7_crossing) {
            z7_neighbors(f, quad, i, j, neighbor_ids);
            std::sort(neighbor_ids.begin(), neighbor_ids.end());
            neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()),
                               neighbor_ids.end());
            out[k] = hexify::cell_id_vector(neighbor_ids);
            continue;
        }

        for (int d = 0; d < 6; d++) {
            long long ni = step_i[d];
            long long nj = step_j[d];

            if (inside[d]) {
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

            auto ll = hexify::face_xy_to_sphere_ll(tri_x, tri_y, tri_face);
            int final_quad;
            long long final_i, final_j;
            frame_locate_projected(f, hexify::snyder_forward_sphere(ll.first, ll.second),
                                   final_quad, final_i, final_j);

            neighbor_ids.push_back(frame_encode(f, final_quad, final_i, final_j));
        }

        // Remove duplicates (can happen at pentagons / boundary)
        std::sort(neighbor_ids.begin(), neighbor_ids.end());
        neighbor_ids.erase(std::unique(neighbor_ids.begin(), neighbor_ids.end()),
                           neighbor_ids.end());
        // Remove self
        const int64_t self_id = static_cast<int64_t>(idx) + 1;
        neighbor_ids.erase(
            std::remove(neighbor_ids.begin(), neighbor_ids.end(), self_id),
            neighbor_ids.end());

        out[k] = hexify::cell_id_vector(neighbor_ids);
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_get_neighbors_isea(NumericVector icosa, Rcpp::NumericVector cell_id,
                                  int resolution, int aperture, IntegerVector ap_seq) {
    activate_grid(icosa);
    return neighbors_in_frame(cell_id, grid_frame(resolution, aperture, ap_seq));
}

// cpp_get_neighbors_isea() with every quad-leaving step sent through the
// forward pipeline, for checking the digit arithmetic against it
// [[Rcpp::export]]
Rcpp::List cpp_get_neighbors_isea_lonlat(NumericVector icosa, Rcpp::NumericVector cell_id,
                                         int resolution, int aperture,
                                         IntegerVector ap_seq) {
    activate_grid(icosa);
    return neighbors_in_frame(cell_id, grid_frame(resolution, aperture, ap_seq),
                              Crossing::LonLat);
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
    std::vector<int64_t> row_nbr;
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
            stop("cell %lld: its %d walls do not pair one to one with its %d neighbours",
                 static_cast<long long>(hexify::cell_id_get(cell_id[k])),
                 static_cast<int>(shapes.size()), static_cast<int>(nb.size()));
        }
        const hexify::UnitVec centre = cell_centre_sphere(g.cells[k]);
        for (R_xlen_t j = 0; j < nb.size(); j++) {
            rows.add(k, hexify::measure_wall(shapes[wall_of[j]], centre, nbr_centres[j]));
            row_nbr.push_back(hexify::cell_id_get(nb[j]));
        }
    }

    List out = List::create(_["perimeter"] = perimeter);
    if (walls) {
        out["walls"] = rows.frame("neighbor_id", hexify::cell_id_vector(row_nbr));
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
    hexify::require_cell_ids(cell_id);
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

