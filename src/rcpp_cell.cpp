// rcpp_cell.cpp
// Rcpp bindings for cell ID conversion and coordinate transforms
//
// This file provides the R interface for:
// - Lon/lat to cell ID conversion
// - Cell ID to lon/lat conversion
// - Cell ID to cell info conversion
// - Quad IJ coordinate conversion
// - Z7 decoding
// - PLANE coordinate conversions
//
// Copyright (c) 2024-2025 hexify authors. MIT License.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include "constants.h"
#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "aperture.h"
#include "aperture_sequence.h"
#include "grid_math.h"
#include "index_z7.h"
#include "ijk_coordinates.h"
#include "coordinate_transforms.h"

using namespace Rcpp;


// ============================================================================
// Z7 Decoding
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_decode_z7(std::string index_body, int aperture) {
    if (aperture != 7) {
        Rcpp::stop("cpp_decode_z7: only aperture 7 is supported");
    }

    try {
        int quadNum;
        long long i, j;

        // Calculate resolution from the Z7 string
        int resolution = index_body.length() - 2;
        if (resolution < 0) {
            resolution = 0;
        }

        // Call Z7 decode implementation
        hexify::z7::decode(index_body, resolution, quadNum, i, j);

        return DataFrame::create(
            Named("quad") = quadNum,
            Named("i") = i,
            Named("j") = j,
            Named("resolution") = resolution
        );

    } catch (const std::exception& e) {
        Rcpp::stop("Error in cpp_decode_z7: %s", e.what());
    }
}

// ============================================================================
// Triangle to Quad Coordinate Conversion
// ============================================================================

// [[Rcpp::export]]
Rcpp::List cpp_icosa_tri_to_quad_ij(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                                     int aperture, int resolution) {
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
Rcpp::List cpp_icosa_tri_to_quad_xy(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y) {
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
Rcpp::List cpp_quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y) {
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
Rcpp::List cpp_lonlat_to_quad_ij(double lon_deg, double lat_deg,
                                  int aperture, int resolution) {
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
// Quad coordinates to Cell ID
// ============================================================================

// ============================================================================
// Grid Pattern Helpers
// ============================================================================
// Different apertures use different grid patterns:
//
// APERTURE 3 (ISEA3H):
//   - Even resolutions: "aligned" grid where all integer (i,j) are valid cells
//   - Odd resolutions: "offset" grid where only 1/3 of cells are valid
//                      (those where (i+j) % 3 == 0)
//   - Cell count: N = 10 * 3^res + 2
//   - Grid dim: sqrt(3)^res for aligned, sqrt(3)^(res+1) for offset
//
// APERTURE 4 (ISEA4H):
//   - Always "aligned" (Class I) - all (i,j) pairs valid
//   - Cell count: N = 10 * 4^res + 2
//   - Grid dim: 2^res
//
// APERTURE 7 (ISEA7H):
//   - Even resolutions: Class III-I
//   - Odd resolutions: Class III-II
//   - Cell count: N = 10 * 7^res + 2
//   - Grid dim: sqrt(7)^res for Class III-I, sqrt(21)^res for Class III-II
// ============================================================================

// Check if aperture 3 resolution uses aligned (even) or offset (odd) grid
static inline bool is_aligned_grid_ap3(int resolution) {
    return (resolution % 2) == 0;
}

// Calculate max grid index for aperture 3
// Class I (even resolution): maxI = sqrt(3)^res - 1
// Class II (odd resolution): maxI = 3 * sqrt(3)^(res-1) - 1
// This corresponds to numI = maxI + 1 cells per dimension,
// with total cells per quad = numI * numI / 3 (for Class II)
static long long calc_max_grid_dim_ap3(int resolution) {
    if (resolution == 0) return 0;

    bool is_class1 = is_aligned_grid_ap3(resolution);

    // Compute sqrt(3)^resolution
    double scale = 1.0;
    for (int r = 1; r <= resolution; r++) {
        scale *= 1.7320508075688772935;  // sqrt(3)
    }

    if (is_class1) {
        // Class I: maxI = sqrt(3)^res - 1
        return static_cast<long long>(scale + 0.000001) - 1;
    } else {
        // Class II: maxI = 3 * sqrt(3)^(res-1) - 1 = sqrt(3)^(res-1) * 3 - 1
        // Since sqrt(3)^res = sqrt(3)^(res-1) * sqrt(3),
        // we have sqrt(3)^(res-1) = scale / sqrt(3)
        // maxI = (scale / sqrt(3)) * 3 - 1 = scale * sqrt(3) - 1
        double maxI = scale * 1.7320508075688772935 - 1.0;
        return static_cast<long long>(maxI + 0.000001);
    }
}

// Calculate grid dimension for aperture 4
static long long calc_max_grid_dim_ap4(int resolution) {
    if (resolution == 0) return 0;
    return (1LL << resolution) - 1;  // 2^res - 1
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

// The 30-degree Class II lattice of aperture 3's odd resolutions: one substrate
// point in three is a cell, those with (i + j) % 3 == 0.
static const SubstrateLattice kOffsetLatticeAp3 = {3, 2};

// Substrate lattice of a pure single-aperture grid. Aperture 7 stores
// surrogates rather than substrate coordinates and packs them aligned.
static SubstrateLattice lattice_for_aperture(int aperture, int resolution) {
    if (aperture == 3 && !is_aligned_grid_ap3(resolution)) return kOffsetLatticeAp3;
    return kAlignedLattice;
}

// Calculate cell count and offset per quad for any aperture: each of the 10
// quads owns aperture^res cells and the two poles bring the total to
// 10 * aperture^res + 2.
static void calc_grid_params(int resolution, int aperture,
                             uint64_t& nCells, uint64_t& offsetPerQuad) {
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }
    nCells = 10;
    for (int r = 0; r < resolution; r++) {
        nCells *= aperture;
    }
    nCells += 2;
    offsetPerQuad = (nCells - 2) / 10;
}

// Grid dimension (per-axis) for a given resolution/aperture. Used by every
// cell-ID <-> (quad,i,j) conversion below to select which calc_max_grid_dim_*
// helper applies. Aperture 7 reports its Class I substrate scale; its cell
// index comes from the surrogate rather than a row-major walk of that box.
static long long grid_dim_for_aperture(int resolution, int aperture) {
    if (aperture == 3) {
        return calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        return calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        return hexify::ap7_classI_scale(resolution);
    }
}

// Decode a validated 1-based cell ID into (quad, i, j). Throws via
// Rcpp::stop() if cell_id_raw is non-finite (NA/NaN/Inf) or outside
// [1, nCells], so callers never index a static lookup table with a garbage
// quad/i/j derived from an out-of-range or NA cell ID.
//
// handle_ap7_south_pole preserves existing per-caller behavior: some call
// sites special-case quad 11 under aperture 7 as the south pole pentagon
// (i=j=0) and some decode it like any other cell; this parameter keeps that
// distinction rather than silently changing either behavior during the
// dedup of this decode logic.
static void decode_cell_id(double cell_id_raw, int resolution, int aperture,
                            long long dim, uint64_t offsetPerQuad, uint64_t nCells,
                            const SubstrateLattice& lat, bool handle_ap7_south_pole,
                            int& quad, long long& i, long long& j) {
    if (!std::isfinite(cell_id_raw) || cell_id_raw < 1.0 ||
        cell_id_raw > static_cast<double>(nCells)) {
        Rcpp::stop("cell_id must be a finite value in [1, %.0f] for resolution %d, aperture %d",
                   static_cast<double>(nCells), resolution, aperture);
    }

    uint64_t idx = static_cast<uint64_t>(cell_id_raw);
    idx--;  // Convert to 0-based

    if (idx == 0) {
        // First cell: quad 0 (north pole), i=0, j=0
        quad = 0;
        i = 0;
        j = 0;
        return;
    }

    // Adjust for quad 0
    idx--;

    // Determine quad
    quad = static_cast<int>(idx / offsetPerQuad) + 1;
    idx -= (quad - 1) * offsetPerQuad;

    if (handle_ap7_south_pole && quad == 11 && aperture == 7) {
        // South pole pentagon
        i = 0;
        j = 0;
    } else if (aperture == 7) {
        // Decode to surrogate (i,j stored as surrogates for ap7)
        hexify::ap7_quad_index_to_surrogate(idx, resolution, i, j);
    } else {
        ij_from_cell_index(idx, dim, lat, i, j);
    }
}

// [[Rcpp::export]]
NumericVector cpp_quad_ij_to_cell(IntegerVector quad, NumericVector i,
                                   NumericVector j, int resolution, int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_quad_ij_to_cell: aperture must be 3, 4, or 7");
    }

    int n = quad.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim = grid_dim_for_aperture(resolution, aperture);

    // Check if using offset grid (only aperture 3 odd resolutions)
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        // A coordinate that has stepped outside its quad names a cell of a
        // neighbouring quad, so re-express it there before packing. The walk
        // up and down the cell hierarchy reaches these: an ancestor of a cell
        // near a quad edge need not lie in the same quad. Coordinates already
        // inside their quad pass through unchanged.
        if (!hexify::quad_ij_canonicalize(q, ii, jj, aperture, resolution)) {
            // Outside every adjacent quad, which is where the icosahedron
            // folds at a vertex. No cell owns the coordinate.
            result[k] = NA_REAL;
            continue;
        }

        uint64_t offset = 0;
        if (q > 0) {
            offset = 1 + (q - 1) * offsetPerQuad;
        }

        // 2D cell index within quad. For ap7 the input (i,j) are surrogates.
        uint64_t bnd2D_idx = (aperture == 7)
            ? hexify::ap7_surrogate_to_quad_index(ii, jj, resolution)
            : cell_index_2d(ii, jj, dim, sub_lat);

        // Final cell ID (1-based)
        uint64_t cell_id = offset + bnd2D_idx + 1;

        result[k] = static_cast<double>(cell_id);
    }

    return result;
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell(NumericVector lon, NumericVector lat,
                                  int resolution, int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_lonlat_to_cell: aperture must be 3, 4, or 7");
    }

    int n = lon.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim = grid_dim_for_aperture(resolution, aperture);

    // Aperture 3 odd resolutions use offset grid; all others use aligned
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        hexify::ProjectionResult fwd = hexify::snyder_forward(lon[k], lat[k]);

        uint64_t bnd2D_seq;
        int quad;

        if (aperture == 7) {
            // AP7: exact-integer quantization straight to the surrogate (the
            // resolution-r cell IJK). icosa_tri_to_quad_ij() now returns the
            // canonical surrogate directly (clean Class I quantize + DGGRID
            // edgeTable canonicalization + exact coarsen).
            long long sur_i, sur_j;
            hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x,
                                         fwd.icosa_triangle_y,
                                         7, resolution, quad, sur_i, sur_j);
            bnd2D_seq = hexify::ap7_surrogate_to_quad_index(sur_i, sur_j, resolution);
        } else {
            // AP3/AP4: standard substrate-based encoding
            long long i, j;
            hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                         aperture, resolution, quad, i, j);
            bnd2D_seq = cell_index_2d(i, j, dim, sub_lat);
        }

        // Calculate cell ID offset within quad
        uint64_t offset = 0;
        if (quad > 0) {
            offset = 1 + (quad - 1) * offsetPerQuad;
        }

        uint64_t cid = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(cid);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_lonlat(NumericVector cell_id, int resolution,
                              int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_lonlat: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    NumericVector lon(n);
    NumericVector lat(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim = grid_dim_for_aperture(resolution, aperture);

    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/true,
                        quad, i, j);

        // Convert to lon/lat
        double quad_x, quad_y;
        if (aperture == 7) {
            // AP7: surrogate IJ → quad XY directly
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, quad_x, quad_y);
        } else {
            // AP3/AP4: substrate IJ → quad XY
            hexify::quad_ij_to_xy(quad, i, j, aperture, resolution, quad_x, quad_y);
        }

        // quad_xy -> icosa triangle coords
        int icosa_triangle_face;
        double icosa_triangle_x, icosa_triangle_y;
        if (!hexify::try_quad_xy_to_icosa_tri(quad, quad_x, quad_y, icosa_triangle_face, icosa_triangle_x, icosa_triangle_y)) {
            // Surrogate center falls in an invalid region (e.g. pentagon gap)
            lon[k] = NA_REAL;
            lat[k] = NA_REAL;
            continue;
        }

        // Step 3: icosa triangle coords -> lon/lat
        auto ll = hexify::face_xy_to_ll(icosa_triangle_x, icosa_triangle_y, icosa_triangle_face);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }

    return DataFrame::create(
        _["lon_deg"] = lon,
        _["lat_deg"] = lat
    );
}

// ============================================================================
// Cell ID to Quad IJ Conversion
// ============================================================================
// Converts cell IDs to Quad IJ coordinates.
// This is the inverse of cpp_quad_ij_to_cell.
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_cell_to_quad_ij(NumericVector cell_id, int resolution,
                               int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_quad_ij: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_i(n);
    NumericVector out_j(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim = grid_dim_for_aperture(resolution, aperture);

    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/true,
                        quad, i, j);

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
DataFrame cpp_cell_to_quad_xy(NumericVector cell_id, int resolution,
                               int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_quad_xy: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_qx(n);
    NumericVector out_qy(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Convert to Quad XY
        double quad_x, quad_y;
        if (aperture == 7) {
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, quad_x, quad_y);
        } else {
            hexify::quad_ij_to_xy(quad, i, j, aperture, resolution, quad_x, quad_y);
        }

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
NumericVector cpp_quad_xy_to_cell(IntegerVector quad, NumericVector quad_x,
                                   NumericVector quad_y, int resolution,
                                   int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_quad_xy_to_cell: aperture must be 3, 4, or 7");
    }

    int n = quad.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);

    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        double qx = quad_x[k];
        double qy = quad_y[k];

        uint64_t bnd2D_seq;
        int out_quad;

        if (aperture == 7) {
            // AP7: exact-integer quantization straight to the surrogate.
            long long sur_i, sur_j;
            hexify::quad_xy_to_ij(q, qx, qy, 7, resolution, out_quad, sur_i, sur_j);
            bnd2D_seq = hexify::ap7_surrogate_to_quad_index(sur_i, sur_j, resolution);
        } else {
            // AP3/AP4: standard substrate path
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            hexify::quad_xy_to_icosa_tri(q, qx, qy, icosa_triangle_face,
                                         icosa_triangle_x, icosa_triangle_y);
            long long i, j;
            hexify::icosa_tri_to_quad_ij(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y,
                                         aperture, resolution, out_quad, i, j);
            bnd2D_seq = cell_index_2d(i, j, dim, sub_lat);
        }

        // Calculate cell ID
        uint64_t offset = 0;
        if (out_quad > 0) {
            offset = 1 + (out_quad - 1) * offsetPerQuad;
        }

        uint64_t cid = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(cid);
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
DataFrame cpp_cell_to_icosa_tri(NumericVector cell_id, int resolution,
                                 int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_icosa_tri: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    IntegerVector out_face(n);
    NumericVector out_tx(n);
    NumericVector out_ty(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Convert to Quad XY
        double quad_x, quad_y;
        if (aperture == 7) {
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, quad_x, quad_y);
        } else {
            hexify::quad_ij_to_xy(quad, i, j, aperture, resolution, quad_x, quad_y);
        }

        // Convert Quad XY → Icosa Triangle
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
DataFrame cpp_quad_ij_to_icosa_tri(IntegerVector quad, NumericVector i,
                                    NumericVector j, int resolution,
                                    int aperture) {
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
// Polygon Corner Generation from Cell ID
// ============================================================================
// Generates hexagon corner coordinates for a vector of cell ID values.
// Uses the same coordinate transformation:
// 1. Compute vertex offsets in quad 2D space
// 2. Convert each vertex through quad_xy → face coords → lon/lat
// ============================================================================

// Class I (flat-top) hex vertex offsets in unscaled 2D grid coordinates
// r_ = 1/sqrt(3), r2 = r_/2
// Vertices counter-clockwise from top: (0, r_), (-0.5, r2), (-0.5, -r2), (0, -r_), (0.5, -r2), (0.5, r2)
constexpr double kHexR = 0.57735026918962576451;  // 1/sqrt(3)
constexpr double kHexR2 = 0.28867513459481288225; // 1/(2*sqrt(3))

// Class I vertex offsets (flat-top hexagon)
static const double kClass1VertexX[6] = { 0.0, -0.5, -0.5,  0.0,  0.5, 0.5};
static const double kClass1VertexY[6] = {kHexR, kHexR2, -kHexR2, -kHexR, -kHexR2, kHexR2};

// Class II vertex offsets (Class I rotated 30° CCW for pointy-top hexagon)
// These offsets work correctly for quad 3 (face 2 with dazh=0).
// TODO: Other quads may need rotation adjustment based on face azimuth.
// Vertex winding order for polygon construction (counter-clockwise)
static const double kClass2VertexX[6] = {-kHexR2, kHexR2, kHexR, kHexR2, -kHexR2, -kHexR};
static const double kClass2VertexY[6] = {-0.5, -0.5, 0.0, 0.5, 0.5, 0.0};

// Pentagon vertex skip indices for icosahedral vertex cells
// All cells at (i=0, j=0) in any quad are pentagon cells (12 total, one per icosahedral vertex)
// The invalid region depends on the quad:
//   Quads 0-5: Region 3 is invalid -> skip vertex 3
//   Quads 6-11: Region 4 is invalid -> skip vertex 2
// These indices are based on the Class I vertex layout where vertices map to regions:
//   Vertex 0 -> Region 0 (upper), Vertex 1 -> Region 5 (upper-left)
//   Vertex 2 -> Region 4 (lower-left), Vertex 3 -> Region 3 (lower)
//   Vertex 4 -> Region 2 (lower-right), Vertex 5 -> Region 1 (upper-right)
constexpr int kPentagonSkipVertexRegion3 = 3;  // Skip vertex 3 for quads 0-5 (region 3 invalid)
constexpr int kPentagonSkipVertexRegion4 = 2;  // Skip vertex 2 for quads 6-11 (region 4 invalid)

// [[Rcpp::export]]
DataFrame cpp_cell_to_polygon(NumericVector cell_id, int resolution,
                               int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_polygon: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();

    // Each cell has 7 vertices (6 corners + 1 to close the polygon)
    int total_vertices = n * 7;

    NumericVector out_cell_id(total_vertices);
    NumericVector out_lon(total_vertices);
    NumericVector out_lat(total_vertices);
    IntegerVector out_order(total_vertices);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim = grid_dim_for_aperture(resolution, aperture);

    // Check if using offset grid (only aperture 3 odd resolutions)
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    // Determine grid class and select vertex offsets
    // Class I = flat-top hexagons, Class II = pointy-top (30° rotated)
    bool is_class1 = (aperture == 4) || (aperture == 3 && resolution % 2 == 0);
    const double* vertexX = is_class1 ? kClass1VertexX : kClass2VertexX;
    const double* vertexY = is_class1 ? kClass1VertexY : kClass2VertexY;

    // Calculate scale factors:
    // For Class II, the grid (i,j) indices are in substrate coordinates (scaled by sqrt(3))
    // relative to the backFrame. Vertex offsets are in backFrame coordinates.
    // So we need two scales: one for the center, one for the vertex offsets.
    double center_scale;   // Scale to convert (i,j) grid coords to Quad XY
    double vertex_scale;   // Scale to convert vertex offsets to Quad XY
    if (aperture == 3) {
        if (is_class1) {
            center_scale = std::pow(hexify::kSqrt3, resolution);
            vertex_scale = center_scale;
        } else {
            // Class II: substrate scale for center, backFrame scale for vertices
            center_scale = std::pow(hexify::kSqrt3, resolution + 1);  // sqrt(3)^(res+1)
            vertex_scale = std::pow(hexify::kSqrt3, resolution);      // sqrt(3)^res
        }
    } else if (aperture == 4) {
        center_scale = std::pow(2.0, resolution);
        vertex_scale = center_scale;
    } else {  // aperture == 7
        // Substrate scale for center coords (decoded substrate i,j)
        double base_scale = std::pow(std::sqrt(7.0), resolution);
        bool is_class3i = (resolution % 2 == 0);
        double substrate_mult = is_class3i ? hexify::kSqrt7 : hexify::kSqrt21;
        center_scale = base_scale * substrate_mult;
        vertex_scale = base_scale;
    }

    int out_idx = 0;

    for (int k = 0; k < n; k++) {
        double orig_cell_id = cell_id[k];

        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Get cell center in Quad XY coordinates
        double qx_center, qy_center;
        if (aperture == 7) {
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, qx_center, qy_center);
        } else {
            double grid_x = static_cast<double>(i) - 0.5 * static_cast<double>(j);
            double grid_y = static_cast<double>(j) * hexify::kSin60;
            qx_center = grid_x / center_scale;
            qy_center = grid_y / center_scale;
        }

        // Generate 6 corners: add vertex offsets (which are in backFrame scale)
        double first_lon = 0.0, first_lat = 0.0;
        for (int c = 0; c < 6; c++) {
            // Vertex offsets are in backFrame coordinates, normalize to quad_xy
            double qx_vertex = qx_center + vertexX[c] / vertex_scale;
            double qy_vertex = qy_center + vertexY[c] / vertex_scale;

            // Convert quad_xy to icosa triangle coords
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            hexify::quad_xy_to_icosa_tri(quad, qx_vertex, qy_vertex, icosa_triangle_face, icosa_triangle_x, icosa_triangle_y);

            // Convert to lon/lat
            auto ll = hexify::face_xy_to_ll(icosa_triangle_x, icosa_triangle_y, icosa_triangle_face);

            out_cell_id[out_idx] = orig_cell_id;
            out_lon[out_idx] = ll.first;
            out_lat[out_idx] = ll.second;
            out_order[out_idx] = c + 1;

            if (c == 0) {
                first_lon = ll.first;
                first_lat = ll.second;
            }
            out_idx++;
        }

        // Close the polygon by repeating the first vertex
        out_cell_id[out_idx] = orig_cell_id;
        out_lon[out_idx] = first_lon;
        out_lat[out_idx] = first_lat;
        out_order[out_idx] = 7;
        out_idx++;
    }

    return DataFrame::create(
        _["hex_id"] = out_cell_id,
        _["lon"] = out_lon,
        _["lat"] = out_lat,
        _["order"] = out_order
    );
}

// Vectorized version that returns a list of coordinate matrices for efficiency
// [[Rcpp::export]]
List cpp_cell_to_corners(NumericVector cell_id, int resolution,
                          int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_corners: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    // Determine grid class and select vertex offsets
    bool is_class1 = (aperture == 4) || (aperture == 3 && resolution % 2 == 0);
    const double* vertexX = is_class1 ? kClass1VertexX : kClass2VertexX;
    const double* vertexY = is_class1 ? kClass1VertexY : kClass2VertexY;

    // Calculate scale factors (same logic as cpp_cell_to_polygon)
    double center_scale;
    double vertex_scale;
    if (aperture == 3) {
        if (is_class1) {
            center_scale = std::pow(hexify::kSqrt3, resolution);
            vertex_scale = center_scale;
        } else {
            center_scale = std::pow(hexify::kSqrt3, resolution + 1);
            vertex_scale = std::pow(hexify::kSqrt3, resolution);
        }
    } else if (aperture == 4) {
        center_scale = std::pow(2.0, resolution);
        vertex_scale = center_scale;
    } else {  // aperture == 7
        double base_scale = std::pow(std::sqrt(7.0), resolution);
        bool is_class3i = (resolution % 2 == 0);
        double substrate_mult = is_class3i ? hexify::kSqrt7 : hexify::kSqrt21;
        center_scale = base_scale * substrate_mult;
        vertex_scale = base_scale;
    }

    // Return a list of n elements, each element is a Nx2 matrix (lon, lat)
    // Hexagons: 7 rows (6 vertices + closing)
    // Pentagons: 6 rows (5 vertices + closing) - at icosahedral vertices
    List result(n);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Check if this is a pentagon cell (at icosahedral vertex)
        bool is_pentagon = (i == 0 && j == 0);
        int skip_vertex = -1;
        if (is_pentagon) {
            skip_vertex = (quad <= 5) ? kPentagonSkipVertexRegion3 : kPentagonSkipVertexRegion4;
        }

        // Get cell center in Quad XY coordinates
        double qx_center, qy_center;
        if (aperture == 7) {
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, qx_center, qy_center);
        } else {
            double grid_x = static_cast<double>(i) - 0.5 * static_cast<double>(j);
            double grid_y = static_cast<double>(j) * hexify::kSin60;
            qx_center = grid_x / center_scale;
            qy_center = grid_y / center_scale;
        }

        // Pentagon: 6 rows (5 vertices + closing)
        // Hexagon: 7 rows (6 vertices + closing)
        int n_vertices = is_pentagon ? 5 : 6;
        int n_rows = n_vertices + 1;  // +1 for closing vertex
        NumericMatrix coords(n_rows, 2);

        int out_idx = 0;  // Output index (may differ from c if skipping a vertex)
        for (int c = 0; c < 6; c++) {
            // Skip the invalid vertex for pentagon cells
            if (c == skip_vertex) {
                continue;
            }

            // Vertex offsets are in backFrame coordinates, normalize to quad_xy
            double qx_vertex = qx_center + vertexX[c] / vertex_scale;
            double qy_vertex = qy_center + vertexY[c] / vertex_scale;

            // Try to convert quad_xy to icosa triangle coords
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            if (hexify::try_quad_xy_to_icosa_tri(quad, qx_vertex, qy_vertex, icosa_triangle_face, icosa_triangle_x, icosa_triangle_y)) {
                // Normal case: convert to lon/lat via Snyder inverse
                auto ll = hexify::face_xy_to_ll(icosa_triangle_x, icosa_triangle_y, icosa_triangle_face);
                coords(out_idx, 0) = ll.first;
                coords(out_idx, 1) = ll.second;
            } else {
                // Edge case: vertex is in invalid region (quad boundary)
                // Fall back to using the center's face and projecting vertex offset
                // This is an approximation but works for edge cells
                int center_icosa_triangle_face;
                double center_icosa_triangle_x, center_icosa_triangle_y;
                if (hexify::try_quad_xy_to_icosa_tri(quad, qx_center, qy_center,
                                                     center_icosa_triangle_face, center_icosa_triangle_x, center_icosa_triangle_y)) {
                    // Project vertex through center's face
                    double vertex_icosa_triangle_x = center_icosa_triangle_x + (qx_vertex - qx_center);
                    double vertex_icosa_triangle_y = center_icosa_triangle_y + (qy_vertex - qy_center);
                    auto ll = hexify::face_xy_to_ll(vertex_icosa_triangle_x, vertex_icosa_triangle_y, center_icosa_triangle_face);
                    coords(out_idx, 0) = ll.first;
                    coords(out_idx, 1) = ll.second;
                } else {
                    // Both center and vertex are invalid - use NaN
                    coords(out_idx, 0) = NA_REAL;
                    coords(out_idx, 1) = NA_REAL;
                }
            }
            out_idx++;
        }

        // Close polygon
        coords(n_vertices, 0) = coords(0, 0);
        coords(n_vertices, 1) = coords(0, 1);

        // Extend polar cells to reach 90°N or 90°S
        // Find ALL vertices above 88° and extend them to the pole
        // Also normalize longitude to 0° at poles (all longitudes = same point)
        for (int v = 0; v <= n_vertices; v++) {
            double lat = coords(v, 1);
            if (!NumericVector::is_na(lat)) {
                if (lat > 88.0) {
                    coords(v, 1) = 90.0;   // Extend to North Pole
                    coords(v, 0) = 0.0;    // Normalize lon (all lons same at pole)
                } else if (lat < -88.0) {
                    coords(v, 1) = -90.0;  // Extend to South Pole
                    coords(v, 0) = 0.0;    // Normalize lon
                }
            }
        }

        colnames(coords) = CharacterVector::create("lon", "lat");
        result[k] = coords;
    }

    return result;
}

// ============================================================================
// Mixed Aperture Sequence Cell ID Conversion
// ============================================================================
// A mixed grid refines by a different aperture at each level, given as the
// sequence ap_seq (see aperture_sequence.h). ISEA43H is aperture 4 for the
// first 'mixed_aperture_level' resolutions and aperture 3 for the rest.
//
// Cell count: N = 10 * (product of the apertures) + 2
// ============================================================================

// Which substrate points a mixed sequence's cells sit on
static SubstrateLattice lattice_for_mixed(const std::vector<int>& ap_seq) {
    return sublattice_of(hexify::hex_form_sequence(ap_seq));
}

// An aperture sequence from R. Entry 0 names the base grid and entries 1.. are
// the refinement steps, so the resolution is one less than the length; the
// entries themselves are validated by hex_form_sequence().
static std::vector<int> as_ap_seq(const IntegerVector& ap_seq, const char* fn) {
    if (ap_seq.size() < 1) {
        Rcpp::stop(std::string(fn) + ": ap_seq must name at least the base grid");
    }
    return std::vector<int>(ap_seq.begin(), ap_seq.end());
}

// Cell count and per-quad offset for a mixed aperture sequence
static void calc_grid_params_mixed(const std::vector<int>& ap_seq,
                                    uint64_t& nCells, uint64_t& offsetPerQuad) {
    int resolution = static_cast<int>(ap_seq.size()) - 1;
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }

    nCells = 10;
    for (size_t k = 1; k < ap_seq.size(); k++) {
        nCells *= static_cast<uint64_t>(ap_seq[k]);
    }

    nCells += 2;
    offsetPerQuad = (nCells - 2) / 10;
}

// [[Rcpp::export]]
double cpp_ap_seq_edge_dim(IntegerVector ap_seq_in) {
    std::vector<int> ap_seq = as_ap_seq(ap_seq_in, "cpp_ap_seq_edge_dim");
    return static_cast<double>(hexify::quad_edge_coord_mixed(ap_seq));
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_seq(NumericVector lon, NumericVector lat,
                                     IntegerVector ap_seq_in) {
    std::vector<int> ap_seq = as_ap_seq(ap_seq_in, "cpp_lonlat_to_cell_seq");

    int n = lon.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_mixed(ap_seq, nCells, offsetPerQuad);

    // Grid dimension
    long long dim = hexify::quad_edge_coord_mixed(ap_seq);

    // Check if using offset grid
    SubstrateLattice sub_lat = lattice_for_mixed(ap_seq);

    for (int k = 0; k < n; k++) {
        hexify::ProjectionResult fwd = hexify::snyder_forward(lon[k], lat[k]);

        int quad_pre;
        double quad_x, quad_y;
        hexify::icosa_tri_to_quad_xy(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                     quad_pre, quad_x, quad_y);

        int quad;
        long long i, j;
        hexify::quad_xy_to_ij_mixed(quad_pre, quad_x, quad_y, ap_seq, quad, i, j);

        // Calculate cell ID offset within quad
        uint64_t offset = 0;
        if (quad > 0) {
            offset = 1 + (quad - 1) * offsetPerQuad;
        }

        // Calculate 2D cell index based on grid pattern
        uint64_t bnd2D_seq = cell_index_2d(i, j, dim, sub_lat);

        uint64_t cid = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(cid);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_lonlat_seq(NumericVector cell_id, IntegerVector ap_seq_in) {
    std::vector<int> ap_seq = as_ap_seq(ap_seq_in, "cpp_cell_to_lonlat_seq");
    int resolution = static_cast<int>(ap_seq.size()) - 1;

    int n = cell_id.size();
    NumericVector lon(n);
    NumericVector lat(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_mixed(ap_seq, nCells, offsetPerQuad);

    // Grid dimension
    long long dim = hexify::quad_edge_coord_mixed(ap_seq);

    // Check if using offset grid
    SubstrateLattice sub_lat = lattice_for_mixed(ap_seq);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        // Any aperture other than 7 selects decode_cell_id's substrate branch,
        // which is the one a mixed sequence stores its cells on.
        decode_cell_id(cell_id[k], resolution, /*aperture=*/3, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Convert quad IJ to lon/lat via quad_xy -> icosa triangle -> lon/lat
        double quad_x, quad_y;
        hexify::quad_ij_to_xy_mixed(quad, i, j, ap_seq, quad_x, quad_y);

        int icosa_triangle_face;
        double icosa_triangle_x, icosa_triangle_y;
        hexify::quad_xy_to_icosa_tri(quad, quad_x, quad_y, icosa_triangle_face, icosa_triangle_x, icosa_triangle_y);

        auto ll = hexify::face_xy_to_ll(icosa_triangle_x, icosa_triangle_y, icosa_triangle_face);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }

    return DataFrame::create(
        _["lon_deg"] = lon,
        _["lat_deg"] = lat
    );
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_quad_ij_seq(NumericVector cell_id, IntegerVector ap_seq_in) {
    std::vector<int> ap_seq = as_ap_seq(ap_seq_in, "cpp_cell_to_quad_ij_seq");
    int resolution = static_cast<int>(ap_seq.size()) - 1;

    int n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_i(n);
    NumericVector out_j(n);
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_mixed(ap_seq, nCells, offsetPerQuad);
    long long dim = hexify::quad_edge_coord_mixed(ap_seq);
    SubstrateLattice sub_lat = lattice_for_mixed(ap_seq);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        // Any aperture other than 7 selects decode_cell_id's substrate branch,
        // which is the one a mixed sequence stores its cells on.
        decode_cell_id(cell_id[k], resolution, /*aperture=*/3, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);
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

// [[Rcpp::export]]
NumericVector cpp_quad_ij_to_cell_seq(IntegerVector quad, NumericVector i,
                                      NumericVector j, IntegerVector ap_seq_in) {
    std::vector<int> ap_seq = as_ap_seq(ap_seq_in, "cpp_quad_ij_to_cell_seq");

    int n = quad.size();
    NumericVector result(n);
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_mixed(ap_seq, nCells, offsetPerQuad);
    long long dim = hexify::quad_edge_coord_mixed(ap_seq);
    SubstrateLattice sub_lat = lattice_for_mixed(ap_seq);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        // Mirror the cell-ID packing of cpp_lonlat_to_cell_seq(): offset by quad,
        // then the 2D boundary sequence index within the quad's substrate.
        uint64_t offset = 0;
        if (q > 0) {
            offset = 1 + (q - 1) * offsetPerQuad;
        }

        uint64_t bnd2D_seq = cell_index_2d(ii, jj, dim, sub_lat);

        uint64_t cid = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(cid);
    }

    return result;
}

// ============================================================================
// Neighbor Finding (v0.7.0)
// ============================================================================

// Helper: convert (quad, i, j) to cell_id
// For ap7: (i,j) are surrogates. For ap3/4: (i,j) are substrates.
static double encode_cell_id(int quad, long long i, long long j,
                              uint64_t offsetPerQuad, long long dim,
                              const SubstrateLattice& lat, int aperture, int resolution) {
    if (quad == 0 && i == 0 && j == 0) {
        return 1.0;
    }

    uint64_t offset = 1 + (quad - 1) * offsetPerQuad;
    uint64_t bnd2D_seq = (aperture == 7)
        ? hexify::ap7_surrogate_to_quad_index(i, j, resolution)
        : cell_index_2d(i, j, dim, lat);

    return static_cast<double>(offset + bnd2D_seq + 1);
}

// Resolution 0 is the 12 base cells, one per icosahedron vertex. Every one of
// them is a pentagon and each quad holds a single cell, so adjacency there is
// the icosahedron's vertex graph rather than a step through a quad frame.
// Row q lists the quads sharing an edge with quad q.
static const int kBaseCellNeighbors[12][5] = {
    { 1,  2,  3,  4,  5},
    { 0,  2,  5,  6, 10},
    { 0,  1,  3,  6,  7},
    { 0,  2,  4,  7,  8},
    { 0,  3,  5,  8,  9},
    { 0,  1,  4,  9, 10},
    { 1,  2,  7, 10, 11},
    { 2,  3,  6,  8, 11},
    { 3,  4,  7,  9, 11},
    { 4,  5,  8, 10, 11},
    { 1,  5,  6,  9, 11},
    { 6,  7,  8,  9, 10}
};

// [[Rcpp::export]]
Rcpp::List cpp_get_neighbors_isea(Rcpp::NumericVector cell_id, int resolution,
                                   int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        Rcpp::stop("cpp_get_neighbors_isea: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    Rcpp::List out(n);

    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);

    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    long long max_ij = hexify::get_max_ij(aperture, resolution);

    // The aligned lattice writes (i, j) as x = i - j/2, y = j * sin60, so the
    // six cells one unit away sit at +/-(1,0), +/-(0,1) and +/-(1,1). Aperture
    // 4, aperture 3's even resolutions and aperture 7's surrogates all live on
    // this lattice.
    static const long long class1_offsets[6][2] = {
        { 1,  0}, { 1,  1}, { 0,  1},
        {-1,  0}, {-1, -1}, { 0, -1}
    };

    // Aperture 3's odd resolutions carry cells on the 30-degree Class II
    // lattice, where one substrate point in three is a cell and adjacent cells
    // stand sqrt(3) substrate units apart.
    static const long long class2_offsets[6][2] = {
        { 2,  1}, { 1,  2}, {-1,  1},
        {-2, -1}, {-1, -2}, { 1, -1}
    };

    for (int k = 0; k < n; k++) {
        double cell_id_raw = cell_id[k];
        if (!std::isfinite(cell_id_raw) || cell_id_raw < 1.0 ||
            cell_id_raw > static_cast<double>(nCells)) {
            Rcpp::stop("cell_id must be a finite value in [1, %.0f] for resolution %d, aperture %d",
                       static_cast<double>(nCells), resolution, aperture);
        }
        uint64_t idx = static_cast<uint64_t>(cell_id_raw);
        idx--;

        std::vector<double> neighbor_ids;
        neighbor_ids.reserve(6);

        if (resolution == 0) {
            Rcpp::NumericVector base(5);
            for (int d = 0; d < 5; d++) {
                base[d] = kBaseCellNeighbors[idx][d] + 1;
            }
            out[k] = base;
            continue;
        }

        int quad;
        long long i, j;

        if (idx == 0 || idx == nCells - 1) {
            // Quads 0 and 11 hold a single cell each -- the icosahedron vertex
            // where five quads meet -- so their own frame carries no offsets to
            // step through. Read the vertex from an adjacent quad, where it
            // sits at the corner coordinate the edge tables fold into the pole,
            // and take the six offsets from there.
            long long edge = (aperture == 7)
                ? hexify::ap7_classI_scale(resolution)
                : hexify::get_max_ij(aperture, resolution) + 1;
            bool north = (idx == 0);
            long long ci = north ? 0 : edge;
            long long cj = north ? edge : 0;
            quad = north ? 1 : 6;
            if (aperture == 7) {
                hexify::ap7_substrate_to_surrogate_ijk(ci, cj, resolution, i, j);
            } else {
                i = ci;
                j = cj;
            }
        } else {
            idx--;
            quad = static_cast<int>(idx / offsetPerQuad) + 1;
            uint64_t within_quad = idx - (quad - 1) * offsetPerQuad;

            if (aperture == 7) {
                hexify::ap7_quad_index_to_surrogate(within_quad, resolution, i, j);
            } else {
                ij_from_cell_index(within_quad, dim, sub_lat, i, j);
            }
        }

        const long long (*offsets)[2] =
            (aperture == 3 && !is_aligned_grid_ap3(resolution))
                ? class2_offsets
                : class1_offsets;

        for (int d = 0; d < 6; d++) {
            long long ni = i + offsets[d][0];
            long long nj = j + offsets[d][1];

            bool in_bounds = (aperture == 7)
                ? hexify::ap7_surrogate_in_quad(ni, nj, resolution)
                : (ni >= 0 && nj >= 0 && ni <= max_ij && nj <= max_ij);

            if (in_bounds) {
                neighbor_ids.push_back(
                    encode_cell_id(quad, ni, nj,
                                   offsetPerQuad, dim, sub_lat,
                                   aperture, resolution));
                continue;
            }

            // The neighbour sits in another quad: send its centre back through
            // the forward pipeline, which names the quad that owns it.
            double nbr_qx, nbr_qy;
            hexify::quad_ij_to_xy(quad, ni, nj, aperture, resolution,
                                   nbr_qx, nbr_qy);

            int tri_face;
            double tri_x, tri_y;
            if (!hexify::try_quad_xy_to_icosa_tri(quad, nbr_qx, nbr_qy,
                                                   tri_face, tri_x, tri_y)) {
                // Under-runs of a quad's own frame have no image in it, so the
                // projection has nowhere to send them. Step to the owning quad
                // through the edge table instead.
                int alt_quad = quad;
                long long alt_i = ni, alt_j = nj;
                if (hexify::quad_ij_canonicalize(alt_quad, alt_i, alt_j,
                                                  aperture, resolution)) {
                    neighbor_ids.push_back(
                        encode_cell_id(alt_quad, alt_i, alt_j,
                                       offsetPerQuad, dim, sub_lat,
                                       aperture, resolution));
                }
                continue;
            }

            auto ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
            auto fwd = hexify::snyder_forward(ll.first, ll.second);
            int final_quad;
            long long final_i, final_j;
            hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x,
                                          fwd.icosa_triangle_y,
                                          aperture, resolution,
                                          final_quad, final_i, final_j);

            neighbor_ids.push_back(
                encode_cell_id(final_quad, final_i, final_j,
                               offsetPerQuad, dim, sub_lat,
                               aperture, resolution));
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
Rcpp::List cpp_get_neighbors_z7(Rcpp::CharacterVector index_ids, int resolution) {
    int n = index_ids.size();
    Rcpp::List out(n);

    for (int k = 0; k < n; k++) {
        if (index_ids[k] == NA_STRING) {
            out[k] = Rcpp::CharacterVector(0);
            continue;
        }

        std::string idx = Rcpp::as<std::string>(index_ids[k]);

        int quadNum;
        long long ci, cj;
        hexify::z7::decode(idx, resolution, quadNum, ci, cj);

        hexify::z7::IVec3D coord(ci, cj);

        std::vector<std::string> neighbors;
        neighbors.reserve(6);

        for (int d = 1; d <= 6; d++) {
            hexify::z7::IVec3D nbr = coord;
            nbr.neighbor(static_cast<hexify::z7::IVec3D::Direction>(d));

            try {
                std::string nbr_idx = hexify::z7::encode(quadNum, nbr.i(), nbr.j(), resolution);
                nbr_idx = hexify::z7::canonical_form(nbr_idx);
                neighbors.push_back(nbr_idx);
            } catch (...) {
                // Skip invalid neighbors (edge/boundary)
            }
        }

        out[k] = Rcpp::wrap(neighbors);
    }

    return out;
}

// ============================================================================
// PLANE Coordinate Conversions
// ============================================================================
// PLANE coordinates represent the unfolded icosahedron in 2D.
// The transformation from Icosa Triangle to PLANE involves:
// 1. Rotate the point by rot60 * 60 degrees
// 2. Translate by the triangle's offset position
//
// This creates a flat map layout ~5.5 x 1.73 units containing all 20 triangles.
// Produces standard ISEA PLANE coordinates for visualization.
// ============================================================================

// [[Rcpp::export]]
DataFrame cpp_icosa_tri_to_plane(IntegerVector icosa_triangle_face,
                                  NumericVector icosa_triangle_x,
                                  NumericVector icosa_triangle_y) {
    int n = icosa_triangle_face.size();
    NumericVector out_px(n);
    NumericVector out_py(n);

    for (int k = 0; k < n; k++) {
        int face = icosa_triangle_face[k];
        if (face < 0 || face >= 20) {
            out_px[k] = NA_REAL;
            out_py[k] = NA_REAL;
            continue;
        }

        double x = icosa_triangle_x[k];
        double y = icosa_triangle_y[k];

        // Get layout parameters for this face
        const hexify::PlaneTriLayout& layout = hexify::kPlaneLayout[face];

        // Rotate by rot60 * 60 degrees
        if (layout.rot60 != 0) {
            double angle_rad = layout.rot60 * 60.0 * hexify::kDegToRad;
            double cos_ang = std::cos(angle_rad);
            double sin_ang = std::sin(angle_rad);
            double x_rot = x * cos_ang - y * sin_ang;
            double y_rot = x * sin_ang + y * cos_ang;
            x = x_rot;
            y = y_rot;
        }

        // Add offset
        out_px[k] = x + layout.offset_x;
        out_py[k] = y + layout.offset_y;
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_plane(NumericVector cell_id, int resolution, int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_cell_to_plane: aperture must be 3, 4, or 7");
    }

    int n = cell_id.size();
    NumericVector out_px(n);
    NumericVector out_py(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim = grid_dim_for_aperture(resolution, aperture);
    SubstrateLattice sub_lat = lattice_for_aperture(aperture, resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, sub_lat, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Convert to Quad XY
        double quad_x, quad_y;
        if (aperture == 7) {
            hexify::surrogate_ij_to_quad_xy_ap7(i, j, resolution, quad_x, quad_y);
        } else {
            hexify::quad_ij_to_xy(quad, i, j, aperture, resolution, quad_x, quad_y);
        }

        // Convert Quad XY → Icosa Triangle
        int tri_face;
        double tri_x, tri_y;
        hexify::quad_xy_to_icosa_tri(quad, quad_x, quad_y, tri_face,
                                     tri_x, tri_y);

        // Convert Icosa Triangle → PLANE
        const hexify::PlaneTriLayout& layout = hexify::kPlaneLayout[tri_face];

        double x = tri_x;
        double y = tri_y;

        // Rotate by rot60 * 60 degrees
        if (layout.rot60 != 0) {
            double angle_rad = layout.rot60 * 60.0 * hexify::kDegToRad;
            double cos_ang = std::cos(angle_rad);
            double sin_ang = std::sin(angle_rad);
            double x_rot = x * cos_ang - y * sin_ang;
            double y_rot = x * sin_ang + y * cos_ang;
            x = x_rot;
            y = y_rot;
        }

        // Add offset
        out_px[k] = x + layout.offset_x;
        out_py[k] = y + layout.offset_y;
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}

// [[Rcpp::export]]
DataFrame cpp_lonlat_to_plane(NumericVector lon, NumericVector lat) {
    int n = lon.size();
    if (lat.size() != n) {
        stop("cpp_lonlat_to_plane: lon and lat must have same length");
    }

    NumericVector out_px(n);
    NumericVector out_py(n);

    for (int k = 0; k < n; k++) {
        // Project to icosahedron
        auto fwd = hexify::snyder_forward(lon[k], lat[k]);

        int face = fwd.face;
        double x = fwd.icosa_triangle_x;
        double y = fwd.icosa_triangle_y;

        // Convert Icosa Triangle → PLANE
        const hexify::PlaneTriLayout& layout = hexify::kPlaneLayout[face];

        // Rotate by rot60 * 60 degrees
        if (layout.rot60 != 0) {
            double angle_rad = layout.rot60 * 60.0 * hexify::kDegToRad;
            double cos_ang = std::cos(angle_rad);
            double sin_ang = std::sin(angle_rad);
            double x_rot = x * cos_ang - y * sin_ang;
            double y_rot = x * sin_ang + y * cos_ang;
            x = x_rot;
            y = y_rot;
        }

        // Add offset
        out_px[k] = x + layout.offset_x;
        out_py[k] = y + layout.offset_y;
    }

    return DataFrame::create(
        _["plane_x"] = out_px,
        _["plane_y"] = out_py
    );
}
