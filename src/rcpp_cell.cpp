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
#include "constants.h"
#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "aperture.h"
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

// 2D cell index for aligned grid (aperture 3 even res, aperture 4 all res)
// Simple row-major ordering where all (i,j) pairs are valid
static uint64_t cell_index_2d_aligned(long long i, long long j, long long dim) {
    return static_cast<uint64_t>(i) * dim + j;
}

// ============================================================================
// Aperture 7: Surrogate-based encoding
// ============================================================================
// Surrogates are the canonical cell coordinates for aperture 7.
// Each ap7 cell corresponds to exactly one surrogate (i,j), unlike substrates
// where the substrate grid has ~7x more positions than actual cells.
//
// The surrogate grid dimension is approximately sqrt(7)^res per axis.
// Due to the 19.1° rotation between surrogate and quad frames, surrogate
// coordinates can be negative. We add an offset to shift to non-negative range.
// ============================================================================

// Offset to shift minimum surrogate coordinate to >= 0
// From rotation analysis: min surrogate ≈ -0.4 * sqrt(7)^res
static long long calc_sur_offset_ap7(int resolution) {
    if (resolution == 0) return 0;
    double S = std::pow(2.6457513110645905905, resolution);  // sqrt(7)^res
    return static_cast<long long>(0.45 * S) + 2;
}

// One-axis dimension of the shifted surrogate grid
// After shifting by offset, surrogates range from [0, dim-1]
// Max surrogate ≈ 1.15 * sqrt(7)^res, so dim = offset + ceil(1.15*S) + 1
static long long calc_sur_dim_ap7(int resolution) {
    if (resolution == 0) return 1;
    double S = std::pow(2.6457513110645905905, resolution);  // sqrt(7)^res
    long long offset = calc_sur_offset_ap7(resolution);
    long long max_sur = static_cast<long long>(1.2 * S) + 2;
    return offset + max_sur + 1;
}

// Encode surrogate (i,j) directly as cell index within a quad
static uint64_t encode_surrogate_ap7(
    long long sur_i, long long sur_j, int resolution, long long sur_dim) {
    long long sur_offset = calc_sur_offset_ap7(resolution);
    return cell_index_2d_aligned(sur_i + sur_offset, sur_j + sur_offset, sur_dim);
}

// Decode cell index to surrogate (i,j) directly
static void decode_surrogate_ap7(
    uint64_t idx, int resolution, long long sur_dim,
    long long& sur_i, long long& sur_j) {
    long long sur_offset = calc_sur_offset_ap7(resolution);
    sur_i = static_cast<long long>(idx / sur_dim) - sur_offset;
    sur_j = static_cast<long long>(idx % sur_dim) - sur_offset;
}

// 2D cell index for offset grid (aperture 3 odd resolutions)
// Only 1/3 of cells valid - those where (i+j) % 3 == 0
// The modulo arithmetic compacts the sparse grid into dense numbering
static uint64_t cell_index_2d_offset_ap3(long long i, long long j, long long dim) {
    uint64_t idx = i * dim / 3;
    switch (i % 3) {
        case 0: idx += j / 3; break;
        case 1: idx += (j - 2) / 3; break;
        case 2: idx += (j - 1) / 3; break;
    }
    return idx;
}

// Inverse: cell index to (i, j) for offset grid (aperture 3)
static void ij_from_cell_index_offset_ap3(uint64_t idx, long long dim,
                                      long long& i, long long& j) {
    i = (idx * 3) / dim;
    j = (idx * 3) % dim;
    switch (i % 3) {
        case 0: break;
        case 1: j += 2; break;
        case 2: j += 1; break;
    }
}

// Calculate cell count and offset per quad for any aperture
// For apertures 3 and 4: nCells = 10 * aperture^res + 2
// For aperture 7: we use surrogate-based encoding where each quad has
// sur_dim² ID slots. Surrogates are the canonical ap7 cell coordinates
// (one surrogate per cell), avoiding the substrate overcount problem.
static void calc_grid_params(int resolution, int aperture,
                             uint64_t& nCells, uint64_t& offsetPerQuad) {
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }
    if (aperture == 7) {
        long long sur_dim = calc_sur_dim_ap7(resolution);
        offsetPerQuad = static_cast<uint64_t>(sur_dim) * sur_dim;
        nCells = 10 * offsetPerQuad + 2;
    } else {
        nCells = 10;
        for (int r = 0; r < resolution; r++) {
            nCells *= aperture;
        }
        nCells += 2;
        offsetPerQuad = (nCells - 2) / 10;
    }
}

// Grid dimension (per-axis) for a given resolution/aperture. Used by every
// cell-ID <-> (quad,i,j) conversion below to select which calc_max_grid_dim_*
// helper applies.
static long long grid_dim_for_aperture(int resolution, int aperture) {
    if (aperture == 3) {
        return calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        return calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        return calc_sur_dim_ap7(resolution);
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
                            bool use_offset_ap3, bool handle_ap7_south_pole,
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
    } else if (use_offset_ap3) {
        ij_from_cell_index_offset_ap3(idx, dim, i, j);
    } else if (aperture == 7) {
        // Decode to surrogate (i,j stored as surrogates for ap7)
        decode_surrogate_ap7(idx, resolution, dim, i, j);
    } else {
        i = idx / dim;
        j = idx % dim;
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
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        uint64_t offset = 0;
        if (q > 0) {
            offset = 1 + (q - 1) * offsetPerQuad;
        }

        // 2D cell index within quad
        uint64_t bnd2D_idx;
        if (use_offset) {
            bnd2D_idx = cell_index_2d_offset_ap3(ii, jj, dim);
        } else if (aperture == 7) {
            // Input (i,j) are surrogate coords for ap7
            bnd2D_idx = encode_surrogate_ap7(ii, jj, resolution, dim);
        } else {
            bnd2D_idx = cell_index_2d_aligned(ii, jj, dim);
        }

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
    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        hexify::ProjectionResult fwd = hexify::snyder_forward(lon[k], lat[k]);

        uint64_t bnd2D_seq;
        int quad;

        if (aperture == 7) {
            // AP7: Go through substrate (Class III quantize) then
            // convert to surrogate. This ensures only valid ap7 cell
            // positions are used (matching DGGRID's architecture).
            long long sub_i, sub_j;
            hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x,
                                         fwd.icosa_triangle_y,
                                         7, resolution, quad, sub_i, sub_j);

            // Convert substrate → surrogate for storage
            long long sur_i, sur_j;
            hexify::substrate_to_surrogate_ap7(sub_i, sub_j, resolution,
                                               sur_i, sur_j);
            bnd2D_seq = encode_surrogate_ap7(sur_i, sur_j, resolution, dim);
        } else {
            // AP3/AP4: standard substrate-based encoding
            long long i, j;
            hexify::icosa_tri_to_quad_ij(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                         aperture, resolution, quad, i, j);
            if (use_offset_ap3) {
                bnd2D_seq = cell_index_2d_offset_ap3(i, j, dim);
            } else {
                bnd2D_seq = cell_index_2d_aligned(i, j, dim);
            }
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

    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, use_offset_ap3, /*handle_ap7_south_pole=*/true,
                        quad, i, j);

        // Handle polar pentagons directly
        if (quad == 0) {
            lon[k] = 0.0;
            lat[k] = 90.0;
            continue;
        }
        if (quad == 11) {
            lon[k] = 0.0;
            lat[k] = -90.0;
            continue;
        }

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

    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, use_offset_ap3, /*handle_ap7_south_pole=*/true,
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
    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, use_offset_ap3, /*handle_ap7_south_pole=*/false,
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

    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        double qx = quad_x[k];
        double qy = quad_y[k];

        uint64_t bnd2D_seq;
        int out_quad;

        if (aperture == 7) {
            // AP7: Class III quantize → substrate → surrogate
            long long sub_i, sub_j;
            hexify::quad_xy_to_ij(q, qx, qy, 7, resolution, out_quad, sub_i, sub_j);

            long long sur_i, sur_j;
            hexify::substrate_to_surrogate_ap7(sub_i, sub_j, resolution, sur_i, sur_j);
            bnd2D_seq = encode_surrogate_ap7(sur_i, sur_j, resolution, dim);
        } else {
            // AP3/AP4: standard substrate path
            int icosa_triangle_face;
            double icosa_triangle_x, icosa_triangle_y;
            hexify::quad_xy_to_icosa_tri(q, qx, qy, icosa_triangle_face,
                                         icosa_triangle_x, icosa_triangle_y);
            long long i, j;
            hexify::icosa_tri_to_quad_ij(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y,
                                         aperture, resolution, out_quad, i, j);
            if (use_offset_ap3) {
                bnd2D_seq = cell_index_2d_offset_ap3(i, j, dim);
            } else {
                bnd2D_seq = cell_index_2d_aligned(i, j, dim);
            }
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
    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, use_offset_ap3, /*handle_ap7_south_pole=*/false,
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
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

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
                        nCells, use_offset, /*handle_ap7_south_pole=*/false,
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
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

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
                        nCells, use_offset, /*handle_ap7_south_pole=*/false,
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
// Mixed Aperture 3/4 (ISEA43H) Cell ID Conversion
// ============================================================================
// ISEA43H uses aperture 4 for the first 'mixed_aperture_level' resolutions
// and aperture 3 for the remaining resolutions.
//
// Cell count: N = 10 * 4^mixed_level * 3^(res - mixed_level) + 2
// ============================================================================

// Calculate grid dimension for mixed aperture 3/4
// The grid class depends on the count of aperture-3 levels at odd positions
static long long calc_max_grid_dim_ap43(int resolution, int mixed_aperture_level) {
    if (resolution == 0) return 0;

    // Count aperture-3 resolutions and determine grid class
    int ap3_count = 0;
    double scale = 1.0;

    for (int r = 1; r <= resolution; r++) {
        if (r <= mixed_aperture_level) {
            // Aperture 4 level: scale by 2
            scale *= 2.0;
        } else {
            // Aperture 3 level: scale by sqrt(3)
            scale *= 1.7320508075688772935;
            ap3_count++;
        }
    }

    // If odd number of aperture-3 levels, we're in "offset" grid class
    // which requires extra sqrt(3) factor for substrate
    bool use_offset = (ap3_count % 2) == 1;
    if (use_offset) {
        scale *= 1.7320508075688772935;
    }

    return static_cast<long long>(scale + 0.000001) - 1;
}

// Check if mixed aperture grid uses offset (Class II) pattern
static inline bool is_offset_grid_ap43(int resolution, int mixed_aperture_level) {
    int ap3_count = 0;
    for (int r = mixed_aperture_level + 1; r <= resolution; r++) {
        ap3_count++;
    }
    return (ap3_count % 2) == 1;
}

// Calculate cell count and offset per quad for mixed aperture
// Formula: nCells = 10 * 4^mixed_level * 3^(res - mixed_level) + 2
static void calc_grid_params_ap43(int resolution, int mixed_aperture_level,
                                   uint64_t& nCells, uint64_t& offsetPerQuad) {
    if (resolution < hexify::kMinResolution || resolution > hexify::kMaxResolution) {
        Rcpp::stop("resolution must be between %d and %d",
                   hexify::kMinResolution, hexify::kMaxResolution);
    }
    nCells = 10;

    // Apply aperture 4 for first mixed_aperture_level resolutions
    for (int r = 0; r < mixed_aperture_level && r < resolution; r++) {
        nCells *= 4;
    }

    // Apply aperture 3 for remaining resolutions
    for (int r = mixed_aperture_level; r < resolution; r++) {
        nCells *= 3;
    }

    nCells += 2;
    offsetPerQuad = (nCells - 2) / 10;
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap43(NumericVector lon, NumericVector lat,
                                       int resolution, int mixed_aperture_level) {
    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
        stop("cpp_lonlat_to_cell_ap43: mixed_aperture_level must be between 0 and resolution");
    }

    int n = lon.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_ap43(resolution, mixed_aperture_level, nCells, offsetPerQuad);

    // Grid dimension
    long long dim = calc_max_grid_dim_ap43(resolution, mixed_aperture_level) + 1;

    // Check if using offset grid
    bool use_offset = is_offset_grid_ap43(resolution, mixed_aperture_level);

    for (int k = 0; k < n; k++) {
        hexify::ProjectionResult fwd = hexify::snyder_forward(lon[k], lat[k]);

        int quad_pre;
        double quad_x, quad_y;
        hexify::icosa_tri_to_quad_xy(fwd.face, fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                                     quad_pre, quad_x, quad_y);

        int quad;
        long long i, j;
        hexify::quad_xy_to_ij_ap43(quad_pre, quad_x, quad_y, resolution, mixed_aperture_level,
                                    quad, i, j);

        // Calculate cell ID offset within quad
        uint64_t offset = 0;
        if (quad > 0) {
            offset = 1 + (quad - 1) * offsetPerQuad;
        }

        // Calculate 2D cell index based on grid pattern
        uint64_t bnd2D_seq;
        if (use_offset) {
            bnd2D_seq = cell_index_2d_offset_ap3(i, j, dim);
        } else {
            bnd2D_seq = cell_index_2d_aligned(i, j, dim);
        }

        uint64_t cid = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(cid);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_cell_to_lonlat_ap43(NumericVector cell_id, int resolution,
                                   int mixed_aperture_level) {
    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
        stop("cpp_cell_to_lonlat_ap43: mixed_aperture_level must be between 0 and resolution");
    }

    int n = cell_id.size();
    NumericVector lon(n);
    NumericVector lat(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params_ap43(resolution, mixed_aperture_level, nCells, offsetPerQuad);

    // Grid dimension
    long long dim = calc_max_grid_dim_ap43(resolution, mixed_aperture_level) + 1;

    // Check if using offset grid
    bool use_offset = is_offset_grid_ap43(resolution, mixed_aperture_level);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        // aperture=3 here is a placeholder for decode_cell_id's ap7 branch,
        // which mixed-aperture (4/3) grids never take.
        decode_cell_id(cell_id[k], resolution, /*aperture=*/3, dim, offsetPerQuad,
                        nCells, use_offset, /*handle_ap7_south_pole=*/false,
                        quad, i, j);

        // Convert quad IJ to lon/lat via quad_xy -> icosa triangle -> lon/lat
        double quad_x, quad_y;
        hexify::quad_ij_to_xy_ap43(quad, i, j, resolution, mixed_aperture_level, quad_x, quad_y);

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
DataFrame cpp_cell_to_quad_ij_ap43(NumericVector cell_id, int resolution,
                                    int mixed_aperture_level) {
    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
        stop("cpp_cell_to_quad_ij_ap43: mixed_aperture_level must be between 0 and resolution");
    }

    int n = cell_id.size();
    IntegerVector out_quad(n);
    NumericVector out_i(n);
    NumericVector out_j(n);

    uint64_t nCells, offsetPerQuad;
    calc_grid_params_ap43(resolution, mixed_aperture_level, nCells, offsetPerQuad);
    long long dim = calc_max_grid_dim_ap43(resolution, mixed_aperture_level) + 1;
    bool use_offset = is_offset_grid_ap43(resolution, mixed_aperture_level);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        // aperture=3 here is a placeholder for decode_cell_id's ap7 branch,
        // which mixed-aperture (4/3) grids never take.
        decode_cell_id(cell_id[k], resolution, /*aperture=*/3, dim, offsetPerQuad,
                        nCells, use_offset, /*handle_ap7_south_pole=*/false,
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
// Neighbor Finding (v0.7.0)
// ============================================================================

// Helper: convert (quad, i, j) to cell_id
// For ap7: (i,j) are surrogates. For ap3/4: (i,j) are substrates.
static double encode_cell_id(int quad, long long i, long long j,
                              uint64_t offsetPerQuad, long long dim,
                              bool use_offset_ap3, int aperture, int resolution) {
    if (quad == 0 && i == 0 && j == 0) {
        return 1.0;
    }

    uint64_t offset = 1 + (quad - 1) * offsetPerQuad;
    uint64_t bnd2D_seq;

    if (use_offset_ap3) {
        bnd2D_seq = cell_index_2d_offset_ap3(i, j, dim);
    } else if (aperture == 7) {
        bnd2D_seq = encode_surrogate_ap7(i, j, resolution, dim);
    } else {
        bnd2D_seq = cell_index_2d_aligned(i, j, dim);
    }

    return static_cast<double>(offset + bnd2D_seq + 1);
}

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

    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    long long max_ij = hexify::get_max_ij(aperture, resolution);

    // Standard hex neighbor offsets in axial coordinates (ap3, ap4)
    static const long long hex_offsets[6][2] = {
        { 1,  0}, { 0,  1}, {-1,  1},
        {-1,  0}, { 0, -1}, { 1, -1}
    };

    // Aperture 7 surrogate neighbor offsets (hexify coordinate convention)
    static const long long sur_offsets[6][2] = {
        {1, 0}, {1, 1}, {0, 1}, {-1, 0}, {-1, -1}, {0, -1}
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

        if (idx == 0) {
            // North pole: pentagon, skip for now
            out[k] = Rcpp::NumericVector(0);
            continue;
        }

        idx--;
        int quad = static_cast<int>(idx / offsetPerQuad) + 1;
        uint64_t within_quad = idx - (quad - 1) * offsetPerQuad;

        if (aperture == 7) {
            // Aperture 7: direct surrogate-integer approach.
            // 1. Decode cell_id → surrogate (i,j) directly
            // 2. Add ±1 neighbor offset in surrogate space
            // 3. Encode directly (or cross-quad via forward pipeline)
            long long sur_i, sur_j;
            decode_surrogate_ap7(within_quad, resolution, dim, sur_i, sur_j);

            for (int d = 0; d < 6; d++) {
                long long nbr_sur_i = sur_i + sur_offsets[d][0];
                long long nbr_sur_j = sur_j + sur_offsets[d][1];

                // Bounds check: use surrogate offset/dim to detect cross-quad
                long long sur_offset = calc_sur_offset_ap7(resolution);
                long long shifted_i = nbr_sur_i + sur_offset;
                long long shifted_j = nbr_sur_j + sur_offset;
                bool in_bounds = (shifted_i >= 0 && shifted_j >= 0 &&
                                   shifted_i < dim && shifted_j < dim);

                if (in_bounds) {
                    // Same quad: encode directly
                    uint64_t nbr_offset = 1 + (quad - 1) * offsetPerQuad;
                    uint64_t nbr_idx = encode_surrogate_ap7(
                        nbr_sur_i, nbr_sur_j, resolution, dim);
                    neighbor_ids.push_back(
                        static_cast<double>(nbr_offset + nbr_idx + 1));
                } else {
                    // Cross-quad boundary: route through forward pipeline
                    double nbr_qx, nbr_qy;
                    hexify::surrogate_ij_to_quad_xy_ap7(
                        nbr_sur_i, nbr_sur_j, resolution, nbr_qx, nbr_qy);

                    int tri_face;
                    double tri_x, tri_y;
                    if (!hexify::try_quad_xy_to_icosa_tri(quad, nbr_qx, nbr_qy,
                                                           tri_face, tri_x, tri_y)) {
                        continue;
                    }

                    auto ll = hexify::face_xy_to_ll(tri_x, tri_y, tri_face);
                    auto fwd = hexify::snyder_forward(ll.first, ll.second);

                    int final_quad;
                    double final_qx, final_qy;
                    hexify::icosa_tri_to_quad_xy(fwd.face, fwd.icosa_triangle_x,
                                                  fwd.icosa_triangle_y,
                                                  final_quad, final_qx, final_qy);
                    long long final_sur_i, final_sur_j;
                    hexify::quad_xy_to_surrogate_ij_ap7(
                        final_qx, final_qy, resolution, final_sur_i, final_sur_j);

                    neighbor_ids.push_back(
                        encode_cell_id(final_quad, final_sur_i, final_sur_j,
                                       offsetPerQuad, dim, use_offset_ap3,
                                       aperture, resolution));
                }
            }
        } else {
            // Aperture 3 and 4: standard hex offset approach
            long long i, j;
            if (use_offset_ap3) {
                ij_from_cell_index_offset_ap3(within_quad, dim, i, j);
            } else {
                i = within_quad / dim;
                j = within_quad % dim;
            }

            for (int d = 0; d < 6; d++) {
                long long ni = i + hex_offsets[d][0];
                long long nj = j + hex_offsets[d][1];

                bool in_bounds = (ni >= 0 && nj >= 0 && ni <= max_ij && nj <= max_ij);

                if (in_bounds) {
                    neighbor_ids.push_back(
                        encode_cell_id(quad, ni, nj,
                                       offsetPerQuad, dim, use_offset_ap3,
                                       aperture, resolution));
                } else {
                    double nbr_qx, nbr_qy;
                    hexify::quad_ij_to_xy(quad, ni, nj, aperture, resolution,
                                           nbr_qx, nbr_qy);

                    int tri_face;
                    double tri_x, tri_y;
                    if (!hexify::try_quad_xy_to_icosa_tri(quad, nbr_qx, nbr_qy,
                                                           tri_face, tri_x, tri_y)) {
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
                                       offsetPerQuad, dim, use_offset_ap3,
                                       aperture, resolution));
                }
            }
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
    bool use_offset_ap3 = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int quad;
        long long i, j;
        decode_cell_id(cell_id[k], resolution, aperture, dim, offsetPerQuad,
                        nCells, use_offset_ap3, /*handle_ap7_south_pole=*/false,
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
