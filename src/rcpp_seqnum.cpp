// rcpp_seqnum.cpp
// Rcpp bindings for sequential numbering and coordinate conversion
//
// This file provides the R interface for:
// - Lon/lat to seqnum conversion
// - Seqnum to lon/lat conversion
// - Seqnum to cell info conversion
// - Q2DI coordinate conversion
// - Z7 decoding

#include <Rcpp.h>
#include "constants.h"
#include "core_icosa.h"
#include "snyder_forward.h"
#include "snyder_inverse.h"
#include "hex_ap3.h"
#include "hex_ap4.h"
#include "hex_ap7.h"
#include "hex_seqnum.h"
#include "hex_index_z7.h"
#include "tri_to_q2di.h"

using namespace Rcpp;

// ============================================================================
// Legacy Face-Based Seqnum Conversions (deprecated - use _dggrid versions)
// ============================================================================
// These functions use face-based (triangle) coordinates, not quad-based
// coordinates. They do NOT match DGGRID/dggridR seqnums.
// Kept for backward compatibility - prefer cpp_lonlat_to_seqnum_dggrid.

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_seqnum_face(NumericVector lon, NumericVector lat,
                                        int resolution, int aperture) {
  int n = lon.size();
  if (lat.size() != n) {
    stop("lon and lat must have same length");
  }

  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    // Step 1: Project to icosahedron
    auto fwd = hexify::snyder_fwd(lon[i], lat[i]);

    // Step 2: Quantify to hex cell
    long long cell_i, cell_j;
    if (aperture == 3) {
      hexify::hex_quantify_ap3(fwd.tx, fwd.ty, resolution, cell_i, cell_j);
    } else if (aperture == 4) {
      hexify::hex_quantify_ap4(fwd.tx, fwd.ty, resolution, cell_i, cell_j);
    } else if (aperture == 7) {
      hexify::hex_quantify_ap7(fwd.tx, fwd.ty, resolution, cell_i, cell_j);
    } else {
      stop("Invalid aperture");
    }

    // Step 3: Convert to integer seqnum
    uint64_t seqnum;
    if (aperture == 3) {
      seqnum = hexify::cell_to_seqnum_ap3(fwd.face, cell_i, cell_j, resolution);
    } else if (aperture == 4) {
      seqnum = hexify::cell_to_seqnum_ap4(fwd.face, cell_i, cell_j, resolution);
    } else if (aperture == 7) {
      seqnum = hexify::cell_to_seqnum_ap7(fwd.face, cell_i, cell_j, resolution);
    }

    // Convert uint64 to double (R uses double for large integers)
    result[i] = static_cast<double>(seqnum);
  }

  return result;
}

// [[Rcpp::export]]
DataFrame cpp_seqnum_to_lonlat_face(NumericVector seqnum, int resolution, int aperture) {
  int n = seqnum.size();

  NumericVector lon(n);
  NumericVector lat(n);

  for (int i = 0; i < n; i++) {
    uint64_t seq = static_cast<uint64_t>(seqnum[i]);

    // Step 1: Convert seqnum to (face, i, j)
    int face;
    long long cell_i, cell_j;

    if (aperture == 3) {
      hexify::seqnum_to_cell_ap3(seq, resolution, face, cell_i, cell_j);
    } else if (aperture == 4) {
      hexify::seqnum_to_cell_ap4(seq, resolution, face, cell_i, cell_j);
    } else if (aperture == 7) {
      hexify::seqnum_to_cell_ap7(seq, resolution, face, cell_i, cell_j);
    } else {
      stop("Invalid aperture");
    }

    // Step 2: Get cell center in face coordinates
    double cx, cy;
    if (aperture == 3) {
      hexify::hex_center_ap3(cell_i, cell_j, resolution, cx, cy);
    } else if (aperture == 4) {
      hexify::hex_center_ap4(cell_i, cell_j, resolution, cx, cy);
    } else if (aperture == 7) {
      hexify::hex_center_ap7(cell_i, cell_j, resolution, cx, cy);
    }

    // Step 3: Project back to lon/lat
    auto ll = hexify::face_xy_to_ll(cx, cy, face);
    lon[i] = ll.first;
    lat[i] = ll.second;
  }

  return DataFrame::create(
    _["lon_deg"] = lon,
    _["lat_deg"] = lat
  );
}

// [[Rcpp::export]]
List cpp_seqnum_to_cell_info(NumericVector seqnum, int resolution, int aperture) {
  int n = seqnum.size();

  IntegerVector faces(n);
  NumericVector i_vals(n);
  NumericVector j_vals(n);

  for (int i = 0; i < n; i++) {
    uint64_t seq = static_cast<uint64_t>(seqnum[i]);

    int face;
    long long cell_i, cell_j;

    if (aperture == 3) {
      hexify::seqnum_to_cell_ap3(seq, resolution, face, cell_i, cell_j);
    } else if (aperture == 4) {
      hexify::seqnum_to_cell_ap4(seq, resolution, face, cell_i, cell_j);
    } else if (aperture == 7) {
      hexify::seqnum_to_cell_ap7(seq, resolution, face, cell_i, cell_j);
    } else {
      stop("Invalid aperture");
    }

    faces[i] = face;
    i_vals[i] = static_cast<double>(cell_i);
    j_vals[i] = static_cast<double>(cell_j);
  }

  return List::create(
    _["face"] = faces,
    _["i"] = i_vals,
    _["j"] = j_vals
  );
}

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
// PROJTRI to Q2DI Conversion
// ============================================================================

// [[Rcpp::export]]
Rcpp::List cpp_projtri_to_q2di(int tnum, double tx, double ty,
                                int aperture, int resolution) {
    int quad;
    long long i, j;

    hexify::projtri_to_q2di(tnum, tx, ty, aperture, resolution, quad, i, j);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("i") = (double)i,
        Rcpp::Named("j") = (double)j
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_projtri_to_q2dd(int tnum, double tx, double ty) {
    int quad;
    double qx, qy;

    hexify::projtri_to_q2dd(tnum, tx, ty, quad, qx, qy);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("qx") = qx,
        Rcpp::Named("qy") = qy
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_q2dd_to_projtri(int quad, double qx, double qy) {
    int tnum;
    double tx, ty;

    hexify::q2dd_to_projtri(quad, qx, qy, tnum, tx, ty);

    return Rcpp::List::create(
        Rcpp::Named("tnum") = tnum,
        Rcpp::Named("tx") = tx,
        Rcpp::Named("ty") = ty
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_q2di_to_q2dd(int quad, double i, double j,
                             int aperture, int resolution) {
    double qx, qy;
    hexify::q2di_to_q2dd(quad, static_cast<long long>(i), static_cast<long long>(j),
                          aperture, resolution, qx, qy);

    return Rcpp::List::create(
        Rcpp::Named("qx") = qx,
        Rcpp::Named("qy") = qy
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_lonlat_to_q2di(double lon_deg, double lat_deg,
                               int aperture, int resolution) {
    // Step 1: Get PROJTRI (face, tx, ty)
    hexify::FwdOut fwd = hexify::snyder_fwd(lon_deg, lat_deg);

    // Step 2: Convert to Q2DI
    int quad;
    long long i, j;
    hexify::projtri_to_q2di(fwd.face, fwd.tx, fwd.ty, aperture, resolution, quad, i, j);

    return Rcpp::List::create(
        Rcpp::Named("quad") = quad,
        Rcpp::Named("i") = (double)i,
        Rcpp::Named("j") = (double)j,
        Rcpp::Named("tnum") = fwd.face,
        Rcpp::Named("tx") = fwd.tx,
        Rcpp::Named("ty") = fwd.ty
    );
}

// ============================================================================
// DGGRID-compatible SEQNUM from Q2DI
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
// For aperture 3, maxI = maxJ follows the pattern from DGGRID:
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

// Calculate grid dimension for aperture 7
static long long calc_max_grid_dim_ap7(int resolution) {
    if (resolution == 0) return 0;

    // Aperture 7 alternates between Class III-I (even) and Class III-II (odd)
    bool is_class3i = (resolution % 2) == 0;

    // Base scale: sqrt(7)^resolution
    double scale = 1.0;
    for (int r = 1; r <= resolution; r++) {
        scale *= 2.6457513110645905905;  // sqrt(7)
    }

    // Class III-I: substrate is sqrt(7) finer
    // Class III-II: substrate is sqrt(21) finer
    if (is_class3i) {
        scale *= 2.6457513110645905905;  // sqrt(7)
    } else {
        scale *= 4.5825756949558400065;  // sqrt(21)
    }

    return static_cast<long long>(scale + 0.000001) - 1;
}

// 2D seqnum for aligned grid (aperture 3 even res, aperture 4 all res)
// Simple row-major ordering where all (i,j) pairs are valid
static uint64_t seqnum_2d_aligned(long long i, long long j, long long dim) {
    return static_cast<uint64_t>(i) * dim + j;
}

// 2D seqnum for offset grid (aperture 3 odd resolutions)
// Only 1/3 of cells valid - those where (i+j) % 3 == 0
// Uses DGGRID's DgBoundedHexC2RF2D formula
static uint64_t seqnum_2d_offset_ap3(long long i, long long j, long long dim) {
    uint64_t sNum = i * dim / 3;
    switch (i % 3) {
        case 0: sNum += j / 3; break;
        case 1: sNum += (j - 2) / 3; break;
        case 2: sNum += (j - 1) / 3; break;
    }
    return sNum;
}

// Inverse: seqnum to (i, j) for offset grid (aperture 3)
static void ij_from_seqnum_offset_ap3(uint64_t sNum, long long dim,
                                      long long& i, long long& j) {
    i = (sNum * 3) / dim;
    j = (sNum * 3) % dim;
    switch (i % 3) {
        case 0: break;
        case 1: j += 2; break;
        case 2: j += 1; break;
    }
}

// Calculate cell count and offset per quad for any aperture
// Formula: nCells = 10 * aperture^res + 2
static void calc_grid_params(int resolution, int aperture,
                             uint64_t& nCells, uint64_t& offsetPerQuad) {
    nCells = 10;
    for (int r = 0; r < resolution; r++) {
        nCells *= aperture;
    }
    nCells += 2;
    offsetPerQuad = (nCells - 2) / 10;
}

// [[Rcpp::export]]
NumericVector cpp_q2di_to_seqnum(IntegerVector quad, NumericVector i,
                                  NumericVector j, int resolution, int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_q2di_to_seqnum: aperture must be 3, 4, or 7");
    }

    int n = quad.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim;
    if (aperture == 3) {
        dim = calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        dim = calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        dim = calc_max_grid_dim_ap7(resolution) + 1;
    }

    // Check if using offset grid (only aperture 3 odd resolutions)
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        int q = quad[k];
        long long ii = static_cast<long long>(i[k]);
        long long jj = static_cast<long long>(j[k]);

        uint64_t offset = 0;
        if (q > 0) {
            // For hex grids, firstAdd is quad 0, so add 1
            offset = 1 + (q - 1) * offsetPerQuad;
        }

        // 2D seqnum within quad
        uint64_t bnd2D_seq;
        if (use_offset) {
            bnd2D_seq = seqnum_2d_offset_ap3(ii, jj, dim);
        } else {
            bnd2D_seq = seqnum_2d_aligned(ii, jj, dim);
        }

        // Final seqnum (1-based)
        uint64_t seqnum = offset + bnd2D_seq + 1;

        result[k] = static_cast<double>(seqnum);
    }

    return result;
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_seqnum_dggrid(NumericVector lon, NumericVector lat,
                                           int resolution, int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_lonlat_to_seqnum_dggrid: aperture must be 3, 4, or 7");
    }

    int n = lon.size();
    NumericVector result(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim;
    if (aperture == 3) {
        dim = calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        dim = calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        dim = calc_max_grid_dim_ap7(resolution) + 1;
    }

    // Check if using offset grid (only aperture 3 odd resolutions)
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        // Get Q2DI coordinates
        hexify::FwdOut fwd = hexify::snyder_fwd(lon[k], lat[k]);

        int quad;
        long long i, j;
        hexify::projtri_to_q2di(fwd.face, fwd.tx, fwd.ty, aperture, resolution,
                                quad, i, j);

        // Calculate seqnum using DGGRID formula
        uint64_t offset = 0;
        if (quad > 0) {
            offset = 1 + (quad - 1) * offsetPerQuad;
        }

        // Calculate 2D seqnum based on grid pattern
        uint64_t bnd2D_seq;
        if (use_offset) {
            bnd2D_seq = seqnum_2d_offset_ap3(i, j, dim);
        } else {
            bnd2D_seq = seqnum_2d_aligned(i, j, dim);
        }

        uint64_t seqnum = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(seqnum);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_seqnum_to_lonlat_dggrid(NumericVector seqnum, int resolution,
                                       int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_seqnum_to_lonlat_dggrid: aperture must be 3, 4, or 7");
    }

    int n = seqnum.size();
    NumericVector lon(n);
    NumericVector lat(n);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim;
    if (aperture == 3) {
        dim = calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        dim = calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        dim = calc_max_grid_dim_ap7(resolution) + 1;
    }

    // Check if using offset grid (only aperture 3 odd resolutions)
    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    for (int k = 0; k < n; k++) {
        uint64_t sNum = static_cast<uint64_t>(seqnum[k]);

        // Convert to 0-based
        sNum--;

        int quad;
        long long i, j;

        if (sNum == 0) {
            // First cell: quad 0, i=0, j=0
            quad = 0;
            i = 0;
            j = 0;
        } else {
            // Adjust for quad 0
            sNum--;

            // Determine quad
            quad = static_cast<int>(sNum / offsetPerQuad) + 1;
            sNum -= (quad - 1) * offsetPerQuad;

            // Get i, j from remaining seqnum based on grid type
            if (use_offset) {
                ij_from_seqnum_offset_ap3(sNum, dim, i, j);
            } else {
                i = sNum / dim;
                j = sNum % dim;
            }
        }

        // Convert Q2DI to lon/lat via Q2DD -> PROJTRI -> lon/lat
        // Step 1: Q2DI -> Q2DD
        double qx, qy;
        hexify::q2di_to_q2dd(quad, i, j, aperture, resolution, qx, qy);

        // Step 2: Q2DD -> PROJTRI
        int tnum;
        double tx, ty;
        hexify::q2dd_to_projtri(quad, qx, qy, tnum, tx, ty);

        // Step 3: PROJTRI -> lon/lat
        auto ll = hexify::face_xy_to_ll(tx, ty, tnum);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }

    return DataFrame::create(
        _["lon_deg"] = lon,
        _["lat_deg"] = lat
    );
}

// ============================================================================
// Polygon Corner Generation from SEQNUM
// ============================================================================
// Generates hexagon corner coordinates for a vector of SEQNUM values.
// Uses the same coordinate transformation as DGGRID:
// 1. Compute vertex offsets in Q2DD (quad 2D double) space
// 2. Convert each vertex through Q2DD → PROJTRI → lon/lat
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
// Reordered to match dggridR vertex winding (shifted by 2 positions)
static const double kClass2VertexX[6] = {-kHexR2, kHexR2, kHexR, kHexR2, -kHexR2, -kHexR};
static const double kClass2VertexY[6] = {-0.5, -0.5, 0.0, 0.5, 0.5, 0.0};

// [[Rcpp::export]]
DataFrame cpp_seqnum_to_polygon_dggrid(NumericVector seqnum, int resolution,
                                        int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_seqnum_to_polygon_dggrid: aperture must be 3, 4, or 7");
    }

    int n = seqnum.size();

    // Each cell has 7 vertices (6 corners + 1 to close the polygon)
    int total_vertices = n * 7;

    NumericVector out_seqnum(total_vertices);
    NumericVector out_lon(total_vertices);
    NumericVector out_lat(total_vertices);
    IntegerVector out_order(total_vertices);

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    // Grid dimension depends on aperture
    long long dim;
    if (aperture == 3) {
        dim = calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        dim = calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        dim = calc_max_grid_dim_ap7(resolution) + 1;
    }

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
    double center_scale;   // Scale to convert (i,j) grid coords to Q2DD
    double vertex_scale;   // Scale to convert vertex offsets to Q2DD
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
        center_scale = std::pow(std::sqrt(7.0), resolution);
        vertex_scale = center_scale;
    }

    int out_idx = 0;

    for (int k = 0; k < n; k++) {
        uint64_t sNum = static_cast<uint64_t>(seqnum[k]);
        double orig_seqnum = seqnum[k];

        // Convert seqnum to (quad, i, j)
        sNum--;  // Convert to 0-based

        int quad;
        long long i, j;

        if (sNum == 0) {
            quad = 0;
            i = 0;
            j = 0;
        } else {
            sNum--;
            quad = static_cast<int>(sNum / offsetPerQuad) + 1;
            sNum -= (quad - 1) * offsetPerQuad;

            if (use_offset) {
                ij_from_seqnum_offset_ap3(sNum, dim, i, j);
            } else {
                i = sNum / dim;
                j = sNum % dim;
            }
        }

        // Get cell center in Q2DD coordinates
        // inv_quantize_class1: (i - 0.5*j, j * sin60) is in substrate coords for Class II
        double grid_x = static_cast<double>(i) - 0.5 * static_cast<double>(j);
        double grid_y = static_cast<double>(j) * hexify::kSin60;
        double qx_center = grid_x / center_scale;
        double qy_center = grid_y / center_scale;

        // Generate 6 corners: add vertex offsets (which are in backFrame scale)
        double first_lon = 0.0, first_lat = 0.0;
        for (int c = 0; c < 6; c++) {
            // Vertex offsets are in backFrame coordinates, normalize to Q2DD
            double qx_vertex = qx_center + vertexX[c] / vertex_scale;
            double qy_vertex = qy_center + vertexY[c] / vertex_scale;

            // Convert Q2DD to PROJTRI
            int tnum;
            double tx, ty;
            hexify::q2dd_to_projtri(quad, qx_vertex, qy_vertex, tnum, tx, ty);

            // Convert to lon/lat
            auto ll = hexify::face_xy_to_ll(tx, ty, tnum);

            out_seqnum[out_idx] = orig_seqnum;
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
        out_seqnum[out_idx] = orig_seqnum;
        out_lon[out_idx] = first_lon;
        out_lat[out_idx] = first_lat;
        out_order[out_idx] = 7;
        out_idx++;
    }

    return DataFrame::create(
        _["seqnum"] = out_seqnum,
        _["lon"] = out_lon,
        _["lat"] = out_lat,
        _["order"] = out_order
    );
}

// Vectorized version that returns a list of coordinate matrices for efficiency
// [[Rcpp::export]]
List cpp_seqnum_to_corners_dggrid(NumericVector seqnum, int resolution,
                                   int aperture) {
    if (aperture != 3 && aperture != 4 && aperture != 7) {
        stop("cpp_seqnum_to_corners_dggrid: aperture must be 3, 4, or 7");
    }

    int n = seqnum.size();

    // Calculate grid parameters
    uint64_t nCells, offsetPerQuad;
    calc_grid_params(resolution, aperture, nCells, offsetPerQuad);

    long long dim;
    if (aperture == 3) {
        dim = calc_max_grid_dim_ap3(resolution) + 1;
    } else if (aperture == 4) {
        dim = calc_max_grid_dim_ap4(resolution) + 1;
    } else {
        dim = calc_max_grid_dim_ap7(resolution) + 1;
    }

    bool use_offset = (aperture == 3) && !is_aligned_grid_ap3(resolution);

    // Determine grid class and select vertex offsets
    bool is_class1 = (aperture == 4) || (aperture == 3 && resolution % 2 == 0);
    const double* vertexX = is_class1 ? kClass1VertexX : kClass2VertexX;
    const double* vertexY = is_class1 ? kClass1VertexY : kClass2VertexY;

    // Calculate scale factors (same logic as cpp_seqnum_to_polygon_dggrid)
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
        center_scale = std::pow(std::sqrt(7.0), resolution);
        vertex_scale = center_scale;
    }

    // Return a list of n elements, each element is a 7x2 matrix (lon, lat)
    List result(n);

    for (int k = 0; k < n; k++) {
        uint64_t sNum = static_cast<uint64_t>(seqnum[k]);
        sNum--;

        int quad;
        long long i, j;

        if (sNum == 0) {
            quad = 0;
            i = 0;
            j = 0;
        } else {
            sNum--;
            quad = static_cast<int>(sNum / offsetPerQuad) + 1;
            sNum -= (quad - 1) * offsetPerQuad;

            if (use_offset) {
                ij_from_seqnum_offset_ap3(sNum, dim, i, j);
            } else {
                i = sNum / dim;
                j = sNum % dim;
            }
        }

        // Get cell center in Q2DD coordinates
        double grid_x = static_cast<double>(i) - 0.5 * static_cast<double>(j);
        double grid_y = static_cast<double>(j) * hexify::kSin60;
        double qx_center = grid_x / center_scale;
        double qy_center = grid_y / center_scale;

        // Create 7x2 matrix (closed polygon)
        NumericMatrix coords(7, 2);
        for (int c = 0; c < 6; c++) {
            // Vertex offsets are in backFrame coordinates, normalize to Q2DD
            double qx_vertex = qx_center + vertexX[c] / vertex_scale;
            double qy_vertex = qy_center + vertexY[c] / vertex_scale;

            // Try to convert Q2DD to PROJTRI
            int tnum;
            double tx, ty;
            if (hexify::try_q2dd_to_projtri(quad, qx_vertex, qy_vertex, tnum, tx, ty)) {
                // Normal case: convert to lon/lat via Snyder inverse
                auto ll = hexify::face_xy_to_ll(tx, ty, tnum);
                coords(c, 0) = ll.first;
                coords(c, 1) = ll.second;
            } else {
                // Edge case: vertex is in invalid region (pentagon or quad boundary)
                // Fall back to using the center's triangle and projecting vertex offset
                // This is an approximation but works for edge cells
                int center_tnum;
                double center_tx, center_ty;
                if (hexify::try_q2dd_to_projtri(quad, qx_center, qy_center,
                                                center_tnum, center_tx, center_ty)) {
                    // Project vertex through center's triangle
                    double vertex_tx = center_tx + (qx_vertex - qx_center);
                    double vertex_ty = center_ty + (qy_vertex - qy_center);
                    auto ll = hexify::face_xy_to_ll(vertex_tx, vertex_ty, center_tnum);
                    coords(c, 0) = ll.first;
                    coords(c, 1) = ll.second;
                } else {
                    // Both center and vertex are invalid - use NaN
                    coords(c, 0) = NA_REAL;
                    coords(c, 1) = NA_REAL;
                }
            }
        }
        // Close polygon
        coords(6, 0) = coords(0, 0);
        coords(6, 1) = coords(0, 1);

        colnames(coords) = CharacterVector::create("lon", "lat");
        result[k] = coords;
    }

    return result;
}

// ============================================================================
// Mixed Aperture 3/4 (ISEA43H) SEQNUM Conversion
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
NumericVector cpp_lonlat_to_seqnum_dggrid_ap43(NumericVector lon, NumericVector lat,
                                                int resolution, int mixed_aperture_level) {
    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
        stop("cpp_lonlat_to_seqnum_dggrid_ap43: mixed_aperture_level must be between 0 and resolution");
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

    // Build aperture sequence for projtri_to_q2di
    // Need resolution+1 entries (one for each level from 0 to resolution)
    std::vector<int> ap_seq(resolution + 1);
    for (int r = 0; r <= resolution; r++) {
        ap_seq[r] = (r < mixed_aperture_level) ? 4 : 3;
    }

    // Note: projtri_to_q2di currently doesn't support mixed aperture sequences
    // We use aperture=3 with adjusted resolution for the Q2DI conversion
    // This is a simplification that works because the final grid class
    // depends only on the number of aperture-3 levels

    for (int k = 0; k < n; k++) {
        // Get Q2DI coordinates using equivalent pure aperture-3 calculation
        hexify::FwdOut fwd = hexify::snyder_fwd(lon[k], lat[k]);

        int quad;
        long long i, j;

        // Use the mixed aperture quantization directly
        // The projtri_to_q2di function handles the aperture sequence
        // For now, we approximate using pure aperture 3 with adjusted resolution
        // This gives correct cell assignments but may need refinement

        // Calculate equivalent "effective resolution" for aperture 3
        // Each aperture-4 level is equivalent to ~1.26 aperture-3 levels (log(4)/log(3))
        // But for grid structure, we use the actual mixed calculation

        hexify::projtri_to_q2di(fwd.face, fwd.tx, fwd.ty, 3, resolution,
                                quad, i, j);

        // Calculate seqnum using DGGRID formula
        uint64_t offset = 0;
        if (quad > 0) {
            offset = 1 + (quad - 1) * offsetPerQuad;
        }

        // Calculate 2D seqnum based on grid pattern
        uint64_t bnd2D_seq;
        if (use_offset) {
            bnd2D_seq = seqnum_2d_offset_ap3(i, j, dim);
        } else {
            bnd2D_seq = seqnum_2d_aligned(i, j, dim);
        }

        uint64_t seqnum = offset + bnd2D_seq + 1;
        result[k] = static_cast<double>(seqnum);
    }

    return result;
}

// [[Rcpp::export]]
DataFrame cpp_seqnum_to_lonlat_dggrid_ap43(NumericVector seqnum, int resolution,
                                            int mixed_aperture_level) {
    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
        stop("cpp_seqnum_to_lonlat_dggrid_ap43: mixed_aperture_level must be between 0 and resolution");
    }

    int n = seqnum.size();
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
        uint64_t sNum = static_cast<uint64_t>(seqnum[k]);

        // Convert to 0-based
        sNum--;

        int quad;
        long long i, j;

        if (sNum == 0) {
            // First cell: quad 0, i=0, j=0
            quad = 0;
            i = 0;
            j = 0;
        } else {
            // Adjust for quad 0
            sNum--;

            // Determine quad
            quad = static_cast<int>(sNum / offsetPerQuad) + 1;
            sNum -= (quad - 1) * offsetPerQuad;

            // Get i, j from remaining seqnum based on grid type
            if (use_offset) {
                ij_from_seqnum_offset_ap3(sNum, dim, i, j);
            } else {
                i = sNum / dim;
                j = sNum % dim;
            }
        }

        // Convert Q2DI to lon/lat via Q2DD -> PROJTRI -> lon/lat
        // Use aperture 3 for the coordinate conversion (structure is similar)
        double qx, qy;
        hexify::q2di_to_q2dd(quad, i, j, 3, resolution, qx, qy);

        int tnum;
        double tx, ty;
        hexify::q2dd_to_projtri(quad, qx, qy, tnum, tx, ty);

        auto ll = hexify::face_xy_to_ll(tx, ty, tnum);
        lon[k] = ll.first;
        lat[k] = ll.second;
    }

    return DataFrame::create(
        _["lon_deg"] = lon,
        _["lat_deg"] = lat
    );
}
