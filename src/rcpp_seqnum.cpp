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
// Integer Seqnum Conversions
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_seqnum(NumericVector lon, NumericVector lat,
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
DataFrame cpp_seqnum_to_lonlat(NumericVector seqnum, int resolution, int aperture) {
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

// Calculate grid dimension for aperture 3
static long long calc_max_grid_dim_ap3(int resolution) {
    if (resolution == 0) return 0;

    bool aligned = is_aligned_grid_ap3(resolution);

    // Scale factor = sqrt(3)^resolution
    double scale = 1.0;
    for (int r = 1; r <= resolution; r++) {
        scale *= 1.7320508075688772935;  // sqrt(3)
    }

    // Offset grids need extra sqrt(3) factor
    if (!aligned) {
        scale *= 1.7320508075688772935;
    }

    return static_cast<long long>(scale + 0.000001) - 1;
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
