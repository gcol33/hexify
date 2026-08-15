// rcpp_aperture.cpp
// Rcpp bindings for aperture-specific quantization functions
//
// This file provides the R interface for:
// - Aperture 3 quantization and center/corner computation
// - Aperture 4 quantization and center/corner computation
// - Aperture 7 quantization and center/corner computation
// - Mixed aperture sequence support
// - Roundtrip testing functions
//
// The ap3/ap4/ap7 exported functions delegate to shared helpers parameterized
// by aperture number, which call the unified hex_quantize/hex_center/hex_corners
// API from aperture.h. The mixed-sequence functions take an aperture vector and
// call the aperture_sequence API.

#include <Rcpp.h>
#include <array>
#include <cmath>
#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "aperture.h"
#include "aperture_sequence.h"

using namespace Rcpp;

// ============================================================================
// Constants
// ============================================================================

constexpr int HEX_VERTICES = 6;

// ============================================================================
// Shared helpers (aperture-parameterized, not exported to R)
// ============================================================================

// Validate an R-supplied i/j grid coordinate before casting to long long.
// Casting NA/NaN/Inf to an integer type is undefined behavior; this turns
// that into a normal R error instead.
static long long checked_ij(double v, const char* label) {
  if (!std::isfinite(v)) {
    Rcpp::stop("'%s' must be a finite number (not NA/NaN/Inf)", label);
  }
  return static_cast<long long>(v);
}

static NumericVector quantize_impl(double icosa_triangle_x, double icosa_triangle_y,
                                   int aperture, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize(icosa_triangle_x, icosa_triangle_y, aperture, resolution, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

static NumericVector center_impl(double i, double j, int aperture, int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center(checked_ij(i, "i"),
                     checked_ij(j, "j"),
                     aperture, resolution, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

static List corners_impl(double i, double j, int aperture, int resolution,
                          double hex_radius) {
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners(checked_ij(i, "i"),
                      checked_ij(j, "j"),
                      aperture, resolution, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

static NumericVector lonlat_to_cell_impl(double lon_deg, double lat_deg,
                                          int aperture, int resolution) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i = 0, j = 0;
  hexify::hex_quantize(fwd.icosa_triangle_x, fwd.icosa_triangle_y,
                       aperture, resolution, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

static NumericVector cell_to_lonlat_impl(int face, double i, double j,
                                          int aperture, int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center(checked_ij(i, "i"),
                     checked_ij(j, "j"),
                     aperture, resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

static bool test_roundtrip_impl(double icosa_triangle_x, double icosa_triangle_y,
                                 int aperture, int resolution) {
  long long i1 = 0, j1 = 0;
  hexify::hex_quantize(icosa_triangle_x, icosa_triangle_y, aperture, resolution, i1, j1);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center(i1, j1, aperture, resolution, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize(cx, cy, aperture, resolution, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

static List batch_test_roundtrip_impl(NumericVector tx_vec, NumericVector ty_vec,
                                       int aperture, int resolution) {
  const int n = tx_vec.size();
  if (ty_vec.size() != n) {
    stop("tx_vec and ty_vec must have same length");
  }

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n), cx(n), cy(n), i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; ++k) {
    long long i1 = 0, j1 = 0;
    hexify::hex_quantize(tx_vec[k], ty_vec[k], aperture, resolution, i1, j1);
    i_orig[k] = static_cast<double>(i1);
    j_orig[k] = static_cast<double>(j1);

    double cx_k = 0.0, cy_k = 0.0;
    hexify::hex_center(i1, j1, aperture, resolution, cx_k, cy_k);
    cx[k] = cx_k;
    cy[k] = cy_k;

    long long i2 = 0, j2 = 0;
    hexify::hex_quantize(cx_k, cy_k, aperture, resolution, i2, j2);
    i_recomp[k] = static_cast<double>(i2);
    j_recomp[k] = static_cast<double>(j2);

    success[k] = (i1 == i2) && (j1 == j2);
  }

  return List::create(
    _["success"] = success,
    _["i_orig"] = i_orig,
    _["j_orig"] = j_orig,
    _["cx"] = cx,
    _["cy"] = cy,
    _["i_recomp"] = i_recomp,
    _["j_recomp"] = j_recomp
  );
}

// ============================================================================
// Aperture 3 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return quantize_impl(icosa_triangle_x, icosa_triangle_y, 3, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap3(double i, double j, int resolution) {
  return center_impl(i, j, 3, resolution);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap3(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  return corners_impl(i, j, 3, resolution, hex_radius);
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap3(double lon_deg, double lat_deg,
                                      int resolution) {
  return lonlat_to_cell_impl(lon_deg, lat_deg, 3, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap3(int face, double i, double j,
                                      int resolution) {
  return cell_to_lonlat_impl(face, i, j, 3, resolution);
}

// ============================================================================
// Aperture 4 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return quantize_impl(icosa_triangle_x, icosa_triangle_y, 4, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap4(double i, double j, int resolution) {
  return center_impl(i, j, 4, resolution);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap4(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  return corners_impl(i, j, 4, resolution, hex_radius);
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap4(double lon_deg, double lat_deg,
                                      int resolution) {
  return lonlat_to_cell_impl(lon_deg, lat_deg, 4, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap4(int face, double i, double j,
                                      int resolution) {
  return cell_to_lonlat_impl(face, i, j, 4, resolution);
}

// ============================================================================
// Aperture 7 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return quantize_impl(icosa_triangle_x, icosa_triangle_y, 7, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap7(double i, double j, int resolution) {
  return center_impl(i, j, 7, resolution);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap7(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  return corners_impl(i, j, 7, resolution, hex_radius);
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap7(double lon_deg, double lat_deg,
                                      int resolution) {
  return lonlat_to_cell_impl(lon_deg, lat_deg, 7, resolution);
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap7(int face, double i, double j,
                                      int resolution) {
  return cell_to_lonlat_impl(face, i, j, 7, resolution);
}

// ============================================================================
// Mixed Aperture Sequence Bindings (separate API, takes an aperture vector)
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_mixed(double icosa_triangle_x, double icosa_triangle_y,
                                     IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i = 0, j = 0;
  hexify::hex_quantize_mixed(icosa_triangle_x, icosa_triangle_y, seq, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_mixed(double i, double j,
                                   IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_mixed(checked_ij(i, "i"),
                           checked_ij(j, "j"),
                           seq, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
List cpp_hex_corners_mixed(double i, double j,
                           IntegerVector ap_seq,
                           double hex_radius = 1.0) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners_mixed(checked_ij(i, "i"),
                            checked_ij(j, "j"),
                            seq, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_mixed(double lon_deg, double lat_deg,
                                       IntegerVector ap_seq) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i = 0, j = 0;
  hexify::hex_quantize_mixed(fwd.icosa_triangle_x, fwd.icosa_triangle_y, seq, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_mixed(int face, double i, double j,
                                       IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_mixed(checked_ij(i, "i"),
                           checked_ij(j, "j"),
                           seq, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Round-trip Test Helpers
// ============================================================================

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return test_roundtrip_impl(icosa_triangle_x, icosa_triangle_y, 3, resolution);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return test_roundtrip_impl(icosa_triangle_x, icosa_triangle_y, 4, resolution);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  return test_roundtrip_impl(icosa_triangle_x, icosa_triangle_y, 7, resolution);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_mixed(double icosa_triangle_x, double icosa_triangle_y, IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i1 = 0, j1 = 0;
  hexify::hex_quantize_mixed(icosa_triangle_x, icosa_triangle_y, seq, i1, j1);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_mixed(i1, j1, seq, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize_mixed(cx, cy, seq, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap3(NumericVector tx_vec, NumericVector ty_vec,
                                  int resolution) {
  return batch_test_roundtrip_impl(tx_vec, ty_vec, 3, resolution);
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap4(NumericVector tx_vec, NumericVector ty_vec,
                                  int resolution) {
  return batch_test_roundtrip_impl(tx_vec, ty_vec, 4, resolution);
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap7(NumericVector tx_vec, NumericVector ty_vec,
                                  int resolution) {
  return batch_test_roundtrip_impl(tx_vec, ty_vec, 7, resolution);
}
