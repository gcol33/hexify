// rcpp_aperture.cpp
// Rcpp bindings for aperture-specific quantization functions
//
// This file provides the R interface for:
// - Aperture 3 quantization and center/corner computation
// - Aperture 4 quantization and center/corner computation
// - Aperture 7 quantization and center/corner computation
// - Mixed aperture 3/4 support
// - Roundtrip testing functions

#include <Rcpp.h>
#include <array>
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
// Aperture 3 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize_ap3(icosa_triangle_x, icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap3(double i, double j, int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap3(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap3(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners_ap3(static_cast<long long>(i),
                          static_cast<long long>(j),
                          resolution, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap3(double lon_deg, double lat_deg,
                                      int resolution) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i = 0, j = 0;
  hexify::hex_quantize_ap3(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap3(int face, double i, double j,
                                      int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap3(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Aperture 4 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize_ap4(icosa_triangle_x, icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap4(double i, double j, int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap4(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap4(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners_ap4(static_cast<long long>(i),
                          static_cast<long long>(j),
                          resolution, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap4(double lon_deg, double lat_deg,
                                      int resolution) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i = 0, j = 0;
  hexify::hex_quantize_ap4(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap4(int face, double i, double j,
                                      int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap4(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Aperture 7 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize_ap7(icosa_triangle_x, icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap7(double i, double j, int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap7(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap7(double i, double j, int resolution,
                         double hex_radius = 1.0) {
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners_ap7(static_cast<long long>(i),
                          static_cast<long long>(j),
                          resolution, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap7(double lon_deg, double lat_deg,
                                      int resolution) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  long long i = 0, j = 0;
  hexify::hex_quantize_ap7(fwd.icosa_triangle_x, fwd.icosa_triangle_y, resolution, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap7(int face, double i, double j,
                                      int resolution) {
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap7(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Mixed Aperture 3/4 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantize_ap34(double icosa_triangle_x, double icosa_triangle_y,
                                    IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i = 0, j = 0;
  hexify::hex_quantize_ap34(icosa_triangle_x, icosa_triangle_y, seq, i, j);
  return NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap34(double i, double j,
                                  IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap34(static_cast<long long>(i),
                          static_cast<long long>(j),
                          seq, cx, cy);
  return NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
List cpp_hex_corners_ap34(double i, double j,
                          IntegerVector ap_seq,
                          double hex_radius = 1.0) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  std::array<double, HEX_VERTICES> xs{};
  std::array<double, HEX_VERTICES> ys{};
  hexify::hex_corners_ap34(static_cast<long long>(i),
                           static_cast<long long>(j),
                           seq, hex_radius, xs.data(), ys.data());
  return List::create(
    _["x"] = NumericVector(xs.begin(), xs.end()),
    _["y"] = NumericVector(ys.begin(), ys.end())
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap34(double lon_deg, double lat_deg,
                                       IntegerVector ap_seq) {
  auto fwd = hexify::snyder_forward(lon_deg, lat_deg);
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i = 0, j = 0;
  hexify::hex_quantize_ap34(fwd.icosa_triangle_x, fwd.icosa_triangle_y, seq, i, j);
  return NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap34(int face, double i, double j,
                                       IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap34(static_cast<long long>(i),
                          static_cast<long long>(j),
                          seq, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Round-trip Test Helpers
// ============================================================================

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i1 = 0, j1 = 0;
  hexify::hex_quantize_ap3(icosa_triangle_x, icosa_triangle_y, resolution, i1, j1);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap3(i1, j1, resolution, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize_ap3(cx, cy, resolution, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap34(double icosa_triangle_x, double icosa_triangle_y, IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i1 = 0, j1 = 0;
  hexify::hex_quantize_ap34(icosa_triangle_x, icosa_triangle_y, seq, i1, j1);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap34(i1, j1, seq, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize_ap34(cx, cy, seq, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize_ap4(icosa_triangle_x, icosa_triangle_y, resolution, i, j);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap4(i, j, resolution, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize_ap4(cx, cy, resolution, i2, j2);
  return (i == i2) && (j == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution) {
  long long i = 0, j = 0;
  hexify::hex_quantize_ap7(icosa_triangle_x, icosa_triangle_y, resolution, i, j);
  double cx = 0.0, cy = 0.0;
  hexify::hex_center_ap7(i, j, resolution, cx, cy);
  long long i2 = 0, j2 = 0;
  hexify::hex_quantize_ap7(cx, cy, resolution, i2, j2);
  return (i == i2) && (j == j2);
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap3(NumericVector tx_vec,
                                  NumericVector ty_vec,
                                  int resolution) {
  const int n = tx_vec.size();
  if (ty_vec.size() != n) {
    stop("tx_vec and ty_vec must have same length");
  }

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n), cx(n), cy(n), i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; ++k) {
    long long i1 = 0, j1 = 0;
    hexify::hex_quantize_ap3(tx_vec[k], ty_vec[k], resolution, i1, j1);

    double cx_k = 0.0, cy_k = 0.0;
    hexify::hex_center_ap3(i1, j1, resolution, cx_k, cy_k);

    long long i2 = 0, j2 = 0;
    hexify::hex_quantize_ap3(cx_k, cy_k, resolution, i2, j2);

    success[k] = (i1 == i2) && (j1 == j2);
    i_orig[k] = static_cast<double>(i1);
    j_orig[k] = static_cast<double>(j1);
    cx[k] = cx_k;
    cy[k] = cy_k;
    i_recomp[k] = static_cast<double>(i2);
    j_recomp[k] = static_cast<double>(j2);
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

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap4(NumericVector tx_vec,
                                  NumericVector ty_vec,
                                  int resolution) {
  const int n = tx_vec.size();
  if (ty_vec.size() != n) {
    stop("tx_vec and ty_vec must have same length");
  }

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n), cx(n), cy(n), i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; ++k) {
    long long i = 0, j = 0;
    hexify::hex_quantize_ap4(tx_vec[k], ty_vec[k], resolution, i, j);
    i_orig[k] = static_cast<double>(i);
    j_orig[k] = static_cast<double>(j);

    double center_x = 0.0, center_y = 0.0;
    hexify::hex_center_ap4(i, j, resolution, center_x, center_y);
    cx[k] = center_x;
    cy[k] = center_y;

    long long i2 = 0, j2 = 0;
    hexify::hex_quantize_ap4(center_x, center_y, resolution, i2, j2);
    i_recomp[k] = static_cast<double>(i2);
    j_recomp[k] = static_cast<double>(j2);

    success[k] = (i == i2) && (j == j2);
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

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap7(NumericVector tx_vec,
                                  NumericVector ty_vec,
                                  int resolution) {
  const int n = tx_vec.size();
  if (ty_vec.size() != n) {
    stop("tx_vec and ty_vec must have same length");
  }

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n), cx(n), cy(n), i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; ++k) {
    long long i = 0, j = 0;
    hexify::hex_quantize_ap7(tx_vec[k], ty_vec[k], resolution, i, j);
    i_orig[k] = static_cast<double>(i);
    j_orig[k] = static_cast<double>(j);

    double center_x = 0.0, center_y = 0.0;
    hexify::hex_center_ap7(i, j, resolution, center_x, center_y);
    cx[k] = center_x;
    cy[k] = center_y;

    long long i2 = 0, j2 = 0;
    hexify::hex_quantize_ap7(center_x, center_y, resolution, i2, j2);
    i_recomp[k] = static_cast<double>(i2);
    j_recomp[k] = static_cast<double>(j2);

    success[k] = (i == i2) && (j == j2);
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
