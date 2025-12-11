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
#include "core_icosa.h"
#include "snyder_forward.h"
#include "snyder_inverse.h"
#include "hex_ap3.h"
#include "hex_ap4.h"
#include "hex_ap7.h"
#include "hex_ap34.h"

using namespace Rcpp;

// ============================================================================
// Aperture 3 - Pure DGGRID (i,j) coordinates
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_quantify_ap3(double tx, double ty, int resolution) {
  long long i, j;
  hexify::hex_quantify_ap3(tx, ty, resolution, i, j);
  return Rcpp::NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_center_ap3(double i, double j, int resolution) {
  double cx, cy;
  hexify::hex_center_ap3(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  return Rcpp::NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
Rcpp::List cpp_hex_corners_ap3(double i, double j, int resolution,
                               double hex_radius = 1.0) {
  double xs[6], ys[6];
  hexify::hex_corners_ap3(static_cast<long long>(i),
                          static_cast<long long>(j),
                          resolution, hex_radius, xs, ys);
  return Rcpp::List::create(
    _["x"] = Rcpp::NumericVector(xs, xs + 6),
    _["y"] = Rcpp::NumericVector(ys, ys + 6)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_lonlat_to_cell_ap3(double lon_deg, double lat_deg,
                                            int resolution) {
  auto fwd = hexify::snyder_fwd(lon_deg, lat_deg);
  long long i, j;
  hexify::hex_quantify_ap3(fwd.tx, fwd.ty, resolution, i, j);
  return Rcpp::NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_cell_to_lonlat_ap3(int face, double i, double j,
                                            int resolution) {
  double cx, cy;
  hexify::hex_center_ap3(static_cast<long long>(i),
                         static_cast<long long>(j),
                         resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return Rcpp::NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Aperture 4 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantify_ap4(double tx, double ty, int resolution) {
  long long i, j;
  hexify::hex_quantify_ap4(tx, ty, resolution, i, j);

  NumericVector result = NumericVector::create(
    Named("i") = (double)i,
    Named("j") = (double)j
  );
  return result;
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap4(double i, double j, int resolution) {
  double cx, cy;
  hexify::hex_center_ap4((long long)i, (long long)j, resolution, cx, cy);

  NumericVector result = NumericVector::create(
    Named("cx") = cx,
    Named("cy") = cy
  );
  return result;
}

// [[Rcpp::export]]
List cpp_hex_corners_ap4(double i, double j, int resolution, double hex_radius = 1.0) {
  double x[6], y[6];
  hexify::hex_corners_ap4((long long)i, (long long)j, resolution, hex_radius, x, y);

  return List::create(
    Named("x") = NumericVector(x, x + 6),
    Named("y") = NumericVector(y, y + 6)
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap4(double lon_deg, double lat_deg, int resolution) {
  auto fwd = hexify::snyder_fwd(lon_deg, lat_deg);
  long long i, j;
  hexify::hex_quantify_ap4(fwd.tx, fwd.ty, resolution, i, j);

  return NumericVector::create(
    Named("face") = fwd.face,
    Named("i") = (double)i,
    Named("j") = (double)j
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap4(int face, double i, double j, int resolution) {
  double cx, cy;
  hexify::hex_center_ap4((long long)i, (long long)j, resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);

  return NumericVector::create(
    Named("lon") = ll.first,
    Named("lat") = ll.second
  );
}

// ============================================================================
// Aperture 7 Bindings
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_hex_quantify_ap7(double tx, double ty, int resolution) {
  long long i, j;
  hexify::hex_quantify_ap7(tx, ty, resolution, i, j);

  return NumericVector::create(
    Named("i") = (double)i,
    Named("j") = (double)j
  );
}

// [[Rcpp::export]]
NumericVector cpp_hex_center_ap7(double i, double j, int resolution) {
  double cx, cy;
  hexify::hex_center_ap7((long long)i, (long long)j, resolution, cx, cy);

  return NumericVector::create(
    Named("cx") = cx,
    Named("cy") = cy
  );
}

// [[Rcpp::export]]
List cpp_hex_corners_ap7(double i, double j, int resolution, double hex_radius = 1.0) {
  double x[6], y[6];
  hexify::hex_corners_ap7((long long)i, (long long)j, resolution, hex_radius, x, y);

  return List::create(
    Named("x") = NumericVector(x, x + 6),
    Named("y") = NumericVector(y, y + 6)
  );
}

// [[Rcpp::export]]
NumericVector cpp_lonlat_to_cell_ap7(double lon_deg, double lat_deg, int resolution) {
  auto fwd = hexify::snyder_fwd(lon_deg, lat_deg);
  long long i, j;
  hexify::hex_quantify_ap7(fwd.tx, fwd.ty, resolution, i, j);

  return NumericVector::create(
    Named("face") = fwd.face,
    Named("i") = (double)i,
    Named("j") = (double)j
  );
}

// [[Rcpp::export]]
NumericVector cpp_cell_to_lonlat_ap7(int face, double i, double j, int resolution) {
  double cx, cy;
  hexify::hex_center_ap7((long long)i, (long long)j, resolution, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);

  return NumericVector::create(
    Named("lon") = ll.first,
    Named("lat") = ll.second
  );
}

// ============================================================================
// Mixed Aperture 3/4
// ============================================================================

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_quantify_ap34(double tx, double ty,
                                          Rcpp::IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i, j;
  hexify::hex_quantify_ap34(tx, ty, seq, i, j);
  return Rcpp::NumericVector::create(
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_center_ap34(double i, double j,
                                        Rcpp::IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx, cy;
  hexify::hex_center_ap34(static_cast<long long>(i),
                          static_cast<long long>(j),
                          seq, cx, cy);
  return Rcpp::NumericVector::create(_["cx"] = cx, _["cy"] = cy);
}

// [[Rcpp::export]]
Rcpp::List cpp_hex_corners_ap34(double i, double j,
                                Rcpp::IntegerVector ap_seq,
                                double hex_radius = 1.0) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double xs[6], ys[6];
  hexify::hex_corners_ap34(static_cast<long long>(i),
                           static_cast<long long>(j),
                           seq, hex_radius, xs, ys);
  return Rcpp::List::create(
    _["x"] = Rcpp::NumericVector(xs, xs + 6),
    _["y"] = Rcpp::NumericVector(ys, ys + 6)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_lonlat_to_cell_ap34(double lon_deg, double lat_deg,
                                             Rcpp::IntegerVector ap_seq) {
  auto fwd = hexify::snyder_fwd(lon_deg, lat_deg);
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i, j;
  hexify::hex_quantify_ap34(fwd.tx, fwd.ty, seq, i, j);
  return Rcpp::NumericVector::create(
    _["face"] = fwd.face,
    _["i"] = static_cast<double>(i),
    _["j"] = static_cast<double>(j)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_cell_to_lonlat_ap34(int face, double i, double j,
                                             Rcpp::IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  double cx, cy;
  hexify::hex_center_ap34(static_cast<long long>(i),
                          static_cast<long long>(j),
                          seq, cx, cy);
  auto ll = hexify::face_xy_to_ll(cx, cy, face);
  return Rcpp::NumericVector::create(_["lon"] = ll.first, _["lat"] = ll.second);
}

// ============================================================================
// Round-trip Test Helpers
// ============================================================================

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap3(double tx, double ty, int resolution) {
  long long i1, j1;
  hexify::hex_quantify_ap3(tx, ty, resolution, i1, j1);
  double cx, cy;
  hexify::hex_center_ap3(i1, j1, resolution, cx, cy);
  long long i2, j2;
  hexify::hex_quantify_ap3(cx, cy, resolution, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap34(double tx, double ty, Rcpp::IntegerVector ap_seq) {
  std::vector<int> seq(ap_seq.begin(), ap_seq.end());
  long long i1, j1;
  hexify::hex_quantify_ap34(tx, ty, seq, i1, j1);
  double cx, cy;
  hexify::hex_center_ap34(i1, j1, seq, cx, cy);
  long long i2, j2;
  hexify::hex_quantify_ap34(cx, cy, seq, i2, j2);
  return (i1 == i2) && (j1 == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap4(double tx, double ty, int resolution) {
  long long i, j;
  hexify::hex_quantify_ap4(tx, ty, resolution, i, j);
  double cx, cy;
  hexify::hex_center_ap4(i, j, resolution, cx, cy);
  long long i2, j2;
  hexify::hex_quantify_ap4(cx, cy, resolution, i2, j2);
  return (i == i2 && j == j2);
}

// [[Rcpp::export]]
bool cpp_test_roundtrip_ap7(double tx, double ty, int resolution) {
  long long i, j;
  hexify::hex_quantify_ap7(tx, ty, resolution, i, j);
  double cx, cy;
  hexify::hex_center_ap7(i, j, resolution, cx, cy);
  long long i2, j2;
  hexify::hex_quantify_ap7(cx, cy, resolution, i2, j2);
  return (i == i2 && j == j2);
}

// [[Rcpp::export]]
Rcpp::List cpp_batch_test_roundtrip_ap3(Rcpp::NumericVector tx_vec,
                                        Rcpp::NumericVector ty_vec,
                                        int resolution) {
  int n = tx_vec.size();
  if (ty_vec.size() != n) Rcpp::stop("tx_vec and ty_vec must have same length");

  Rcpp::LogicalVector success(n);
  Rcpp::NumericVector i_orig(n), j_orig(n), cx(n), cy(n), i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; ++k) {
    long long i1, j1;
    hexify::hex_quantify_ap3(tx_vec[k], ty_vec[k], resolution, i1, j1);
    double cx_k, cy_k;
    hexify::hex_center_ap3(i1, j1, resolution, cx_k, cy_k);
    long long i2, j2;
    hexify::hex_quantify_ap3(cx_k, cy_k, resolution, i2, j2);

    success[k] = (i1 == i2) && (j1 == j2);
    i_orig[k] = static_cast<double>(i1);
    j_orig[k] = static_cast<double>(j1);
    cx[k] = cx_k;
    cy[k] = cy_k;
    i_recomp[k] = static_cast<double>(i2);
    j_recomp[k] = static_cast<double>(j2);
  }

  return Rcpp::List::create(_["success"] = success, _["i_orig"] = i_orig,
    _["j_orig"] = j_orig, _["cx"] = cx, _["cy"] = cy,
    _["i_recomp"] = i_recomp, _["j_recomp"] = j_recomp);
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap4(NumericVector tx_vec, NumericVector ty_vec,
                                   int resolution) {
  int n = tx_vec.size();

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n);
  NumericVector cx(n), cy(n);
  NumericVector i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; k++) {
    double tx = tx_vec[k];
    double ty = ty_vec[k];

    long long i, j;
    hexify::hex_quantify_ap4(tx, ty, resolution, i, j);
    i_orig[k] = (double)i;
    j_orig[k] = (double)j;

    double center_x, center_y;
    hexify::hex_center_ap4(i, j, resolution, center_x, center_y);
    cx[k] = center_x;
    cy[k] = center_y;

    long long i2, j2;
    hexify::hex_quantify_ap4(center_x, center_y, resolution, i2, j2);
    i_recomp[k] = (double)i2;
    j_recomp[k] = (double)j2;

    success[k] = (i == i2 && j == j2);
  }

  return List::create(
    Named("success") = success,
    Named("i_orig") = i_orig,
    Named("j_orig") = j_orig,
    Named("cx") = cx,
    Named("cy") = cy,
    Named("i_recomp") = i_recomp,
    Named("j_recomp") = j_recomp
  );
}

// [[Rcpp::export]]
List cpp_batch_test_roundtrip_ap7(NumericVector tx_vec, NumericVector ty_vec,
                                   int resolution) {
  int n = tx_vec.size();

  LogicalVector success(n);
  NumericVector i_orig(n), j_orig(n);
  NumericVector cx(n), cy(n);
  NumericVector i_recomp(n), j_recomp(n);

  for (int k = 0; k < n; k++) {
    double tx = tx_vec[k];
    double ty = ty_vec[k];

    long long i, j;
    hexify::hex_quantify_ap7(tx, ty, resolution, i, j);
    i_orig[k] = (double)i;
    j_orig[k] = (double)j;

    double center_x, center_y;
    hexify::hex_center_ap7(i, j, resolution, center_x, center_y);
    cx[k] = center_x;
    cy[k] = center_y;

    long long i2, j2;
    hexify::hex_quantify_ap7(center_x, center_y, resolution, i2, j2);
    i_recomp[k] = (double)i2;
    j_recomp[k] = (double)j2;

    success[k] = (i == i2 && j == j2);
  }

  return List::create(
    Named("success") = success,
    Named("i_orig") = i_orig,
    Named("j_orig") = j_orig,
    Named("cx") = cx,
    Named("cy") = cy,
    Named("i_recomp") = i_recomp,
    Named("j_recomp") = j_recomp
  );
}
