// rcpp_projection.cpp
// Rcpp bindings for Snyder projection functions
//
// This file provides the R interface for:
// - Icosahedron construction
// - Forward projection (lon/lat → face/tx/ty)
// - Inverse projection (face/x/y → lon/lat)
// - Precision control

#include <Rcpp.h>
#include <cmath>
#include "icosahedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"

using namespace Rcpp;

// ============================================================================
// Icosahedron Construction
// ============================================================================

// [[Rcpp::export]]
void cpp_build_icosa(double vert0_lon_deg = 11.25,
                     double vert0_lat_deg = 58.28252559,
                     double azimuth_deg   = 0.0) {
  hexify::build_icosa_full(vert0_lon_deg, vert0_lat_deg, azimuth_deg);
}

// [[Rcpp::export]]
int cpp_which_face(double lon_deg, double lat_deg) {
  return hexify::which_face(lon_deg, lat_deg);
}

// [[Rcpp::export]]
DataFrame cpp_face_centers() {
  const auto& C = hexify::face_centers();
  NumericVector lon(20), lat(20);
  for (int i = 0; i < 20; ++i) {
    lon[i] = C[i].lon;
    lat[i] = C[i].lat;
  }
  return DataFrame::create(_["lon"] = lon, _["lat"] = lat);
}

// ============================================================================
// Forward Projection
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_snyder_forward(double lon_deg, double lat_deg) {
  auto out = hexify::snyder_forward(lon_deg, lat_deg);
  return NumericVector::create(_["face"] = out.face,
                               _["icosa_triangle_x"] = out.icosa_triangle_x,
                               _["icosa_triangle_y"] = out.icosa_triangle_y);
}

// [[Rcpp::export]]
NumericVector cpp_project_to_icosa_triangle(int face, double lon_deg, double lat_deg) {
  auto xy = hexify::snyder_forward_to_face(face, lon_deg, lat_deg);
  return NumericVector::create(_["icosa_triangle_x"] = xy.first,
                               _["icosa_triangle_y"] = xy.second);
}

// ============================================================================
// Inverse Projection
// ============================================================================

// [[Rcpp::export]]
void cpp_snyder_inv_set_precision(std::string mode = "",
                                  Rcpp::Nullable<double> tol = R_NilValue,
                                  Rcpp::Nullable<int>    max_iters = R_NilValue) {
  double tol_v = -1.0;
  int    mi_v  = -1;
  if (tol.isNotNull())       tol_v = Rcpp::as<double>(tol);
  if (max_iters.isNotNull()) mi_v  = Rcpp::as<int>(max_iters);
  hexify::snyder_inv_set_precision(mode, tol_v, mi_v);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_snyder_inv_get_precision() {
  auto pr = hexify::snyder_inv_get_precision();
  return Rcpp::NumericVector::create(_["tol"] = pr.first,
                                     _["max_iters"] = pr.second);
}

// [[Rcpp::export]]
void cpp_snyder_inv_set_verbose(bool v = true) {
  hexify::snyder_inv_set_verbose(v);
}

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_snyder_inv_get_stats_and_reset() {
  auto t = hexify::snyder_inv_get_stats_and_reset();
  return Rcpp::IntegerVector::create(
    _["calls"]       = std::get<0>(t),
    _["iters_total"] = std::get<1>(t),
    _["iters_max"]   = std::get<2>(t),
    _["capped"]      = std::get<3>(t)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_face_xy_to_ll(double x, double y, int face,
                                      Rcpp::Nullable<double> tol = R_NilValue,
                                      Rcpp::Nullable<int>    max_iters = R_NilValue) {
  double tol_v = -1.0;
  int    mi_v  = -1;
  if (tol.isNotNull())       tol_v = Rcpp::as<double>(tol);
  if (max_iters.isNotNull()) mi_v  = Rcpp::as<int>(max_iters);
  auto ll = hexify::face_xy_to_ll(x, y, face, tol_v, mi_v);
  return Rcpp::NumericVector::create(_["lon"] = ll.first,
                                     _["lat"] = ll.second);
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_icosa_face_params(int face) {
  if (face < 0 || face >= 20) Rcpp::stop("face must be 0..19");
  const auto& C = hexify::face_centers();
  return Rcpp::NumericVector::create(
    _["cen_lat"] = C[face].lat,
    _["cen_lon"] = C[face].lon,
    _["face_azimuth_offset"] = hexify::snyder_get_face_azimuth_offset(face)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_index_face_to_lonlat(double x, double y,
                                                 double cen_lat, double cen_lon,
                                                 double face_azimuth_offset,
                                                 bool degrees = true,
                                                 Rcpp::Nullable<double> tol = R_NilValue,
                                                 Rcpp::Nullable<int>    max_iters = R_NilValue) {
  const auto& S = hexify::ico();
  int face = 0;
  double best = 1e300;
  for (int f = 0; f < 20; ++f) {
    double d = std::fabs(S.centers[f].lat - cen_lat)
             + std::fabs(S.centers[f].lon - cen_lon)
             + std::fabs(S.face_azimuth_offset[f] - face_azimuth_offset);
    if (d < best) { best = d; face = f; }
  }

  double tol_v = -1.0;
  int    mi_v  = -1;
  if (tol.isNotNull())       tol_v = Rcpp::as<double>(tol);
  if (max_iters.isNotNull()) mi_v  = Rcpp::as<int>(max_iters);

  auto ll_deg = hexify::face_xy_to_ll(x, y, face, tol_v, mi_v);

  if (!degrees) {
    const double lon_rad = hexify::deg2rad(ll_deg.first);
    const double lat_rad = hexify::deg2rad(ll_deg.second);
    return Rcpp::NumericVector::create(lon_rad, lat_rad);
  }
  return Rcpp::NumericVector::create(ll_deg.first, ll_deg.second);
}
