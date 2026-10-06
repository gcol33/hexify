// rcpp_projection.cpp
// Rcpp bindings for the face projections (ISEA and Fuller)
//
// This file provides the R interface for:
// - Icosahedron construction and the faces of each solid
// - Forward projection (lon/lat → face/tx/ty)
// - Inverse projection (face/x/y → lon/lat)
// - Precision control

#include <Rcpp.h>
#include <cmath>
#include "polyhedron.h"
#include "projection_forward.h"
#include "projection_inverse.h"
#include "snyder_triangle.h"
#include "rcpp_icosa.h"

using namespace Rcpp;

// ============================================================================
// Icosahedron Construction
// ============================================================================

// [[Rcpp::export]]
void cpp_build_icosa(double vert0_lon_deg = 11.25,
                     double vert0_lat_deg = 58.282525588538995,
                     double azimuth_deg   = 0.0) {
  hexify::build_icosa_full(vert0_lon_deg, vert0_lat_deg, azimuth_deg);
}

// The solid an icosa argument names: its name, its numbers of faces, vertices
// and diamond quads, whether it carries a hexagonal grid, the valence of each
// vertex (= quad), Snyder's g and G in degrees, and the arc of an edge in
// degrees.
// [[Rcpp::export]]
List cpp_solid_info(NumericVector icosa) {
  activate_icosa(icosa);
  const hexify::SolidTopology& t = hexify::topo();
  IntegerVector valence(t.n_verts);
  for (int v = 0; v < t.n_verts; ++v) valence[v] = t.valence[v];
  // The edge arc is twice the angle between a face centre and an edge
  // midpoint, read from the face's right spherical triangle: cos(arc / 2) =
  // cos(g) / cos(angle from centre to edge midpoint), with
  // tan(centre to midpoint) = tan(g) cos(60 deg).
  const double g = t.snyder.el_angle;
  const double mid = std::atan(std::tan(g) * 0.5);
  const double edge_arc = 2.0 * std::acos(std::cos(g) / std::cos(mid));
  return List::create(_["name"] = std::string(t.name),
                      _["n_faces"] = t.n_faces,
                      _["n_verts"] = t.n_verts,
                      _["n_diamonds"] = t.n_diamonds(),
                      _["has_grid"] = t.has_quads,
                      _["valence"] = valence,
                      _["g_deg"] = hexify::rad2deg(g),
                      _["G_deg"] = hexify::rad2deg(t.snyder.g_angle),
                      _["edge_arc_deg"] = hexify::rad2deg(edge_arc));
}

// [[Rcpp::export]]
int cpp_which_face(NumericVector icosa, double lon_deg, double lat_deg) {
  activate_icosa(icosa);
  return hexify::which_face(lon_deg, lat_deg);
}

// [[Rcpp::export]]
DataFrame cpp_face_centers(NumericVector icosa) {
  activate_icosa(icosa);
  const auto& C = hexify::face_centers();
  const int n = hexify::poly().n_faces();
  NumericVector lon(n), lat(n);
  for (int i = 0; i < n; ++i) {
    lon[i] = C[i].lon;
    lat[i] = C[i].lat;
  }
  return DataFrame::create(_["lon"] = lon, _["lat"] = lat);
}

// ============================================================================
// Forward Projection
// ============================================================================

// [[Rcpp::export]]
NumericVector cpp_icosa_forward(NumericVector icosa, double lon_deg, double lat_deg) {
  activate_icosa(icosa);
  auto out = hexify::snyder_forward(lon_deg, lat_deg);
  return NumericVector::create(_["face"] = out.face,
                               _["icosa_triangle_x"] = out.icosa_triangle_x,
                               _["icosa_triangle_y"] = out.icosa_triangle_y);
}

// [[Rcpp::export]]
NumericVector cpp_project_to_icosa_triangle(NumericVector icosa, int face,
                                            double lon_deg, double lat_deg) {
  activate_icosa(icosa);
  auto xy = hexify::snyder_forward_to_face(face, lon_deg, lat_deg);
  return NumericVector::create(_["icosa_triangle_x"] = xy.first,
                               _["icosa_triangle_y"] = xy.second);
}

// ============================================================================
// Inverse Projection
// ============================================================================

// 'newton' solves the ISEA inverse by Newton's method instead of in closed
// form, to check the closed form against.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_face_xy_to_ll(NumericVector icosa, double x, double y, int face,
                                      bool newton = false) {
  activate_icosa(icosa);
  auto ll = hexify::face_xy_to_ll(x, y, face,
                                  newton ? hexify::InverseSolver::Newton
                                         : hexify::InverseSolver::Closed);
  return Rcpp::NumericVector::create(_["lon"] = ll.first,
                                     _["lat"] = ll.second);
}

// Snyder's map on the spherical triangle (v0, v1, v2): unit vectors (rows of
// 'v') -> barycentric (b1, b2), and back.
// [[Rcpp::export]]
NumericMatrix cpp_snyder_triangle_forward(NumericVector v0, NumericVector v1,
                                          NumericVector v2, NumericMatrix v) {
  const hexify::SnyderTriangle t = hexify::make_snyder_triangle(v0.begin(), v1.begin(), v2.begin());
  NumericMatrix out(v.nrow(), 2);
  for (int i = 0; i < v.nrow(); ++i) {
    const double p[3] = {v(i, 0), v(i, 1), v(i, 2)};
    hexify::snyder_triangle_forward(t, p, out(i, 0), out(i, 1));
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_snyder_triangle_inverse(NumericVector v0, NumericVector v1,
                                          NumericVector v2, NumericMatrix b) {
  const hexify::SnyderTriangle t = hexify::make_snyder_triangle(v0.begin(), v1.begin(), v2.begin());
  NumericMatrix out(b.nrow(), 3);
  for (int i = 0; i < b.nrow(); ++i) {
    double p[3];
    hexify::snyder_triangle_inverse(t, b(i, 0), b(i, 1), p);
    for (int k = 0; k < 3; ++k) out(i, k) = p[k];
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_icosa_face_params(NumericVector icosa, int face) {
  activate_icosa(icosa);
  if (face < 0 || face >= hexify::poly().n_faces()) Rcpp::stop("face out of range for the solid");
  const auto& C = hexify::face_centers();
  return Rcpp::NumericVector::create(
    _["cen_lat"] = C[face].lat,
    _["cen_lon"] = C[face].lon,
    _["face_azimuth_offset"] = hexify::snyder_get_face_azimuth_offset(face)
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_hex_index_face_to_lonlat(NumericVector icosa, double x, double y,
                                                 double cen_lat, double cen_lon,
                                                 double face_azimuth_offset,
                                                 bool degrees = true) {
  activate_icosa(icosa);
  const auto& S = hexify::poly();
  int face = 0;
  double best = 1e300;
  for (int f = 0; f < S.n_faces(); ++f) {
    double d = std::fabs(S.centers[f].lat - cen_lat)
             + std::fabs(S.centers[f].lon - cen_lon)
             + std::fabs(S.face_azimuth_offset[f] - face_azimuth_offset);
    if (d < best) { best = d; face = f; }
  }

  auto ll_deg = hexify::face_xy_to_ll(x, y, face);

  if (!degrees) {
    const double lon_rad = hexify::deg2rad(ll_deg.first);
    const double lat_rad = hexify::deg2rad(ll_deg.second);
    return Rcpp::NumericVector::create(lon_rad, lat_rad);
  }
  return Rcpp::NumericVector::create(ll_deg.first, ll_deg.second);
}
