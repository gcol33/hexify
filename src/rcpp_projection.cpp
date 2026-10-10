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
#include "authalic.h"

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
    lat[i] = hexify::to_geodetic_lat(C[i].lat);
  }
  return DataFrame::create(_["lon"] = lon, _["lat"] = lat);
}

// ============================================================================
// Earth model
// ============================================================================

// Latitudes in degrees taken from geodetic to the sphere (authalic latitude)
// on the icosa argument's ellipsoid, or back with 'inverse'; unchanged on the
// sphere. 'exact' evaluates the closed forms the conversions are fitted to.
// [[Rcpp::export]]
NumericVector cpp_sphere_latitude(NumericVector icosa, NumericVector lat,
                                  bool inverse = false, bool exact = false) {
  activate_icosa(icosa);
  const double f = hexify::active_flattening();
  const R_xlen_t n = lat.size();
  NumericVector out(n);
  for (R_xlen_t k = 0; k < n; ++k) {
    const double x = lat[k];
    if (ISNAN(x) || f == 0.0) {
      out[k] = x;
    } else if (exact) {
      const double r = hexify::deg2rad(x);
      out[k] = hexify::rad2deg(inverse ? hexify::geodetic_lat_exact(r, f)
                                       : hexify::authalic_lat_exact(r, f));
    } else {
      out[k] = inverse ? hexify::to_geodetic_lat_deg(x) : hexify::to_sphere_lat_deg(x);
    }
  }
  return out;
}

// The authalic radius of an ellipsoid over its semi-major axis.
// [[Rcpp::export]]
double cpp_authalic_radius_ratio(double flattening) {
  if (!(flattening >= 0.0 && flattening < 1.0)) stop("flattening must lie in [0, 1)");
  return hexify::authalic_radius_ratio(flattening);
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

// 'newton' solves the inverse (Snyder's or Fuller's) by Newton's method
// instead of in closed form, to check the closed form against.
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

// ============================================================================
// Distortion
// ============================================================================

// The derivative of the map from the active ellipsoid onto its authalic
// sphere at a point `geo` of the sphere, composed into the face projection's
// derivative s. A tangent vector of the ellipsoid with ground components
// (north, east) reaches the sphere as (north / k, east * k), k the parallel
// scale; s reads vectors in the frame (away from the face centre, a quarter
// turn clockwise), whose first axis lies at azimuth alpha.
static void compose_authalic_scale(const hexify::Geo& geo, int face,
                                   hexify::FaceScale& s) {
  const double f = hexify::active_flattening();
  if (f == 0.0) return;
  const hexify::PolyData& P = hexify::poly();
  const double k = hexify::authalic_parallel_scale(hexify::to_geodetic_lat(geo.lat), f);
  const double dlon = geo.lon - P.center_lon[face];
  const double y = std::sin(dlon) * P.center_coslat[face];
  const double x = P.center_coslat[face] * std::sin(geo.lat) * std::cos(dlon) -
                   P.center_sinlat[face] * std::cos(geo.lat);
  // Azimuth at geo of the great circle from the face centre through it
  const double alpha = std::atan2(y, x);
  const double ca = std::cos(alpha), sa = std::sin(alpha);
  // m maps (north, east) on the ellipsoid to the face-projection frame
  const double m[2][2] = {{ca / k, sa * k}, {-sa / k, ca * k}};
  hexify::FaceScale out;
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 2; ++c) out.j[r][c] = s.j[r][0] * m[0][c] + s.j[r][1] * m[1][c];
  }
  s = out;
}

// The step of Lambert's construction an R call names: 0 the Lambert point,
// 1 the nudged point, 2 the face projection (hexify::ConstructionStage).
static hexify::ConstructionStage construction_stage(int code) {
  if (code < 0 || code > 2) stop("stage must be 0 (Lambert), 1 (nudge) or 2 (face)");
  return static_cast<hexify::ConstructionStage>(code);
}

// Tissot's indicatrix of the map onto the plane of `stage` at a point of a
// face: the derivative j (face_scale) is a turn by beta, a stretch by (a, b)
// and a turn, so a small circle on the sphere maps to an ellipse with
// semi-axes a >= b along the plane directions beta and beta + 90 degrees. On
// an ellipsoid the circle is drawn on the ellipsoid.
static void tissot_row(const hexify::Geo& geo, int face, hexify::ConstructionStage stage,
                       double& a, double& b, double& beta) {
  hexify::FaceScale s = hexify::face_scale(geo, face, stage);
  compose_authalic_scale(geo, face, s);
  const double e = 0.5 * (s.j[0][0] + s.j[1][1]);
  const double f = 0.5 * (s.j[0][0] - s.j[1][1]);
  const double g = 0.5 * (s.j[1][0] + s.j[0][1]);
  const double h = 0.5 * (s.j[1][0] - s.j[0][1]);
  const double q = std::hypot(e, h);
  const double r = std::hypot(f, g);
  a = q + r;
  b = std::fabs(q - r);
  beta = 0.5 * (std::atan2(g, f) + std::atan2(h, e));
}

static DataFrame tissot_frame(const IntegerVector& face, const NumericVector& a,
                              const NumericVector& b, const NumericVector& beta) {
  return DataFrame::create(_["face"] = face, _["a"] = a, _["b"] = b,
                           _["angle"] = beta);
}

// Tissot's indicatrix at points given in lon/lat, each read on the face it
// lies on (or on 'face' where that is not NA): the scale factors a >= b and
// the plane direction of a, in radians from the face's x axis, of the map
// onto the plane of construction step 'stage' (construction_stage()).
// [[Rcpp::export]]
DataFrame cpp_lonlat_tissot(NumericVector icosa, NumericVector lon,
                            NumericVector lat, IntegerVector face, int stage = 2) {
  activate_icosa(icosa);
  const hexify::ConstructionStage st = construction_stage(stage);
  const R_xlen_t n = lon.size();
  IntegerVector f(n);
  NumericVector a(n), b(n), beta(n);
  for (R_xlen_t k = 0; k < n; ++k) {
    f[k] = face[k] == NA_INTEGER ? hexify::which_face(lon[k], lat[k]) : face[k];
    if (f[k] < 0 || f[k] >= hexify::poly().n_faces()) stop("face out of range for the solid");
    const hexify::Geo g(hexify::deg2rad(lon[k]), hexify::to_sphere_lat(hexify::deg2rad(lat[k])));
    tissot_row(g, f[k], st, a[k], b[k], beta[k]);
  }
  return tissot_frame(f, a, b, beta);
}

// The same at points of one face given in its triangle coordinates.
// [[Rcpp::export]]
DataFrame cpp_face_tri_tissot(NumericVector icosa, int face, NumericVector tx,
                              NumericVector ty, int stage = 2) {
  activate_icosa(icosa);
  const hexify::ConstructionStage st = construction_stage(stage);
  if (face < 0 || face >= hexify::poly().n_faces()) stop("face out of range for the solid");
  const R_xlen_t n = tx.size();
  IntegerVector f(n, face);
  NumericVector a(n), b(n), beta(n);
  for (R_xlen_t k = 0; k < n; ++k) {
    const auto ll = hexify::face_xy_to_sphere_ll(tx[k], ty[k], face);
    const hexify::Geo g(hexify::deg2rad(ll.first), hexify::deg2rad(ll.second));
    tissot_row(g, face, st, a[k], b[k], beta[k]);
  }
  return tissot_frame(f, a, b, beta);
}

// ============================================================================
// Lambert's construction of Snyder's projection
// ============================================================================

// Points given in lon/lat through Lambert's construction of Snyder's
// projection (hexify::construction_point()), each onto the face the forward
// projection puts it on, or onto 'face' where that is not NA. Angles in
// radians, lengths in units of the unit sphere.
// [[Rcpp::export]]
List cpp_lonlat_construction(NumericVector icosa, NumericVector lon,
                                  NumericVector lat, IntegerVector face) {
  activate_icosa(icosa);
  const R_xlen_t n = lon.size();
  IntegerVector f(n);
  NumericVector z(n), az(n), lu(n), lv(n), azp(n), fs(n), pu(n), pv(n), tx(n), ty(n);
  for (R_xlen_t k = 0; k < n; ++k) {
    f[k] = face[k] == NA_INTEGER ? hexify::snyder_forward(lon[k], lat[k]).face : face[k];
    if (f[k] < 0 || f[k] >= hexify::poly().n_faces()) stop("face out of range for the solid");
    const hexify::Geo g(hexify::deg2rad(lon[k]), hexify::to_sphere_lat(hexify::deg2rad(lat[k])));
    const hexify::ConstructionPoint c = hexify::construction_point(g, f[k]);
    z[k] = c.z;
    az[k] = c.az;
    lu[k] = c.lambert_u;
    lv[k] = c.lambert_v;
    azp[k] = c.az_prime;
    fs[k] = c.f;
    pu[k] = c.plane_u;
    pv[k] = c.plane_v;
    tx[k] = c.tx;
    ty[k] = c.ty;
  }
  return List::create(_["face"] = f, _["z"] = z, _["az"] = az,
                      _["lambert_u"] = lu, _["lambert_v"] = lv,
                      _["az_prime"] = azp, _["f"] = fs,
                      _["plane_u"] = pu, _["plane_v"] = pv,
                      _["tx"] = tx, _["ty"] = ty,
                      _["r1"] = hexify::topo().snyder.r1);
}
