#include "projection_inverse.h"
#include "projection_forward.h"
#include "projection_fuller.h"
#include "projection_ivea.h"
#include "polyhedron.h"
#include "snyder_triangle.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace {

using hexify::kPi;
using hexify::kTwoPi;
using hexify::k2PiOver3;
using hexify::kEpsBranch;
using hexify::safe_denom;
using hexify::SnyderParams;

inline double wrap_lon_rad(double L) {
  double t = std::fmod(L + kPi, kTwoPi);
  if (t < 0) t += kTwoPi;
  return t - kPi;
}

// cos and sin of k * 120 degrees, k = 0, 1, 2
constexpr double kSectorCos[3] = {1.0, -0.5, -0.5};
constexpr double kSectorSin[3] = {0.0, 0.86602540378443864676, -0.86602540378443864676};

// Snyder's inverse on a face in closed form: face-plane (x, y) -> distance z
// from the face centre and azimuth from the face's first vertex. The point's
// 120-degree sector is the plane triangle (centre, vertex k, vertex k + 1);
// turned back by k * 120 degrees it is the sector triangle of SnyderParams,
// and its barycentric coordinates there give the point on the sphere.
std::pair<double,double> snyder_face_polar(const SnyderParams& sp, double x, double y) {
  const double X = x - 0.5;
  const double Y = y - sp.origin_y_off / sp.edge;

  // Snyder's azimuth reads atan2(x, y)
  double az = std::atan2(X, Y);
  if (az < 0.0) az += kTwoPi;
  const int k = hexify::azimuth_sector(az);

  // Turned into sector 0, whose vertices sit at (0, 1/sqrt(3)) and
  // (1/2, -1/(2 sqrt(3))) about the centre.
  const double Xr = X * kSectorCos[k] - Y * kSectorSin[k];
  const double Yr = Y * kSectorCos[k] + X * kSectorSin[k];
  const double b2 = 2.0 * Xr;
  const double b1 = std::sqrt(3.0) * Yr + Xr;

  // The point lies on the arc from the face centre through p, so it has p's
  // azimuth.
  double p[3];
  const double z = hexify::snyder_triangle_inverse_ray(sp.sector, b1, b2, p);
  return {z, std::atan2(p[0], p[1]) + k * k2PiOver3};
}

// Snyder's inverse by Newton's method on the auxiliary azimuth, run to a fixed
// tolerance: the check on the closed form.
constexpr double kNewtonTol = 1e-14;
constexpr int kNewtonMaxIters = 100;

double solve_snyder_azimuth(const SnyderParams& sp, double azimuth_initial) {
  if (std::abs(azimuth_initial) <= kEpsBranch) return 0.0;

  const double agh = (sp.r1_squared * sp.tan_el * sp.tan_el) /
                     (2.0 * (1.0 / std::tan(azimuth_initial) + sp.cot_30));

  // f(az) = agh - az - G + (pi - h), h = acos(sin(az) sin(G) cos(g) - cos(az) cos(G))
  double azimuth = azimuth_initial;
  for (int iter = 0; iter < kNewtonMaxIters; ++iter) {
    const double s = std::sin(azimuth);
    const double c = std::cos(azimuth);
    const double h = std::acos(hexify::clampd(s * sp.sin_g * sp.cos_el - c * sp.cos_g, -1.0, 1.0));
    const double residual = agh - azimuth - sp.g_angle + (kPi - h);
    const double derivative = (c * sp.sin_g * sp.cos_el + s * sp.cos_g) / safe_denom(std::sin(h)) - 1.0;
    const double delta = -residual / derivative;
    azimuth += delta;
    if (std::abs(delta) <= kNewtonTol) break;
  }
  return azimuth;
}

std::pair<double,double> snyder_face_polar_newton(const SnyderParams& sp, double x, double y) {
  const double px = x * sp.edge - sp.origin_x_off;
  const double py = y * sp.edge - sp.origin_y_off;
  const double rho = std::hypot(px, py);

  double azimuth_transformed = std::atan2(px, py);
  if (azimuth_transformed < 0.0) azimuth_transformed += kTwoPi;
  const int sector = hexify::azimuth_sector(azimuth_transformed);
  azimuth_transformed -= sector * k2PiOver3;

  const double azimuth = solve_snyder_azimuth(sp, azimuth_transformed);

  // z from rho through Snyder's dz and f
  const double dz_angle = std::atan2(sp.tan_el, std::cos(azimuth) + sp.cot_30 * std::sin(azimuth));
  const double denom = safe_denom(std::cos(azimuth_transformed) + sp.cot_30 * std::sin(azimuth_transformed));
  const double sin_half_dz = safe_denom(std::sin(dz_angle / 2.0));
  const double f_scale = sp.tan_el / (2.0 * denom * sin_half_dz);
  const double arg = hexify::clampd(rho / (2.0 * sp.r1 * f_scale), -1.0, 1.0);
  const double z = 2.0 * std::asin(arg);

  return {z, azimuth + sector * k2PiOver3};
}

} // anon

namespace hexify {

std::pair<double,double> face_xy_to_ll(double x, double y, int face, InverseSolver solver)
{
  const PolyData& P = poly();
  if (face < 0 || face >= P.n_faces()) throw std::runtime_error("face out of range for the solid");
  const SnyderParams& sp = P.topo->snyder;

  // Face centers are in radians
  const auto& C = P.centers;
  const double center_lon = C[face].lon;
  const double center_lat = C[face].lat;
  const double center_sinlat = std::sin(center_lat);
  const double center_coslat = std::cos(center_lat);

  // Exact face center shortcut
  if (std::abs(x * sp.edge - sp.origin_x_off) < kEpsBranch &&
      std::abs(y * sp.edge - sp.origin_y_off) < kEpsBranch) {
    return { rad2deg(wrap_lon_rad(center_lon)), rad2deg(center_lat) };
  }

  const bool newton = solver == InverseSolver::Newton;
  const VertexGcParams& vgc = P.topo->vgc;
  std::pair<double,double> polar;
  switch (active_projection()) {
    case FaceProjection::Fuller:
      polar = newton ? fuller_face_polar_newton(x, y) : fuller_face_polar(x, y);
      break;
    case FaceProjection::IVEA:
      polar = newton ? ivea_face_polar_newton(vgc, x, y) : ivea_face_polar(vgc, x, y);
      break;
    case FaceProjection::ISEA:
      polar = newton ? snyder_face_polar_newton(sp, x, y) : snyder_face_polar(sp, x, y);
      break;
  }
  const auto [z, face_az] = polar;

  // Add the per-face azimuth bias (radians)
  double azimuth = face_az + snyder_get_face_azimuth_offset(face);
  while (azimuth <= -kPi) azimuth += kTwoPi;
  while (azimuth >   kPi) azimuth -= kTwoPi;

  // Great circle from the face centre, in the centre's east/north frame taken
  // at its longitude, which stays defined when the centre lies at a pole.
  const double north = std::sin(z) * std::cos(azimuth);
  const double east  = std::sin(z) * std::sin(azimuth);
  const double up    = std::cos(z);
  const double radial = up * center_coslat - north * center_sinlat;
  const double sl = std::sin(center_lon), cl = std::cos(center_lon);
  const double px = radial * cl - east * sl;
  const double py = radial * sl + east * cl;
  const double pz = up * center_sinlat + north * center_coslat;
  const double rho = std::hypot(px, py);
  const double lat = std::atan2(pz, rho);
  // A pole has no longitude; it keeps the centre's.
  const double lon = rho < 1e-12 ? wrap_lon_rad(center_lon) : std::atan2(py, px);

  return { rad2deg(lon), rad2deg(lat) };
}

} // namespace hexify
