#include "projection_inverse.h"
#include "projection_forward.h"
#include "projection_fuller.h"
#include "polyhedron.h"
#include "constants.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace {

using hexify::kPi;
using hexify::kTwoPi;
using hexify::k2PiOver3;
using hexify::kEpsBranch;
using hexify::safe_denom;
using hexify::SnyderParams;

struct PrecCfg { double tol; int max_iters; };
const PrecCfg MODE_FAST    { 1e-10,  25 };
const PrecCfg MODE_DEFAULT { 1e-12,  40 };
const PrecCfg MODE_HIGH    { 1e-14,  80 };
const PrecCfg MODE_ULTRA   { 1e-15, 120 };

PrecCfg CFG = MODE_DEFAULT;
bool    VERBOSE = false;

int ST_calls = 0, ST_iters_total = 0, ST_iters_max = 0, ST_capped = 0;

inline double wrap_lon_rad(double L) {
  double t = std::fmod(L + kPi, kTwoPi);
  if (t < 0) t += kTwoPi;
  return t - kPi;
}

// =============================================================================
// Newton-Raphson solver for Snyder auxiliary angle
// =============================================================================

struct NewtonResult {
  double azimuth;    // converged azimuth angle
  int    iterations; // number of iterations performed
  bool   converged;  // true if converged within tolerance
};

/**
 * Computes f(azimuth) and f'(azimuth) for Newton-Raphson iteration.
 *
 * The residual function is: f(azimuth) = agh - azimuth - G + (π - h)
 * where h = acos(sin(azimuth)*sin(G)*cos(EL) - cos(azimuth)*cos(G))
 *
 * @param azimuth Current azimuth estimate
 * @param agh Pre-computed auxiliary constant
 * @return pair<residual, derivative>
 */
inline std::pair<double, double> newton_residual_and_derivative(const SnyderParams& sp,
                                                                double azimuth, double agh) {
  const double sin_azimuth = std::sin(azimuth);
  const double cos_azimuth = std::cos(azimuth);

  // Compute h = acos(sin(azimuth)*sin(G)*cos(EL) - cos(azimuth)*cos(G))
  double h_arg = sin_azimuth * sp.sin_g * sp.cos_el - cos_azimuth * sp.cos_g;
  h_arg = hexify::clampd(h_arg, -1.0, 1.0);
  const double h = std::acos(h_arg);

  // Residual: f(azimuth) = agh - azimuth - G + (π - h)
  const double residual = agh - azimuth - sp.g_angle + (kPi - h);

  // Derivative: f'(azimuth) = (cos(azimuth)*sin(G)*cos(EL) + sin(azimuth)*cos(G)) / sin(h) - 1
  const double sin_h = safe_denom(std::sin(h));

  const double derivative = ((cos_azimuth * sp.sin_g * sp.cos_el + sin_azimuth * sp.cos_g) / sin_h) - 1.0;

  return {residual, derivative};
}

/**
 * Solves for the Snyder auxiliary azimuth angle using Newton-Raphson iteration.
 *
 * @param azimuth_initial Initial azimuth estimate (reduced to [0, 120°) sector)
 * @param cfg Precision configuration (tolerance and max iterations)
 * @return NewtonResult with converged angle, iteration count, and convergence status
 */
NewtonResult solve_snyder_azimuth(const SnyderParams& sp, double azimuth_initial,
                                  const PrecCfg& cfg) {
  // Special case: azimuth near zero (radial line through face center)
  if (std::abs(azimuth_initial) <= kEpsBranch) {
    return {0.0, 0, true};
  }

  // Pre-compute the auxiliary constant agh
  const double agh = (sp.r1_squared * sp.tan_el * sp.tan_el) / (2.0 * (1.0 / std::tan(azimuth_initial) + sp.cot_30));

  double azimuth = azimuth_initial;
  for (int iter = 0; iter < cfg.max_iters; ++iter) {
    auto [residual, derivative] = newton_residual_and_derivative(sp, azimuth, agh);

    const double delta = -residual / derivative;
    azimuth += delta;

    if (std::abs(delta) <= cfg.tol) {
      return {azimuth, iter + 1, true};
    }
  }

  // Did not converge within max iterations
  return {azimuth, cfg.max_iters, false};
}

} // anon

namespace hexify {

void snyder_inv_set_precision(const std::string& mode,
                              double tol_override,
                              int    max_iters_override) {
  if (!mode.empty()) {
    if      (mode == "fast")    CFG = MODE_FAST;
    else if (mode == "default") CFG = MODE_DEFAULT;
    else if (mode == "high")    CFG = MODE_HIGH;
    else if (mode == "ultra")   CFG = MODE_ULTRA;
    else throw std::runtime_error("Unknown precision mode: " + mode);
  }
  if (tol_override       >= 0.0) CFG.tol       = tol_override;
  if (max_iters_override >= 0  ) CFG.max_iters = max_iters_override;
}

std::pair<double,double> snyder_inv_get_precision() {
  return {CFG.tol, static_cast<double>(CFG.max_iters)};
}

void snyder_inv_set_verbose(bool v) { VERBOSE = v; }

std::tuple<int,int,int,int> snyder_inv_get_stats_and_reset() {
  auto out = std::make_tuple(ST_calls, ST_iters_total, ST_iters_max, ST_capped);
  ST_calls = ST_iters_total = ST_iters_max = ST_capped = 0;
  return out;
}

} // namespace hexify

namespace {

// Snyder's inverse on a face: face-plane (x, y) -> distance z from the face
// centre and azimuth from the face's first vertex.
std::pair<double,double> snyder_face_polar(const SnyderParams& sp, double x, double y,
                                           const PrecCfg& cfg) {
  const double px = x * sp.edge - sp.origin_x_off;
  const double py = y * sp.edge - sp.origin_y_off;

  // Radial distance in face plane (Snyder notation: ρ)
  const double rho   = std::hypot(px, py);

  // Snyder quirk: azimuth uses atan2(x, y) (not atan2(y, x))
  double azimuth_transformed = std::atan2(px, py);
  if (azimuth_transformed < 0.0) azimuth_transformed += kTwoPi;
  if (azimuth_transformed >= kTwoPi) azimuth_transformed -= kTwoPi;

  // Reduce to its 120° sector for iteration, then restore later
  const int sector = hexify::azimuth_sector(azimuth_transformed);
  azimuth_transformed -= sector * k2PiOver3;

  // Solve for azimuth using Newton-Raphson iteration
  NewtonResult newton = solve_snyder_azimuth(sp, azimuth_transformed, cfg);
  double azimuth = newton.azimuth;

  // Update statistics
  ++ST_calls;
  ST_iters_total += newton.iterations;
  if (newton.iterations > ST_iters_max) ST_iters_max = newton.iterations;
  if (!newton.converged) ++ST_capped;

  // Recover z (great-circle distance from face center) from radial distance
  // Snyder's auxiliary angle for the sector (Snyder notation: δ_z)
  const double dz_angle = std::atan2(sp.tan_el, std::cos(azimuth) + sp.cot_30 * std::sin(azimuth));
  const double denom = safe_denom(std::cos(azimuth_transformed) + sp.cot_30 * std::sin(azimuth_transformed));
  const double sin_half_dz = safe_denom(std::sin(dz_angle / 2.0));

  // Snyder's 'f' scale factor (Snyder notation: f)
  const double f_scale = sp.tan_el / (2.0 * denom * sin_half_dz);
  double arg = (rho / (2.0 * sp.r1 * f_scale));
  arg = hexify::clampd(arg, -1.0, 1.0);
  // Great-circle distance z from face center (Snyder notation: z)
  const double z = 2.0 * std::asin(arg);

  // Restore original 120° sector
  return {z, azimuth + sector * k2PiOver3};
}

// Fuller's inverse on a face, with the Snyder solver's statistics.
std::pair<double,double> fuller_polar(double x, double y, const PrecCfg& cfg) {
  int iters = 0;
  const auto za = hexify::fuller_face_polar(x, y, cfg.tol, cfg.max_iters, &iters);
  ++ST_calls;
  ST_iters_total += iters;
  if (iters > ST_iters_max) ST_iters_max = iters;
  if (iters >= cfg.max_iters) ++ST_capped;
  return za;
}

} // anon

namespace hexify {

std::pair<double,double> face_xy_to_ll(double x, double y, int face,
                                       double tol_override,
                                       int    max_iters_override)
{
  const PolyData& P = poly();
  if (face < 0 || face >= P.n_faces()) throw std::runtime_error("face out of range for the solid");
  const SnyderParams& sp = P.topo->snyder;

  // per-call precision
  PrecCfg cfg = CFG;
  if (tol_override       >= 0.0) cfg.tol       = tol_override;
  if (max_iters_override >= 0  ) cfg.max_iters = max_iters_override;

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

  const auto [z, face_az] = active_projection() == FaceProjection::Fuller
    ? fuller_polar(x, y, cfg) : snyder_face_polar(sp, x, y, cfg);

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
