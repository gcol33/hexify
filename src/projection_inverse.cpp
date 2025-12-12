#include "projection_inverse.h"
#include "projection_forward.h"
#include "icosahedron.h"
#include "constants.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace {

using hexify::kPi;
using hexify::kTwoPi;
using hexify::kPiOver6;
using hexify::k2PiOver3;
using hexify::k4PiOver3;
using hexify::kSnyderR1;
using hexify::kSnyderR1Squared;
using hexify::kSnyderElAngle;
using hexify::kSnyderGAngle;
using hexify::kSnyderOriginXOff;
using hexify::kSnyderOriginYOff;
using hexify::kSnyderIcosaEdge;
using hexify::kEpsBranch;
using hexify::safe_denom;

// Derived trigonometric values (not constexpr because std::tan isn't constexpr)
static const double COT_30 = 1.0 / std::tan(kPiOver6);  // cot(30°) = 1/tan(π/6)
static const double TAN_EL = std::tan(kSnyderElAngle);

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

// Derived trigonometric values for Newton iteration
static const double SIN_G = std::sin(kSnyderGAngle);
static const double COS_G = std::cos(kSnyderGAngle);
static const double COS_EL = std::cos(kSnyderElAngle);

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
inline std::pair<double, double> newton_residual_and_derivative(double azimuth, double agh) {
  const double sin_azimuth = std::sin(azimuth);
  const double cos_azimuth = std::cos(azimuth);

  // Compute h = acos(sin(azimuth)*sin(G)*cos(EL) - cos(azimuth)*cos(G))
  double h_arg = sin_azimuth * SIN_G * COS_EL - cos_azimuth * COS_G;
  h_arg = hexify::clampd(h_arg, -1.0, 1.0);
  const double h = std::acos(h_arg);

  // Residual: f(azimuth) = agh - azimuth - G + (π - h)
  const double residual = agh - azimuth - kSnyderGAngle + (kPi - h);

  // Derivative: f'(azimuth) = (cos(azimuth)*sin(G)*cos(EL) + sin(azimuth)*cos(G)) / sin(h) - 1
  const double sin_h = safe_denom(std::sin(h));

  const double derivative = ((cos_azimuth * SIN_G * COS_EL + sin_azimuth * COS_G) / sin_h) - 1.0;

  return {residual, derivative};
}

/**
 * Solves for the Snyder auxiliary azimuth angle using Newton-Raphson iteration.
 *
 * @param azimuth_initial Initial azimuth estimate (reduced to [0, 120°) sector)
 * @param cfg Precision configuration (tolerance and max iterations)
 * @return NewtonResult with converged angle, iteration count, and convergence status
 */
NewtonResult solve_snyder_azimuth(double azimuth_initial, const PrecCfg& cfg) {
  // Special case: azimuth near zero (radial line through face center)
  if (std::abs(azimuth_initial) <= kEpsBranch) {
    return {0.0, 0, true};
  }

  // Pre-compute the auxiliary constant agh
  const double agh = (kSnyderR1Squared * TAN_EL * TAN_EL) / (2.0 * (1.0 / std::tan(azimuth_initial) + COT_30));

  double azimuth = azimuth_initial;
  for (int iter = 0; iter < cfg.max_iters; ++iter) {
    auto [residual, derivative] = newton_residual_and_derivative(azimuth, agh);

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

std::pair<double,double> face_xy_to_ll(double x, double y, int face,
                                       double tol_override,
                                       int    max_iters_override)
{
  if (face < 0 || face >= 20) throw std::runtime_error("face must be 0..19");

  // per-call precision
  PrecCfg cfg = CFG;
  if (tol_override       >= 0.0) cfg.tol       = tol_override;
  if (max_iters_override >= 0  ) cfg.max_iters = max_iters_override;

  // Face centers are in radians
  const auto& C = face_centers();
  const double center_lon = C[face].lon;
  const double center_lat = C[face].lat;
  const double center_sinlat = std::sin(center_lat);
  const double center_coslat = std::cos(center_lat);

  // Per-face azimuth bias (radians)
  const double face_azimuth = snyder_get_face_azimuth_offset(face);

  // Convert Snyder face-plane (tx,ty) back to px,py (same basis as forward)
  const double px = x * kSnyderIcosaEdge - kSnyderOriginXOff;
  const double py = y * kSnyderIcosaEdge - kSnyderOriginYOff;

  // Exact face center shortcut
  if (std::abs(px) < kEpsBranch && std::abs(py) < kEpsBranch) {
    return { rad2deg(wrap_lon_rad(center_lon)), rad2deg(center_lat) };
  }

  // Radial distance in face plane (Snyder notation: ρ)
  const double rho   = std::hypot(px, py);

  // Snyder quirk: azimuth uses atan2(x, y) (not atan2(y, x))
  double azimuth_transformed = std::atan2(px, py);
  if (azimuth_transformed < 0.0) azimuth_transformed += kTwoPi;
  const double azimuth_original = azimuth_transformed;

  // Reduce to [0,120°) sector for iteration, then restore later
  if (azimuth_transformed > k2PiOver3 && azimuth_transformed <= k4PiOver3) azimuth_transformed -= k2PiOver3;
  if (azimuth_transformed > k4PiOver3) azimuth_transformed -= k4PiOver3;

  // Solve for azimuth using Newton-Raphson iteration
  NewtonResult newton = solve_snyder_azimuth(azimuth_transformed, cfg);
  double azimuth = newton.azimuth;

  // Update statistics
  ++ST_calls;
  ST_iters_total += newton.iterations;
  if (newton.iterations > ST_iters_max) ST_iters_max = newton.iterations;
  if (!newton.converged) ++ST_capped;

  // Recover z (great-circle distance from face center) from radial distance
  // Snyder's auxiliary angle for the sector (Snyder notation: δ_z)
  const double dz_angle = std::atan2(TAN_EL, std::cos(azimuth) + COT_30 * std::sin(azimuth));
  const double denom = safe_denom(std::cos(azimuth_transformed) + COT_30 * std::sin(azimuth_transformed));
  const double sin_half_dz = safe_denom(std::sin(dz_angle / 2.0));

  // Snyder's 'f' scale factor (Snyder notation: f)
  const double f_scale = TAN_EL / (2.0 * denom * sin_half_dz);
  double arg = (rho / (2.0 * kSnyderR1 * f_scale));
  arg = clampd(arg, -1.0, 1.0);
  // Great-circle distance z from face center (Snyder notation: z)
  const double z = 2.0 * std::asin(arg);

  // Restore original 120° sector and add per-face azimuth
  if (azimuth_original >= k2PiOver3 && azimuth_original < k4PiOver3) azimuth += k2PiOver3;
  if (azimuth_original >= k4PiOver3) azimuth += k4PiOver3;

  azimuth += face_azimuth;
  while (azimuth <= -kPi) azimuth += kTwoPi;
  while (azimuth >   kPi) azimuth -= kTwoPi;

  // Great-circle from face center
  double sinlat = center_sinlat * std::cos(z) + center_coslat * std::sin(z) * std::cos(azimuth);
  sinlat = clampd(sinlat, -1.0, 1.0);
  const double lat = std::asin(sinlat);

  double lon;
  if (std::abs(std::abs(lat) - (kPi/2.0)) < 1e-12) {
    lon = center_lon; // poles: azimuth undefined, keep center longitude
  } else {
    double sinlon = std::sin(azimuth) * std::sin(z) / std::cos(lat);
    double coslon = (std::cos(z) - center_sinlat * std::sin(lat)) / (center_coslat * std::cos(lat));
    sinlon = clampd(sinlon, -1.0, 1.0);
    coslon = clampd(coslon, -1.0, 1.0);
    lon = wrap_lon_rad(center_lon + std::atan2(sinlon, coslon));
  }

  return { rad2deg(lon), rad2deg(lat) };
}

} // namespace hexify
