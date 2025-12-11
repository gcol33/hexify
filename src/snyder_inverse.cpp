#include "snyder_inverse.h"
#include "snyder_forward.h"
#include "core_icosa.h"
#include "constants.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <stdexcept>

namespace {

using hexify::kPi;
using hexify::kTwoPi;
using hexify::kDegToRad;
using hexify::k2PiOver3;
using hexify::k4PiOver3;

// Snyder projection constants
constexpr double R1    = 0.9103832815;
constexpr double R1S   = R1 * R1;
constexpr double DH    = 37.37736814 * kDegToRad;
constexpr double GH    = 36.0 * kDegToRad;
constexpr double cot30 = 1.0 / std::tan(30.0 * kDegToRad);
constexpr double tanDH = std::tan(DH);

// Snyder face-plane normalization (same as forward)
constexpr double originXOff = 0.6022955029;
constexpr double originYOff = 0.3477354707;
constexpr double icosaEdge  = 2.0 * originXOff;

struct PrecCfg { double tol; int max_iters; };
const PrecCfg MODE_FAST    { 1e-10,  25 };
const PrecCfg MODE_DEFAULT { 1e-12,  40 };
const PrecCfg MODE_HIGH    { 1e-14,  80 };
const PrecCfg MODE_ULTRA   { 1e-15, 120 };

PrecCfg CFG = MODE_DEFAULT;
bool    VERBOSE = false;

int ST_calls = 0, ST_iters_total = 0, ST_iters_max = 0, ST_capped = 0;
constexpr double EPS_BRANCH = 1e-15;

inline double wrap_lon_rad(double L) {
  double t = std::fmod(L + kPi, kTwoPi);
  if (t < 0) t += kTwoPi;
  return t - kPi;
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
  const double cent_lon = C[face].lon;
  const double cent_lat = C[face].lat;
  const double cent_sin = std::sin(cent_lat);
  const double cent_cos = std::cos(cent_lat);

  // Per-face azimuth bias (radians)
  const double ddazh = snyder_face_dazh_rad(face);

  // Convert Snyder face-plane (tx,ty) back to px,py (same basis as forward)
  const double px = x * icosaEdge - originXOff;
  const double py = y * icosaEdge - originYOff;

  // Exact face center shortcut
  if (std::abs(px) < EPS_BRANCH && std::abs(py) < EPS_BRANCH) {
    return { rad2deg(wrap_lon_rad(cent_lon)), rad2deg(cent_lat) };
  }

  const double ph   = std::hypot(px, py);

  // Snyder quirk: azimuth uses atan2(x, y) (not atan2(y, x))
  double azh1 = std::atan2(px, py);
  if (azh1 < 0.0) azh1 += kTwoPi;
  const double azh0 = azh1;

  // Reduce to [0,120°) sector for iteration, then restore later
  if (azh1 > k2PiOver3 && azh1 <= k4PiOver3) azh1 -= k2PiOver3;
  if (azh1 > k4PiOver3)                      azh1 -= k4PiOver3;

  double azh = azh1;

  // Newton iteration to solve Snyder auxiliary angle
  int iters = 0;
  if (std::abs(azh1) > EPS_BRANCH) {
    const double agh = (R1S * tanDH * tanDH) / (2.0 * (1.0 / std::tan(azh1) + cot30));
    while (true) {
      double h_arg = std::sin(azh) * std::sin(GH) * std::cos(DH)
                   - std::cos(azh) * std::cos(GH);
      h_arg = clampd(h_arg, -1.0, 1.0); // from core_icosa.h
      const double h = std::acos(h_arg);

      const double fazh  = agh - azh - GH + (kPi - h);

      double denom = std::sin(h);
      if (std::abs(denom) < 1e-18) denom = 1e-18;

      const double flazh = (( std::cos(azh) * std::sin(GH) * std::cos(DH)
                            + std::sin(azh) * std::cos(GH)) / denom) - 1.0;

      const double dazh = -fazh / flazh;
      azh += dazh;
      ++iters;

      if (std::abs(dazh) <= cfg.tol) break;
      if (iters >= cfg.max_iters) { ++ST_capped; break; }
    }
  } else {
    azh = 0.0; // radial line through the face center
  }

  ++ST_calls;
  ST_iters_total += iters;
  if (iters > ST_iters_max) ST_iters_max = iters;

  // Recover z from radial distance
  const double dz    = std::atan2(tanDH, std::cos(azh) + cot30 * std::sin(azh));
  double denom_azh1  = (std::cos(azh1) + cot30 * std::sin(azh1));
  if (std::abs(denom_azh1) < 1e-18) denom_azh1 = 1e-18;

  double sdz2 = std::sin(dz / 2.0);
  if (std::abs(sdz2) < 1e-18) sdz2 = 1e-18;

  const double fh  = tanDH / (2.0 * denom_azh1 * sdz2);
  double arg = (ph / (2.0 * R1 * fh));
  arg = clampd(arg, -1.0, 1.0);
  const double z   = 2.0 * std::asin(arg);

  // Restore original 120° sector and add per-face azimuth
  if (azh0 >= k2PiOver3 && azh0 < k4PiOver3) azh += k2PiOver3;
  if (azh0 >= k4PiOver3)                     azh += k4PiOver3;

  azh += ddazh;
  while (azh <= -kPi) azh += kTwoPi;
  while (azh >   kPi) azh -= kTwoPi;

  // Great-circle from face center (cent_lon, cent_lat)
  double sinlat = cent_sin * std::cos(z) + cent_cos * std::sin(z) * std::cos(azh);
  sinlat = clampd(sinlat, -1.0, 1.0);
  const double lat = std::asin(sinlat);

  double lon;
  if (std::abs(std::abs(lat) - (kPi/2.0)) < 1e-12) {
    lon = cent_lon; // poles: azimuth undefined, keep center longitude
  } else {
    double sinlon = std::sin(azh) * std::sin(z) / std::cos(lat);
    double coslon = (std::cos(z) - cent_sin * std::sin(lat)) / (cent_cos * std::cos(lat));
    sinlon = clampd(sinlon, -1.0, 1.0);
    coslon = clampd(coslon, -1.0, 1.0);
    lon = wrap_lon_rad(cent_lon + std::atan2(sinlon, coslon));
  }

  return { rad2deg(lon), rad2deg(lat) };
}

} // namespace hexify
