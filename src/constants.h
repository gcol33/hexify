// constants.h - Shared mathematical constants for hexify
//
// All constants computed to full IEEE 754 double precision (15-17 significant digits).
// Values verified against Wolfram Alpha / mpfr where applicable.
//
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_CONSTANTS_H
#define HEXIFY_CONSTANTS_H

#include <cmath>

namespace hexify {

// =============================================================================
// Fundamental Mathematical Constants
// =============================================================================

// Pi and multiples (full double precision)
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 6.283185307179586476925286766559005768;
constexpr double kPiOver2 = 1.570796326794896619231321691639751442;
constexpr double kPiOver3 = 1.047197551196597746154214461093167628;
constexpr double kPiOver6 = 0.523598775598298873077107230546583814;

// =============================================================================
// Square Roots
// =============================================================================

constexpr double kSqrt3 = 1.732050807568877293527446341505872367;
constexpr double kSqrt7 = 2.645751311064590590501615753639260426;
constexpr double kSqrt21 = 4.582575694955840006588047193728008489;  // sqrt(3 * 7)

// =============================================================================
// Trigonometric Values
// =============================================================================

constexpr double kSin60 = 0.866025403784438646763723170752936183;  // sqrt(3)/2
constexpr double kCos60 = 0.5;
constexpr double kSin30 = 0.5;
constexpr double kCos30 = 0.866025403784438646763723170752936183;  // sqrt(3)/2

// =============================================================================
// Degree/Radian Conversion
// =============================================================================

constexpr double kDegToRad = 0.017453292519943295769236907684886127;  // pi/180
constexpr double kRadToDeg = 57.29577951308232087679815481410517033;  // 180/pi

// =============================================================================
// ISEA Projection Constants
// =============================================================================

// Aperture 7 rotation angle: arctan(sqrt(3)/5) in degrees
// Exact: atan(sqrt(3)/5) = 19.10660535086909...°
// Cross-checked against DGGRID's M_AP7_ROT_DEGS (src/lib/dglib/include/dglib/DgConstants.h)
constexpr double kAp7RotDeg = 19.106605350869094394517474740130082234976075229;

// =============================================================================
// Snyder Projection Sector Angles
// =============================================================================

// Triangle sector boundaries (radians) - used for azimuth reduction
constexpr double k2PiOver3 = 2.094395102393195492308428922186335256;   // 120° = 2π/3
constexpr double k4PiOver3 = 4.188790204786390984616857844372670512;   // 240° = 4π/3

// The 120° sector of a face holding an azimuth in [0, 2π): 0, 1 or 2. Each
// sector is half-open, [k 120°, (k + 1) 120°), so an azimuth on a sector
// boundary belongs to exactly one sector, the same in both directions of the
// projection. Reduce an azimuth with `az - sector * k2PiOver3` and restore it
// with `+ sector * k2PiOver3`.
inline int azimuth_sector(double az) {
  if (az >= k4PiOver3) return 2;
  if (az >= k2PiOver3) return 1;
  return 0;
}

// =============================================================================
// Fuller Projection Constants (Gray 1995, Crider 2008)
// =============================================================================
// References: Gray, R.W. (1995). "Exact Transformation Equations for Fuller's
// World Map". Cartographica 32(3): 17-25. Crider, J.E. (2008). "Exact
// Equations for Fuller's Map Projection and Inverse". Cartographica 43(1):
// 67-72.
//
// Fuller unfolds each spherical face onto a plane triangle whose edges keep
// their arc length, so the plane triangle's edge is the icosahedron's edge arc.

// ARC: arc length of an icosahedron edge on the unit sphere, atan(2).
inline const double kFullerArc = std::atan(2.0);

// alpha: half the edge arc.
inline const double kFullerAlpha = 0.5 * kFullerArc;

// EL: chord length of an icosahedron edge, sqrt(8) / sqrt(5 + sqrt(5)).
inline const double kFullerEL = std::sqrt(8.0) / std::sqrt(5.0 + std::sqrt(5.0));

// DVE: distance from the centre to an edge midpoint, sqrt(3 + sqrt(5)) /
// sqrt(5 + sqrt(5)).
inline const double kFullerDVE = std::sqrt(3.0 + std::sqrt(5.0)) / std::sqrt(5.0 + std::sqrt(5.0));

// Height of the face plane above the centre: sqrt(5 + 2 sqrt(5)) / sqrt(15).
inline const double kFullerZ0 = std::sqrt(5.0 + 2.0 * std::sqrt(5.0)) / std::sqrt(15.0);

// tan(alpha) = EL / (2 DVE).
inline const double kFullerTanAlpha = kFullerEL / (2.0 * kFullerDVE);

// Latitude of the ten icosahedron vertices off the poles in the standard
// orientation: atan(1/2).
inline const double kIcosaVertexLatDeg = std::atan(0.5) * kRadToDeg;

// =============================================================================
// Grid Bounds
// =============================================================================
// Mirrors R/constants.R's MIN_RESOLUTION/MAX_RESOLUTION. Resolutions outside
// this range are rejected before they can reach shift-overflow or
// scale-to-infinity arithmetic in the grid-dimension calculations.

constexpr int kMinResolution = 0;
constexpr int kMaxResolution = 30;

// =============================================================================
// Numerical Precision Constants
// =============================================================================

// Minimum denominator value to prevent division by zero in floating-point math
constexpr double kMinDenom = 1e-18;

// Epsilon for branching decisions (very small values treated as zero)
constexpr double kEpsBranch = 1e-15;

// =============================================================================
// Inline Utility Functions
// =============================================================================

/**
 * Returns a non-zero denominator suitable for division.
 * If abs(x) < kMinDenom, returns kMinDenom with the original sign (or positive if x=0).
 * This prevents division-by-zero while preserving the sign of near-zero values.
 */
inline double safe_denom(double x) noexcept {
  if (x >= 0.0) {
    return (x < kMinDenom) ? kMinDenom : x;
  } else {
    return (x > -kMinDenom) ? -kMinDenom : x;
  }
}

} // namespace hexify

#endif // HEXIFY_CONSTANTS_H
