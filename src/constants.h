// constants.h - Shared mathematical constants for hexify
//
// All constants computed to full IEEE 754 double precision (15-17 significant digits).
// Values verified against Wolfram Alpha / mpfr where applicable.
//
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_CONSTANTS_H
#define HEXIFY_CONSTANTS_H

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

// Aperture 7 rotation angle: arctan(sqrt(3/7)) in degrees
// Exact: atan(sqrt(3/7)) = 19.10660535003926...°
constexpr double kAp7RotDeg = 19.10660535003926406149339781619697490;

// =============================================================================
// Snyder Projection Sector Angles
// =============================================================================

// Triangle sector boundaries (radians) - used for azimuth reduction
constexpr double k2PiOver3 = 2.094395102393195492308428922186335256;   // 120° = 2π/3
constexpr double k4PiOver3 = 4.188790204786390984616857844372670512;   // 240° = 4π/3

} // namespace hexify

#endif // HEXIFY_CONSTANTS_H
