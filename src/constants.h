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

// =============================================================================
// Snyder ISEA Projection Constants (from Snyder 1992)
// =============================================================================
// Reference: Snyder, J.P. (1992). "An Equal-Area Map Projection for Polyhedral Globes"
// Cartographica 29(1): 10-21.

// R1: Scale factor for equal-area property (Snyder notation: R')
constexpr double kSnyderR1 = 0.9103832815;
constexpr double kSnyderR1Squared = kSnyderR1 * kSnyderR1;

// SNYDER_EL_ANGLE: θ - spherical angle from face center to edge midpoint (37.377...°)
// This is the "el" angle in Snyder's notation
constexpr double kSnyderElAngleDeg = 37.37736814;
constexpr double kSnyderElAngle = kSnyderElAngleDeg * kDegToRad;

// SNYDER_G_ANGLE: G - angle at icosahedron vertices (36° exactly)
constexpr double kSnyderGAngleDeg = 36.0;
constexpr double kSnyderGAngle = kSnyderGAngleDeg * kDegToRad;

// Face-plane origin offsets (for normalizing projected coordinates)
constexpr double kSnyderOriginXOff = 0.6022955029;
constexpr double kSnyderOriginYOff = 0.3477354707;
constexpr double kSnyderIcosaEdge = 2.0 * kSnyderOriginXOff;

// =============================================================================
// PLANE Coordinate Layout Table (Icosahedron Unfolding)
// =============================================================================
// Each triangle has a rotation (in 60° increments) and an offset position
// in the unfolded PLANE coordinate system. Standard ISEA icosahedron layout.

// Structure to hold triangle transformation parameters
struct PlaneTriLayout {
    int rot60;      // Rotation in 60° increments (0, 3, etc.)
    double offset_x;  // X offset in PLANE coordinates
    double offset_y;  // Y offset in PLANE coordinates
};

// M_SIN60 = sqrt(3)/2 ≈ 0.866025...
// The layout creates a 5.5 × ~1.73 unit rectangle containing all 20 triangles

constexpr PlaneTriLayout kPlaneLayout[20] = {
    // Upper row (faces 0-4): rot60=0, y = 2*sin60
    {0, 0.0, 2.0 * kSin60},  // face 0
    {0, 1.0, 2.0 * kSin60},  // face 1
    {0, 2.0, 2.0 * kSin60},  // face 2
    {0, 3.0, 2.0 * kSin60},  // face 3
    {0, 4.0, 2.0 * kSin60},  // face 4
    // Upper row continued (faces 5-9): rot60=3, y = 2*sin60
    {3, 1.0, 2.0 * kSin60},  // face 5
    {3, 2.0, 2.0 * kSin60},  // face 6
    {3, 3.0, 2.0 * kSin60},  // face 7
    {3, 4.0, 2.0 * kSin60},  // face 8
    {3, 5.0, 2.0 * kSin60},  // face 9
    // Lower row (faces 10-14): rot60=0, y = sin60
    {0, 0.5, kSin60},        // face 10
    {0, 1.5, kSin60},        // face 11
    {0, 2.5, kSin60},        // face 12
    {0, 3.5, kSin60},        // face 13
    {0, 4.5, kSin60},        // face 14
    // Lower row continued (faces 15-19): rot60=3, y = sin60
    {3, 1.5, kSin60},        // face 15
    {3, 2.5, kSin60},        // face 16
    {3, 3.5, kSin60},        // face 17
    {3, 4.5, kSin60},        // face 18
    {3, 5.5, kSin60}         // face 19
};

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
