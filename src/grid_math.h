// grid_math.h - Shared hexagonal grid mathematics
//
// This module provides the fundamental mathematical operations for hexagonal
// discrete global grid systems (DGGS). It consolidates operations that were
// previously duplicated across aperture-specific files.
//
// ============================================================================
// HEXAGON ROTATION CLASSES
// ============================================================================
//
// Hexagonal grids come in different orientations called "Rotation Classes"
// (terminology from Sahr et al. 2003):
//
//   ROTATION CLASS I (0-degree, flat-top)
//   --------------------------------------
//   - Standard orientation (0 degrees)
//   - Used by aperture-4 always
//   - Flat edge on top
//
//   ROTATION CLASS II (30-degree, pointy-top)
//   ------------------------------------------
//   - Rotated 30 degrees from Rotation Class I
//   - Alternates with Rotation Class I in aperture-3
//   - Pointed vertex on top
//
//   ROTATION CLASS III (Aperture-7 only)
//   -------------------------------------
//   Rotated by arctan(sqrt(3/7)) = ~19.1 degrees from the base orientation.
//   Two variants alternate by resolution:
//     - Rotation Class III-A (even res): Rotation Class I + 19.1 deg = ~19 deg
//     - Rotation Class III-B (odd res):  Rotation Class II + 19.1 deg = ~49 deg
//
// ============================================================================
// COORDINATE SYSTEMS
// ============================================================================
//
// CUBE COORDINATES (q, r, s)
// --------------------------
// Three-axis system where q + r + s = 0. Provides elegant nearest-neighbor
// rounding via the "round and fix" algorithm. Used internally for quantization.
//
//             +s
//              |
//              |
//       +q ----+---- -q
//              |
//              |
//             -s
//
// OFFSET COORDINATES (i, j)
// -------------------------
// Two-axis system output by quantization. Maps directly to cell indices.
// In our implementation: i = q, j = r (from cube coordinates).
//
// ============================================================================
// SURROGATE-SUBSTRATE PATTERN
// ============================================================================
//
// Non-Rotation-Class-I grids use a "surrogate-substrate" quantization pattern:
//
//   1. ROTATE input point to align with a "surrogate" Rotation Class I grid
//   2. QUANTIZE in the surrogate grid (using Rotation Class I math)
//   3. GET CENTER of the surrogate cell
//   4. ROTATE BACK to original orientation
//   5. SCALE UP to a finer "substrate" grid and re-quantize
//
// This pattern produces coordinates compatible with hierarchical ISEA grids.
//
// Scale factors by rotation class:
//   - Rotation Class II:     sqrt(3) = ~1.732  (from Rotation Class I surrogate)
//   - Rotation Class III-A:  sqrt(7) = ~2.646  (aperture-7, even resolutions)
//   - Rotation Class III-B:  sqrt(21) = ~4.583 (aperture-7, odd resolutions)
//
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_GRID_MATH_H
#define HEXIFY_GRID_MATH_H

#include "cube_coordinates.h"
#include "constants.h"
#include <cmath>

namespace hexify {

// ============================================================================
// Rotation Utilities
// ============================================================================

/**
 * Rotate a 2D point by an angle (in radians).
 *
 * @param x      Input/output X coordinate
 * @param y      Input/output Y coordinate
 * @param angle  Rotation angle in radians (positive = counter-clockwise)
 */
inline void rotate_point(double& x, double& y, double angle_rad) {
    double c = std::cos(angle_rad);
    double s = std::sin(angle_rad);
    double new_x = x * c - y * s;
    double new_y = x * s + y * c;
    x = new_x;
    y = new_y;
}

/**
 * Rotate a 2D point using pre-computed sin/cos values.
 * More efficient when the same rotation is applied many times.
 *
 * @param x       Input/output X coordinate
 * @param y       Input/output Y coordinate
 * @param cos_a   Cosine of the rotation angle
 * @param sin_a   Sine of the rotation angle
 */
inline void rotate_point_precomputed(double& x, double& y,
                                      double cos_a, double sin_a) {
    double new_x = x * cos_a - y * sin_a;
    double new_y = x * sin_a + y * cos_a;
    x = new_x;
    y = new_y;
}

/**
 * Rotate a 2D point by the inverse angle (negate sin).
 * Used for "rotate back" operations in surrogate-substrate pattern.
 */
inline void rotate_point_inverse(double& x, double& y,
                                  double cos_a, double sin_a) {
    double new_x = x * cos_a + y * sin_a;
    double new_y = -x * sin_a + y * cos_a;
    x = new_x;
    y = new_y;
}

// ============================================================================
// Rotation Class I (0-Degree, Flat-Top) Hexagon Quantization
// ============================================================================
//
// This is the fundamental quantization algorithm. All other rotation classes
// use this as a building block via the surrogate-substrate pattern.

/**
 * Quantize a point to the nearest Rotation Class I (0-degree, flat-top) hexagon.
 *
 * Algorithm:
 *   1. Convert Cartesian (x, y) to cube coordinates (q, r, s)
 *   2. Round each cube coordinate to nearest integer
 *   3. Fix rounding to maintain q + r + s = 0 constraint
 *   4. Output offset coordinates (i, j) = (q, r)
 *
 * @param x       X coordinate in hex grid space
 * @param y       Y coordinate in hex grid space
 * @param out_i   Output: column index (q from cube coords)
 * @param out_j   Output: row index (r from cube coords)
 */
inline void quantize_rotation_classI(double x, double y,
                                     long long& out_i, long long& out_j) {
    // Convert to cube coordinates using flat-top layout
    CubeCoord cube = cartesian_to_cube(x, y, kSqrt3);

    // Round to nearest hex center (maintains q + r + s = 0)
    cube.round_to_nearest();

    // Extract offset coordinates
    out_i = static_cast<long long>(cube.q);
    out_j = static_cast<long long>(cube.r);
}

/**
 * Get the Cartesian center of a Rotation Class I (0-degree, flat-top) hexagon.
 *
 * @param i       Column index
 * @param j       Row index
 * @param out_x   Output: X coordinate of cell center
 * @param out_y   Output: Y coordinate of cell center
 */
inline void center_rotation_classI(long long i, long long j,
                                   double& out_x, double& out_y) {
    cube_to_cartesian(static_cast<double>(i), static_cast<double>(j),
                      out_x, out_y, kSin60);
}

// ============================================================================
// Rotation Class II (30-Degree, Pointy-Top) Hexagon Quantization
// ============================================================================
//
// Rotation Class II hexagons are rotated 30 degrees from Rotation Class I.
// Uses the surrogate-substrate pattern.

/**
 * Quantize a point to the nearest Rotation Class II (30-degree, pointy-top) hexagon.
 *
 * Uses surrogate-substrate pattern:
 *   1. Rotate by -30 deg to align with Rotation Class I surrogate
 *   2. Quantize in Rotation Class I grid
 *   3. Get surrogate center, rotate back by +30 deg
 *   4. Scale by sqrt(3) to substrate, re-quantize in Rotation Class I grid
 *
 * @param x       X coordinate in hex grid space
 * @param y       Y coordinate in hex grid space
 * @param out_i   Output: column index in substrate coordinates
 * @param out_j   Output: row index in substrate coordinates
 */
inline void quantize_rotation_classII(double x, double y,
                                      long long& out_i, long long& out_j) {
    // Pre-computed rotation constants for -30 degrees
    constexpr double cos_neg30 = kCos30;   // cos(-30) = cos(30)
    constexpr double sin_neg30 = -kSin30;  // sin(-30) = -sin(30) = -0.5

    // Step 1: Rotate to Rotation Class I surrogate frame (-30 degrees)
    double sur_x = x * cos_neg30 - y * sin_neg30;
    double sur_y = x * sin_neg30 + y * cos_neg30;

    // Step 2: Quantize in Rotation Class I surrogate
    long long sur_i, sur_j;
    quantize_rotation_classI(sur_x, sur_y, sur_i, sur_j);

    // Step 3: Get surrogate center
    double cen_x, cen_y;
    center_rotation_classI(sur_i, sur_j, cen_x, cen_y);

    // Step 4: Rotate center back to original frame (+30 degrees)
    // Inverse rotation: cos same, sin negated
    double back_x = cen_x * cos_neg30 + cen_y * sin_neg30;
    double back_y = -cen_x * sin_neg30 + cen_y * cos_neg30;

    // Step 5: Scale to substrate (sqrt(3) finer) and re-quantize
    quantize_rotation_classI(back_x * kSqrt3, back_y * kSqrt3, out_i, out_j);
}

// ============================================================================
// Rotation Class III (Aperture-7) Hexagon Quantization
// ============================================================================
//
// Aperture-7 uses hexagons rotated by arctan(sqrt(3/7)) = ~19.1 degrees.
// Two variants alternate by resolution:
//   - Rotation Class III-A (~19 deg, even res): Rotation Class I + 19.1 deg
//   - Rotation Class III-B (~49 deg, odd res):  Rotation Class II + 19.1 deg

// Pre-computed constants for ~19.1 degree rotation
namespace detail {
    // arctan(sqrt(3/7)) in radians
    constexpr double kAp7RotRad = kAp7RotDeg * kDegToRad;
    constexpr double kCos19 = 0.9449111825230680440492389263705078;  // cos(19.106...°)
    constexpr double kSin19 = 0.3273268353539885718950317563490135;  // sin(19.106...°)
}

/**
 * Quantize to Rotation Class III-A (~19-degree) hexagon (aperture-7, even resolutions).
 *
 * Surrogate: Rotation Class I grid rotated by -19.1 degrees
 * Substrate: sqrt(7) times finer than surrogate
 */
inline void quantize_rotation_classIII_A(double x, double y,
                                         long long& out_i, long long& out_j) {
    using namespace detail;

    // Step 1: Rotate to Rotation Class I surrogate frame (-19.1 degrees)
    double sur_x = x * kCos19 + y * kSin19;   // cos(-a) = cos(a)
    double sur_y = -x * kSin19 + y * kCos19;  // sin(-a) = -sin(a)

    // Step 2: Quantize in Rotation Class I surrogate
    long long sur_i, sur_j;
    quantize_rotation_classI(sur_x, sur_y, sur_i, sur_j);

    // Step 3: Get surrogate center
    double cen_x, cen_y;
    center_rotation_classI(sur_i, sur_j, cen_x, cen_y);

    // Step 4: Rotate back to original frame (+19.1 degrees)
    double back_x = cen_x * kCos19 - cen_y * kSin19;
    double back_y = cen_x * kSin19 + cen_y * kCos19;

    // Step 5: Scale to substrate (sqrt(7) finer) and re-quantize
    quantize_rotation_classI(back_x * kSqrt7, back_y * kSqrt7, out_i, out_j);
}

/**
 * Quantize to Rotation Class III-B (~49-degree) hexagon (aperture-7, odd resolutions).
 *
 * Surrogate: Rotation Class II grid rotated by -19.1 degrees
 *            (equivalent to Rotation Class I grid rotated by -49.1 degrees)
 * Substrate: sqrt(21) times finer than surrogate
 */
inline void quantize_rotation_classIII_B(double x, double y,
                                         long long& out_i, long long& out_j) {
    using namespace detail;

    // Step 1: Rotate to surrogate frame (-19.1 degrees)
    double sur_x = x * kCos19 + y * kSin19;
    double sur_y = -x * kSin19 + y * kCos19;

    // Step 2: The surrogate is Rotation Class II (rotated 30 deg from Rotation Class I).
    // To quantize, first rotate to Rotation Class I orientation.
    constexpr double cos_neg30 = kCos30;
    constexpr double sin_neg30 = -0.5;

    double rotated_x = sur_x * cos_neg30 - sur_y * sin_neg30;
    double rotated_y = sur_x * sin_neg30 + sur_y * cos_neg30;

    // Step 3: Quantize in Rotation Class I grid
    long long sur_i, sur_j;
    quantize_rotation_classI(rotated_x, rotated_y, sur_i, sur_j);

    // Step 4: Get Rotation Class I center
    double cen_x, cen_y;
    center_rotation_classI(sur_i, sur_j, cen_x, cen_y);

    // Step 5: Rotate back to Rotation Class II orientation (+30 degrees)
    double rotated_back_x = cen_x * cos_neg30 + cen_y * sin_neg30;
    double rotated_back_y = -cen_x * sin_neg30 + cen_y * cos_neg30;

    // Step 6: Rotate back to original frame (+19.1 degrees)
    double back_x = rotated_back_x * kCos19 - rotated_back_y * kSin19;
    double back_y = rotated_back_x * kSin19 + rotated_back_y * kCos19;

    // Step 7: Scale to substrate (sqrt(21) finer) and re-quantize
    quantize_rotation_classI(back_x * kSqrt21, back_y * kSqrt21, out_i, out_j);
}

// ============================================================================
// Hexagon Corner Generation
// ============================================================================

/**
 * Generate the 6 corners of a hexagon given its center and radius.
 *
 * @param cx              Center X coordinate
 * @param cy              Center Y coordinate
 * @param radius          Distance from center to each corner
 * @param rotation_deg    Rotation offset in degrees (0 = flat-top, 30 = pointy-top)
 * @param out_x           Output array of 6 X coordinates
 * @param out_y           Output array of 6 Y coordinates
 */
inline void generate_hex_corners(double cx, double cy, double radius,
                                  double rotation_deg,
                                  double* out_x, double* out_y) {
    double rotation_rad = rotation_deg * kDegToRad;

    // Generate 6 vertices starting at top, counter-clockwise
    for (int k = 0; k < 6; ++k) {
        double angle = kPiOver2 + rotation_rad + k * kPiOver3;
        out_x[k] = cx + radius * std::cos(angle);
        out_y[k] = cy + radius * std::sin(angle);
    }
}

// ============================================================================
// Aperture Scale Factors
// ============================================================================

/**
 * Get the cumulative scale factor for a given aperture and resolution.
 *
 * @param aperture    Aperture type (3, 4, or 7)
 * @param resolution  Grid resolution level
 * @return            Scale factor to multiply coordinates by
 */
inline double aperture_scale(int aperture, int resolution) {
    switch (aperture) {
        case 3: return std::pow(kSqrt3, resolution);
        case 4: return std::pow(2.0, resolution);
        case 7: return std::pow(kSqrt7, resolution);
        default: return 1.0;
    }
}

/**
 * Get the substrate scale multiplier for the current orientation class.
 *
 * @param aperture    Aperture type
 * @param resolution  Grid resolution level
 * @return            Additional scale factor for substrate coordinates
 */
inline double substrate_multiplier(int aperture, int resolution) {
    bool is_even = (resolution % 2 == 0);

    switch (aperture) {
        case 3:
            // Rotation Class I (even): no extra scale
            // Rotation Class II (odd): sqrt(3) substrate
            return is_even ? 1.0 : kSqrt3;
        case 4:
            // Always Rotation Class I, no extra scale
            return 1.0;
        case 7:
            // Rotation Class III-A (even): sqrt(7) substrate
            // Rotation Class III-B (odd): sqrt(21) substrate
            return is_even ? kSqrt7 : kSqrt21;
        default:
            return 1.0;
    }
}

/**
 * Get the rotation offset in degrees for hexagon corners.
 *
 * @param aperture    Aperture type
 * @param resolution  Grid resolution level
 * @return            Rotation offset in degrees
 */
inline double corner_rotation_deg(int aperture, int resolution) {
    bool is_even = (resolution % 2 == 0);

    switch (aperture) {
        case 3:
            // Rotation Class I: 0 deg, Rotation Class II: 30 deg
            return is_even ? 0.0 : 30.0;
        case 4:
            // Always Rotation Class I
            return 0.0;
        case 7:
            // Rotation Class III-A: ~19.1 deg, Rotation Class III-B: ~49.1 deg
            return is_even ? kAp7RotDeg : (kAp7RotDeg + 30.0);
        default:
            return 0.0;
    }
}

} // namespace hexify

#endif // HEXIFY_GRID_MATH_H
