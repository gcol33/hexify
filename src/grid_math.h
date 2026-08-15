// grid_math.h - Shared hexagonal grid mathematics
//
// This module provides the fundamental mathematical operations for hexagonal
// discrete global grid systems (DGGS). It consolidates operations that were
// previously duplicated across aperture-specific files.
//
// ============================================================================
// GRID ORIENTATION
// ============================================================================
//
// Refining a hexagonal lattice by aperture a replaces it with a sublattice of
// index a. For a hexagonal lattice such a sublattice exists only when a is the
// norm of an Eisenstein integer z = m + n*w, where w = exp(i*pi/3) and the norm
// is m^2 + m*n + n^2. The refinement scales lengths by sqrt(a) = |z| and
// rotates the lattice by arg(z):
//
//   aperture 3:  z = 1 + w     norm 3, arg 30 deg
//   aperture 4:  z = 2         norm 4, arg 0 deg
//   aperture 7:  z = 2 + w     norm 7, arg 19.106... deg  (kAp7RotDeg)
//                z = 1 + 2w    norm 7, arg 40.894... deg  (the mirror generator)
//
// Refining sends a lattice L to z^-1 L, so a grid's orientation is the product
// of the conjugates of its refinement steps' generators. The classical rotation
// classes (Sahr et al. 2003) are the small cases: norm 1 is Class I (0 degrees,
// flat-top), norm 3 is Class II (30 degrees, pointy-top).
//
// Aperture-7 steps alternate between the two norm-7 generators. Their product
// 7w is a rational integer times a unit, so orientation returns to the base
// every two aperture-7 levels rather than drifting, and the substrate stays
// bounded. Aperture 4's generator is the rational integer 2, which scales the
// lattice without rotating it.
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
// The substrate is sqrt(norm(z)) times finer than the grid itself, so the
// factors that occur are 1, sqrt(3), sqrt(7) and sqrt(21).
//
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_GRID_MATH_H
#define HEXIFY_GRID_MATH_H

#include "cube_coordinates.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>
#include <vector>

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
    // Guard against NaN/Inf inputs to avoid UB in CubeCoord::round_to_nearest()'s
    // std::llround() call. Matches coordinate_transforms.cpp::quantize_class1.
    if (!std::isfinite(x) || !std::isfinite(y)) {
        out_i = 0;
        out_j = 0;
        return;
    }

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
 * Linear scale factor contributed by one refinement step.
 *
 * @param aperture    Aperture type (3, 4, or 7)
 */
inline double aperture_linear_scale(int aperture) {
    switch (aperture) {
        case 3: return kSqrt3;
        case 4: return 2.0;
        case 7: return kSqrt7;
        default: throw std::runtime_error("aperture_linear_scale: aperture must be 3, 4, or 7");
    }
}

/**
 * Get the cumulative scale factor for a given aperture and resolution.
 *
 * @param aperture    Aperture type (3, 4, or 7)
 * @param resolution  Grid resolution level
 * @return            Scale factor to multiply coordinates by
 */
inline double aperture_scale(int aperture, int resolution) {
    return std::pow(aperture_linear_scale(aperture), resolution);
}

// ============================================================================
// Grid Form (Scale + Lattice Orientation)
// ============================================================================

/**
 * How far a grid has been refined from the base grid, and how its lattice sits
 * relative to the Class I (unrotated) lattice.
 *
 * `scale` is the product of sqrt(aperture) over the refinement steps.
 * (m, n) is the canonical generator z = m + n*w of the orientation.
 */
struct HexGridForm {
    double scale;
    long long m;
    long long n;
};

/** Norm of the Eisenstein integer m + n*w. */
inline long long eisenstein_norm(long long m, long long n) {
    return m * m + m * n + n * n;
}

/** Multiply m + n*w by (c + d*w), using w^2 = w - 1. */
inline void eisenstein_multiply(long long& m, long long& n, long long c, long long d) {
    long long new_m = m * c - n * d;
    long long new_n = m * d + n * c + n * d;
    m = new_m;
    n = new_n;
}

/**
 * Reduce a generator to the unique associate representing its orientation:
 * divide out the rational-integer content (which scales without rotating),
 * then multiply by units until arg(z) lies in [0, 60) degrees.
 *
 * Multiplication by the unit w sends (m, n) to (-n, m + n).
 */
inline void eisenstein_canonical(long long& m, long long& n) {
    long long a = m < 0 ? -m : m;
    long long b = n < 0 ? -n : n;
    while (b != 0) {
        long long t = a % b;
        a = b;
        b = t;
    }
    if (a > 1) {
        m /= a;
        n /= a;
    }

    // arg(z) is in [0, 60) exactly when n >= 0 and m > 0; six units to try.
    for (int k = 0; k < 6 && !(n >= 0 && m > 0); ++k) {
        long long new_m = -n;
        long long new_n = m + n;
        m = new_m;
        n = new_n;
    }
}

/**
 * Apply one refinement step to a grid orientation.
 *
 * Refining by a generator z sends the lattice L to z^-1 L, so orientation
 * composes with the conjugate of z (conj(m + n*w) = (m + n) - n*w):
 *
 *   aperture 3:  conj(1 + w)  = 2 - w
 *   aperture 4:  conj(2)      = 2
 *   aperture 7:  conj(1 + 2w) = 3 - 2w   and   conj(2 + w) = 3 - w
 *
 * @param ap7_step  Number of aperture-7 steps already applied. Successive
 *                  aperture-7 steps alternate between the two norm-7
 *                  generators, so orientation returns to the base every two
 *                  levels. The first step is the one that puts odd pure
 *                  aperture-7 resolutions at kAp7RotDeg.
 */
inline void orientation_refine(long long& m, long long& n, int aperture, int ap7_step) {
    switch (aperture) {
        case 3: eisenstein_multiply(m, n, 2, -1); break;
        case 4: eisenstein_multiply(m, n, 2, 0); break;
        case 7:
            if (ap7_step % 2 == 0) {
                eisenstein_multiply(m, n, 3, -2);
            } else {
                eisenstein_multiply(m, n, 3, -1);
            }
            break;
        default:
            throw std::runtime_error("orientation_refine: aperture must be 3, 4, or 7");
    }
    eisenstein_canonical(m, n);
}

/** Form of a pure single-aperture grid at the given resolution. */
inline HexGridForm hex_form_pure(int aperture, int resolution) {
    if (resolution < 0) {
        throw std::runtime_error("hex_form_pure: resolution must be >= 0");
    }
    HexGridForm form;
    form.scale = aperture_scale(aperture, resolution);
    form.m = 1;
    form.n = 0;
    for (int r = 0; r < resolution; ++r) {
        orientation_refine(form.m, form.n, aperture, r);
    }
    return form;
}

/**
 * Form of a grid built from a mixed aperture sequence.
 *
 * ap_seq[0] names the base grid, which is not the result of a refinement;
 * the refinement steps are ap_seq[1..].
 */
inline HexGridForm hex_form_sequence(const std::vector<int>& ap_seq) {
    if (ap_seq.empty()) {
        throw std::runtime_error("hex_form_sequence: ap_seq cannot be empty");
    }
    HexGridForm form;
    form.scale = 1.0;
    form.m = 1;
    form.n = 0;
    aperture_linear_scale(ap_seq[0]);  // validate the base entry
    int ap7_step = 0;
    for (size_t k = 1; k < ap_seq.size(); ++k) {
        int aperture = ap_seq[k];
        form.scale *= aperture_linear_scale(aperture);
        orientation_refine(form.m, form.n, aperture, ap7_step);
        if (aperture == 7) {
            ap7_step++;
        }
    }
    return form;
}

/** Cosine and sine of the lattice rotation. */
inline void form_rotation(const HexGridForm& form, double& cos_a, double& sin_a) {
    if (form.n == 0) {
        cos_a = 1.0;
        sin_a = 0.0;
        return;
    }
    if (form.m == 1 && form.n == 1) {
        cos_a = kCos30;
        sin_a = kSin30;
        return;
    }
    double re = static_cast<double>(form.m) + 0.5 * static_cast<double>(form.n);
    double im = kSin60 * static_cast<double>(form.n);
    double len = std::sqrt(re * re + im * im);
    cos_a = re / len;
    sin_a = im / len;
}

/** Lattice rotation in degrees, in [0, 60). */
inline double form_rotation_deg(const HexGridForm& form) {
    if (form.n == 0) return 0.0;
    if (form.m == 1 && form.n == 1) return 30.0;
    double re = static_cast<double>(form.m) + 0.5 * static_cast<double>(form.n);
    double im = kSin60 * static_cast<double>(form.n);
    return std::atan2(im, re) * kRadToDeg;
}

/** How much finer the Class I substrate is than the grid itself. */
inline double form_substrate(const HexGridForm& form) {
    long long norm = eisenstein_norm(form.m, form.n);
    if (norm == 1) return 1.0;
    return std::sqrt(static_cast<double>(norm));
}

/**
 * Quantize a point to the nearest cell of a grid with the given form.
 *
 * Unrotated grids quantize directly. Rotated grids use the surrogate-substrate
 * pattern: rotate into the Class I surrogate frame, quantize there, take the
 * surrogate centre, rotate back, then re-quantize on the substrate.
 *
 * @param out_i  Output: column index in substrate coordinates
 * @param out_j  Output: row index in substrate coordinates
 */
inline void quantize_form(const HexGridForm& form, double x, double y,
                          long long& out_i, long long& out_j) {
    double grid_x = x * form.scale;
    double grid_y = y * form.scale;

    if (form.n == 0) {
        quantize_rotation_classI(grid_x, grid_y, out_i, out_j);
        return;
    }

    double cos_a, sin_a;
    form_rotation(form, cos_a, sin_a);

    double sur_x = grid_x * cos_a + grid_y * sin_a;
    double sur_y = -grid_x * sin_a + grid_y * cos_a;

    long long sur_i, sur_j;
    quantize_rotation_classI(sur_x, sur_y, sur_i, sur_j);

    double cen_x, cen_y;
    center_rotation_classI(sur_i, sur_j, cen_x, cen_y);

    double back_x = cen_x * cos_a - cen_y * sin_a;
    double back_y = cen_x * sin_a + cen_y * cos_a;

    double substrate = form_substrate(form);
    quantize_rotation_classI(back_x * substrate, back_y * substrate, out_i, out_j);
}

/** Centre of cell (i, j), whose coordinates are on the substrate. */
inline void center_form(const HexGridForm& form, long long i, long long j,
                        double& out_cx, double& out_cy) {
    center_rotation_classI(i, j, out_cx, out_cy);
    double total_scale = form.scale * form_substrate(form);
    out_cx /= total_scale;
    out_cy /= total_scale;
}

/** The 6 corners of cell (i, j). */
inline void corners_form(const HexGridForm& form, long long i, long long j,
                         double hex_radius, double* out_x, double* out_y) {
    double cx, cy;
    center_form(form, i, j, cx, cy);
    generate_hex_corners(cx, cy, hex_radius / form.scale, form_rotation_deg(form),
                         out_x, out_y);
}

} // namespace hexify

#endif // HEXIFY_GRID_MATH_H
