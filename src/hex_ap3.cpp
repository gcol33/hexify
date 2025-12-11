// hex_ap3.cpp - Aperture-3 hierarchical hexagon grid
//
// Implementation based on axial/cube coordinate systems for hexagonal grids.
// Mathematical foundation from:
//   - Sahr, White, Kimerling (2003) "Geodesic Discrete Global Grid Systems"
//   - Red Blob Games hexagonal grid reference (public domain concepts)
//
// The approach uses cube coordinates (q, r, s where q + r + s = 0) which
// provide elegant rounding to nearest hex center.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "hex_ap3.h"
#include "hex_coord.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace hexify {

namespace {

// Convert Cartesian (x, y) to cube coordinates for flat-top hex grid
CubeCoord cartesian_to_cube_flat(double x, double y) {
    return cartesian_to_cube(x, y, kSqrt3);
}

// Convert cube coordinates back to Cartesian for flat-top hex
void cube_to_cartesian_flat(double q, double r, double& x, double& y) {
    cube_to_cartesian(q, r, x, y, kSin60);
}

// Rotate point by angle (radians)
void rotate_point(double& x, double& y, double angle_rad) {
    double c = std::cos(angle_rad);
    double s = std::sin(angle_rad);
    double nx = x * c - y * s;
    double ny = x * s + y * c;
    x = nx;
    y = ny;
}

} // anonymous namespace

// ============================================================================
// Class I Quantization (Even Resolutions)
// ============================================================================
//
// Class I hex grids have flat-top orientation. We quantize by:
// 1. Convert Cartesian to cube coordinates
// 2. Round to nearest hex center
// 3. Extract (i, j) offset coordinates

static void quantify_class1(double x, double y, long long& out_i, long long& out_j) {
    CubeCoord cube = cartesian_to_cube_flat(x, y);
    cube.round_to_nearest();

    // Our cube-to-offset mapping: i = q, j = r
    out_i = static_cast<long long>(cube.q);
    out_j = static_cast<long long>(cube.r);
}

static void center_class1(long long i, long long j, double& x, double& y) {
    cube_to_cartesian_flat(static_cast<double>(i), static_cast<double>(j), x, y);
}

// ============================================================================
// Class II Quantization (Odd Resolutions)
// ============================================================================
//
// Class II hex grids are rotated 30° relative to Class I. The standard approach:
// 1. Rotate input point by -30° to align with a virtual Class I grid
// 2. Quantize in that rotated space
// 3. Express result in a finer "substrate" grid (scaled by sqrt(3))
//
// This gives coordinates compatible with hierarchical ISEA grids.

static void quantify_class2(double x, double y, long long& out_i, long long& out_j) {
    constexpr double angle = -kPi / 6.0;  // -30 degrees

    // Rotate to surrogate Class I orientation
    double rx = x, ry = y;
    rotate_point(rx, ry, angle);

    // Quantize in surrogate grid
    long long sur_i, sur_j;
    quantify_class1(rx, ry, sur_i, sur_j);

    // Get surrogate center
    double sur_x, sur_y;
    center_class1(sur_i, sur_j, sur_x, sur_y);

    // Rotate center back to original orientation
    rotate_point(sur_x, sur_y, -angle);

    // Scale to substrate grid (sqrt(3) finer) and re-quantize
    double sub_x = sur_x * kSqrt3;
    double sub_y = sur_y * kSqrt3;

    quantify_class1(sub_x, sub_y, out_i, out_j);
}

static void center_class2(long long i, long long j, double& x, double& y) {
    // Class II centers use substrate (Class I) positioning
    center_class1(i, j, x, y);
}

// ============================================================================
// Public API
// ============================================================================

void hex_quantify_ap3(double tx, double ty, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantify_ap3: resolution must be >= 0");
    }

    bool is_class1 = (resolution % 2 == 0);
    double scale = std::pow(kSqrt3, resolution);

    double scaled_x = tx * scale;
    double scaled_y = ty * scale;

    if (is_class1) {
        quantify_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        quantify_class2(scaled_x, scaled_y, out_i, out_j);
    }
}

void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap3: resolution must be >= 0");
    }

    bool is_class1 = (resolution % 2 == 0);
    double x, y;

    if (is_class1) {
        center_class1(i, j, x, y);
        double scale = std::pow(kSqrt3, resolution);
        out_cx = x / scale;
        out_cy = y / scale;
    } else {
        center_class2(i, j, x, y);
        // Substrate is sqrt(3) finer than the resolution scale
        double scale = std::pow(kSqrt3, resolution + 1);
        out_cx = x / scale;
        out_cy = y / scale;
    }
}

void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap3: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap3(i, j, resolution, cx, cy);

    double scale = std::pow(kSqrt3, resolution);
    double scaled_radius = hex_radius / scale;

    // Class II hexagons are rotated 30° relative to Class I
    bool is_class1 = (resolution % 2 == 0);
    double rotation_offset = is_class1 ? 0.0 : (kPi / 6.0);

    // Generate 6 vertices starting at top, counter-clockwise
    for (int k = 0; k < 6; ++k) {
        double angle = (kPi / 2.0) + rotation_offset + k * (kPi / 3.0);
        out_x[k] = cx + scaled_radius * std::cos(angle);
        out_y[k] = cy + scaled_radius * std::sin(angle);
    }
}

} // namespace hexify
