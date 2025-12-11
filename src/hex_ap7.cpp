// hex_ap7.cpp - Aperture-7 hierarchical hexagon grid
//
// Aperture 7 uses Class III hexagons with alternating variants:
//   - Class III-I (even res): substrate is sqrt(7)× finer
//   - Class III-II (odd res): substrate is sqrt(21)× finer
//
// Implementation based on cube coordinate systems for hexagonal grids.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "hex_ap7.h"
#include "hex_coord.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// Class I quantization
void quantify_class1(double x, double y, long long& out_i, long long& out_j) {
    CubeCoord cube = cartesian_to_cube(x, y, kSqrt3);
    cube.round_to_nearest();
    out_i = static_cast<long long>(cube.q);
    out_j = static_cast<long long>(cube.r);
}

void center_class1(long long i, long long j, double& x, double& y) {
    cube_to_cartesian(static_cast<double>(i), static_cast<double>(j), x, y, kSin60);
}

// Class II quantization (30° rotated)
void quantify_class2(double x, double y, long long& out_i, long long& out_j) {
    constexpr double angle = -kPi / 6.0;  // -30°
    double c = std::cos(angle);
    double s = std::sin(angle);

    double rx = x * c - y * s;
    double ry = x * s + y * c;

    long long sur_i, sur_j;
    quantify_class1(rx, ry, sur_i, sur_j);

    double sur_x, sur_y;
    center_class1(sur_i, sur_j, sur_x, sur_y);

    double back_x = sur_x * c + sur_y * s;
    double back_y = -sur_x * s + sur_y * c;

    quantify_class1(back_x * kSqrt3, back_y * kSqrt3, out_i, out_j);
}

// ============================================================================
// Class III Quantization
// ============================================================================

// Class III-I (even resolutions): Class I surrogate rotated by kAp7RotDeg
void quantify_class3i(double x, double y, long long& out_i, long long& out_j) {
    const double angle = -kAp7RotDeg * kPi / 180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);

    // Rotate to surrogate
    double rx = x * c - y * s;
    double ry = x * s + y * c;

    // Quantize in Class I surrogate
    long long sur_i, sur_j;
    quantify_class1(rx, ry, sur_i, sur_j);

    // Get surrogate center and rotate back
    double sur_x, sur_y;
    center_class1(sur_i, sur_j, sur_x, sur_y);

    double back_x = sur_x * c + sur_y * s;
    double back_y = -sur_x * s + sur_y * c;

    // Scale to substrate (sqrt(7)× finer)
    quantify_class1(back_x * kSqrt7, back_y * kSqrt7, out_i, out_j);
}

// Class III-II (odd resolutions): Class II surrogate rotated by kAp7RotDeg
void quantify_class3ii(double x, double y, long long& out_i, long long& out_j) {
    const double angle = -kAp7RotDeg * kPi / 180.0;
    double c = std::cos(angle);
    double s = std::sin(angle);

    // Rotate to surrogate
    double rx = x * c - y * s;
    double ry = x * s + y * c;

    // Quantize in Class II surrogate
    long long sur_i, sur_j;
    quantify_class2(rx, ry, sur_i, sur_j);

    // Get surrogate center
    double sur_x, sur_y;
    center_class1(sur_i, sur_j, sur_x, sur_y);

    // Rotate back
    double back_x = sur_x * c + sur_y * s;
    double back_y = -sur_x * s + sur_y * c;

    // Scale to substrate (sqrt(7)× finer than Class II base)
    quantify_class1(back_x * kSqrt7, back_y * kSqrt7, out_i, out_j);
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

void hex_quantify_ap7(double tx, double ty, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantify_ap7: resolution must be >= 0");
    }

    double scale = std::pow(kSqrt7, resolution);
    double scaled_x = tx * scale;
    double scaled_y = ty * scale;

    bool is_class3i = (resolution % 2 == 0);

    if (is_class3i) {
        quantify_class3i(scaled_x, scaled_y, out_i, out_j);
    } else {
        quantify_class3ii(scaled_x, scaled_y, out_i, out_j);
    }
}

void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap7: resolution must be >= 0");
    }

    double x, y;
    center_class1(i, j, x, y);

    double base_scale = std::pow(kSqrt7, resolution);

    bool is_class3i = (resolution % 2 == 0);
    double substrate_multiplier = is_class3i ? kSqrt7 : kSqrt21;

    double substrate_scale = base_scale * substrate_multiplier;

    out_cx = x / substrate_scale;
    out_cy = y / substrate_scale;
}

void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap7: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap7(i, j, resolution, cx, cy);

    double scale = std::pow(kSqrt7, resolution);
    double scaled_radius = hex_radius / scale;

    bool is_class3i = (resolution % 2 == 0);
    double rotation_offset = is_class3i
        ? (kAp7RotDeg * kPi / 180.0)
        : ((kAp7RotDeg + 30.0) * kPi / 180.0);

    for (int k = 0; k < 6; ++k) {
        double angle = (kPi / 2.0) + rotation_offset + k * (kPi / 3.0);
        out_x[k] = cx + scaled_radius * std::cos(angle);
        out_y[k] = cy + scaled_radius * std::sin(angle);
    }
}

} // namespace hexify
