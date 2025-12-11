// hex_ap4.cpp - Aperture-4 hexagon grid
//
// Aperture 4 always uses Class I (flat-top) hexagons with scale factor 2.
// Implementation based on cube coordinate systems for hexagonal grids.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "hex_ap4.h"
#include "hex_coord.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

constexpr double kAp4Scale = 2.0;  // Aperture 4 scale factor

// Class I quantization (aperture 4 is always Class I)
void quantify_class1(double x, double y, long long& out_i, long long& out_j) {
    CubeCoord cube = cartesian_to_cube(x, y, kSqrt3);
    cube.round_to_nearest();
    out_i = static_cast<long long>(cube.q);
    out_j = static_cast<long long>(cube.r);
}

void center_class1(long long i, long long j, double& x, double& y) {
    cube_to_cartesian(static_cast<double>(i), static_cast<double>(j), x, y, kSin60);
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

void hex_quantify_ap4(double tx, double ty, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantify_ap4: resolution must be >= 0");
    }

    // Aperture 4: scale factor is 2^resolution
    double scale = std::pow(kAp4Scale, resolution);

    double scaled_x = tx * scale;
    double scaled_y = ty * scale;

    quantify_class1(scaled_x, scaled_y, out_i, out_j);
}

void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap4: resolution must be >= 0");
    }

    double x, y;
    center_class1(i, j, x, y);

    double scale = std::pow(kAp4Scale, resolution);
    out_cx = x / scale;
    out_cy = y / scale;
}

void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius,
                     double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap4: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap4(i, j, resolution, cx, cy);

    double scale = std::pow(kAp4Scale, resolution);
    double scaled_radius = hex_radius / scale;

    // Flat-top hexagon (Class I orientation)
    for (int k = 0; k < 6; ++k) {
        double angle = (kPi / 2.0) + k * (kPi / 3.0);
        out_x[k] = cx + scaled_radius * std::cos(angle);
        out_y[k] = cy + scaled_radius * std::sin(angle);
    }
}

} // namespace hexify
