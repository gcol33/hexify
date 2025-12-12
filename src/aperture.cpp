// aperture.cpp - Unified aperture implementation for ISEA hexagonal grids
//
// Implements quantization, center computation, and corner generation for
// aperture 3, 4, and 7 grids. All aperture-specific logic is consolidated here.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "aperture.h"
#include "grid_math.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

// ============================================================================
// Unified API Implementation
// ============================================================================

void hex_quantize(double icosa_triangle_x, double icosa_triangle_y, int aperture, int resolution,
                  long long& out_i, long long& out_j) {
    switch (aperture) {
        case 3:
            hex_quantize_ap3(icosa_triangle_x, icosa_triangle_y, resolution, out_i, out_j);
            break;
        case 4:
            hex_quantize_ap4(icosa_triangle_x, icosa_triangle_y, resolution, out_i, out_j);
            break;
        case 7:
            hex_quantize_ap7(icosa_triangle_x, icosa_triangle_y, resolution, out_i, out_j);
            break;
        default:
            throw std::runtime_error("hex_quantize: aperture must be 3, 4, or 7");
    }
}

void hex_center(long long i, long long j, int aperture, int resolution,
                double& out_cx, double& out_cy) {
    switch (aperture) {
        case 3:
            hex_center_ap3(i, j, resolution, out_cx, out_cy);
            break;
        case 4:
            hex_center_ap4(i, j, resolution, out_cx, out_cy);
            break;
        case 7:
            hex_center_ap7(i, j, resolution, out_cx, out_cy);
            break;
        default:
            throw std::runtime_error("hex_center: aperture must be 3, 4, or 7");
    }
}

void hex_corners(long long i, long long j, int aperture, int resolution,
                 double hex_radius, double* out_x, double* out_y) {
    switch (aperture) {
        case 3:
            hex_corners_ap3(i, j, resolution, hex_radius, out_x, out_y);
            break;
        case 4:
            hex_corners_ap4(i, j, resolution, hex_radius, out_x, out_y);
            break;
        case 7:
            hex_corners_ap7(i, j, resolution, hex_radius, out_x, out_y);
            break;
        default:
            throw std::runtime_error("hex_corners: aperture must be 3, 4, or 7");
    }
}

// ============================================================================
// Aperture 3 Implementation
// ============================================================================
//
// Aperture-3 grids alternate between Rotation Class I and Rotation Class II:
//   - Even resolutions: Rotation Class I (flat-top, 0 degrees)
//   - Odd resolutions:  Rotation Class II (pointy-top, 30 degrees)
//
// Cell count formula: N = 10 * 3^res + 2  (includes 12 pentagons)

void hex_quantize_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantize_ap3: resolution must be >= 0");
    }

    double scale = aperture_scale(3, resolution);
    double grid_x = icosa_triangle_x * scale;
    double grid_y = icosa_triangle_y * scale;

    bool is_even_res = (resolution % 2 == 0);
    if (is_even_res) {
        quantize_rotation_classI(grid_x, grid_y, out_i, out_j);
    } else {
        quantize_rotation_classII(grid_x, grid_y, out_i, out_j);
    }
}

void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap3: resolution must be >= 0");
    }

    center_rotation_classI(i, j, out_cx, out_cy);

    double scale = aperture_scale(3, resolution);
    double substrate_mult = substrate_multiplier(3, resolution);
    double total_scale = scale * substrate_mult;

    out_cx /= total_scale;
    out_cy /= total_scale;
}

void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap3: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap3(i, j, resolution, cx, cy);

    double scale = aperture_scale(3, resolution);
    double scaled_radius = hex_radius / scale;
    double rotation_deg = corner_rotation_deg(3, resolution);

    generate_hex_corners(cx, cy, scaled_radius, rotation_deg, out_x, out_y);
}

// ============================================================================
// Aperture 4 Implementation
// ============================================================================
//
// Aperture-4 always uses Rotation Class I (flat-top) hexagons with scale factor 2.
// Unlike aperture-3, there is no alternation between rotation classes.
//
// Cell count formula: N = 10 * 4^res + 2  (includes 12 pentagons)

void hex_quantize_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantize_ap4: resolution must be >= 0");
    }

    double scale = aperture_scale(4, resolution);
    double grid_x = icosa_triangle_x * scale;
    double grid_y = icosa_triangle_y * scale;

    quantize_rotation_classI(grid_x, grid_y, out_i, out_j);
}

void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap4: resolution must be >= 0");
    }

    center_rotation_classI(i, j, out_cx, out_cy);

    double scale = aperture_scale(4, resolution);
    out_cx /= scale;
    out_cy /= scale;
}

void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap4: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap4(i, j, resolution, cx, cy);

    double scale = aperture_scale(4, resolution);
    double scaled_radius = hex_radius / scale;

    generate_hex_corners(cx, cy, scaled_radius, 0.0, out_x, out_y);
}

// ============================================================================
// Aperture 7 Implementation
// ============================================================================
//
// Aperture-7 uses Rotation Class III hexagons rotated by arctan(sqrt(3/7)) = ~19.1 degrees.
// Two variants alternate by resolution:
//   - Rotation Class III-A (even res): Rotation Class I base + 19.1 deg
//   - Rotation Class III-B (odd res):  Rotation Class II base + 19.1 deg
//
// Cell count formula: N = 10 * 7^res + 2  (includes 12 pentagons)

void hex_quantize_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    if (resolution < 0) {
        throw std::runtime_error("hex_quantize_ap7: resolution must be >= 0");
    }

    double scale = aperture_scale(7, resolution);
    double grid_x = icosa_triangle_x * scale;
    double grid_y = icosa_triangle_y * scale;

    bool is_even_res = (resolution % 2 == 0);
    if (is_even_res) {
        quantize_rotation_classIII_A(grid_x, grid_y, out_i, out_j);
    } else {
        quantize_rotation_classIII_B(grid_x, grid_y, out_i, out_j);
    }
}

void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    if (resolution < 0) {
        throw std::runtime_error("hex_center_ap7: resolution must be >= 0");
    }

    center_rotation_classI(i, j, out_cx, out_cy);

    double scale = aperture_scale(7, resolution);
    double substrate_mult = substrate_multiplier(7, resolution);
    double total_scale = scale * substrate_mult;

    out_cx /= total_scale;
    out_cy /= total_scale;
}

void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    if (resolution < 0) {
        throw std::runtime_error("hex_corners_ap7: resolution must be >= 0");
    }

    double cx, cy;
    hex_center_ap7(i, j, resolution, cx, cy);

    double scale = aperture_scale(7, resolution);
    double scaled_radius = hex_radius / scale;
    double rotation_deg = corner_rotation_deg(7, resolution);

    generate_hex_corners(cx, cy, scaled_radius, rotation_deg, out_x, out_y);
}

} // namespace hexify
