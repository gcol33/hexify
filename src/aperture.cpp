// aperture.cpp - Unified aperture implementation for ISEA hexagonal grids
//
// Implements quantization, center computation, and corner generation for
// aperture 3, 4, and 7 grids. Scale and lattice orientation come from the
// shared grid-form model in grid_math.h, which mixed sequences also use.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "aperture.h"
#include "grid_math.h"

namespace hexify {

// ============================================================================
// Unified API Implementation
// ============================================================================

void hex_quantize(double icosa_triangle_x, double icosa_triangle_y, int aperture, int resolution,
                  long long& out_i, long long& out_j) {
    HexGridForm form = hex_form_pure(aperture, resolution);
    quantize_form(form, icosa_triangle_x, icosa_triangle_y, out_i, out_j);
}

void hex_center(long long i, long long j, int aperture, int resolution,
                double& out_cx, double& out_cy) {
    HexGridForm form = hex_form_pure(aperture, resolution);
    center_form(form, i, j, out_cx, out_cy);
}

void hex_corners(long long i, long long j, int aperture, int resolution,
                 double hex_radius, double* out_x, double* out_y) {
    HexGridForm form = hex_form_pure(aperture, resolution);
    corners_form(form, i, j, hex_radius, out_x, out_y);
}

// ============================================================================
// Aperture 3 (ISEA3H)
// ============================================================================
//
// Scale factor sqrt(3) per level; orientation alternates between Class I (even
// resolutions) and Class II (odd resolutions).
// Cell count: N = 10 * 3^res + 2  (includes 12 pentagons)

void hex_quantize_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    hex_quantize(icosa_triangle_x, icosa_triangle_y, 3, resolution, out_i, out_j);
}

void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    hex_center(i, j, 3, resolution, out_cx, out_cy);
}

void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    hex_corners(i, j, 3, resolution, hex_radius, out_x, out_y);
}

// ============================================================================
// Aperture 4 (ISEA4H)
// ============================================================================
//
// Scale factor 2 per level. The generator is the rational integer 2, so every
// resolution stays Class I.
// Cell count: N = 10 * 4^res + 2  (includes 12 pentagons)

void hex_quantize_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    hex_quantize(icosa_triangle_x, icosa_triangle_y, 4, resolution, out_i, out_j);
}

void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    hex_center(i, j, 4, resolution, out_cx, out_cy);
}

void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    hex_corners(i, j, 4, resolution, hex_radius, out_x, out_y);
}

// ============================================================================
// Aperture 7 (ISEA7H)
// ============================================================================
//
// Scale factor sqrt(7) per level. Successive levels rotate by kAp7RotDeg in
// alternating directions, so even resolutions are Class I and odd resolutions
// sit at kAp7RotDeg on a sqrt(7) substrate.
// Cell count: N = 10 * 7^res + 2  (includes 12 pentagons)

void hex_quantize_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j) {
    hex_quantize(icosa_triangle_x, icosa_triangle_y, 7, resolution, out_i, out_j);
}

void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy) {
    hex_center(i, j, 7, resolution, out_cx, out_cy);
}

void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y) {
    hex_corners(i, j, 7, resolution, hex_radius, out_x, out_y);
}

} // namespace hexify
