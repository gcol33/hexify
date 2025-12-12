// aperture.h - Unified aperture interface for ISEA hexagonal grids
//
// This module provides a single API for all supported aperture types (3, 4, 7).
// Each aperture defines how cells subdivide at each resolution level:
//
//   Aperture 3 (ISEA3H):
//     - Scale factor: sqrt(3) per level
//     - Alternates between Rotation Class I (even res) and Rotation Class II (odd res)
//     - Cell count: 10 * 3^res + 2
//
//   Aperture 4 (ISEA4H):
//     - Scale factor: 2 per level
//     - Always uses Rotation Class I (flat-top)
//     - Cell count: 10 * 4^res + 2
//
//   Aperture 7 (ISEA7H):
//     - Scale factor: sqrt(7) per level
//     - Uses Rotation Class III variants (rotated ~19.1 degrees)
//     - Cell count: 10 * 7^res + 2
//
// For mixed aperture sequences (e.g., ISEA43H), see aperture_sequence.h
//
// Copyright (c) 2024 hexify authors. MIT License.

#pragma once

namespace hexify {

// ============================================================================
// Unified Aperture API
// ============================================================================

/**
 * Quantize a point in triangle coordinates to (i, j) cell indices.
 *
 * @param icosa_triangle_x  X coordinate in triangle space (from Snyder projection)
 * @param icosa_triangle_y  Y coordinate in triangle space
 * @param aperture          Aperture type (3, 4, or 7)
 * @param resolution        Grid resolution level (>= 0)
 * @param out_i             Output: column index
 * @param out_j             Output: row index
 * @throws std::runtime_error if aperture or resolution is invalid
 */
void hex_quantize(double icosa_triangle_x, double icosa_triangle_y, int aperture, int resolution,
                  long long& out_i, long long& out_j);

/**
 * Get the center coordinates of a cell in triangle space.
 *
 * @param i           Column index
 * @param j           Row index
 * @param aperture    Aperture type (3, 4, or 7)
 * @param resolution  Grid resolution level (>= 0)
 * @param out_cx      Output: center X coordinate
 * @param out_cy      Output: center Y coordinate
 * @throws std::runtime_error if aperture or resolution is invalid
 */
void hex_center(long long i, long long j, int aperture, int resolution,
                double& out_cx, double& out_cy);

/**
 * Get the 6 corner coordinates of a cell in triangle space.
 *
 * @param i           Column index
 * @param j           Row index
 * @param aperture    Aperture type (3, 4, or 7)
 * @param resolution  Grid resolution level (>= 0)
 * @param hex_radius  Radius from center to each corner
 * @param out_x       Output: array of 6 X coordinates
 * @param out_y       Output: array of 6 Y coordinates
 * @throws std::runtime_error if aperture or resolution is invalid
 */
void hex_corners(long long i, long long j, int aperture, int resolution,
                 double hex_radius, double* out_x, double* out_y);

// ============================================================================
// Aperture-Specific Functions (for backwards compatibility)
// ============================================================================

// Aperture 3
void hex_quantize_ap3(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j);
void hex_center_ap3(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);
void hex_corners_ap3(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y);

// Aperture 4
void hex_quantize_ap4(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j);
void hex_center_ap4(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);
void hex_corners_ap4(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y);

// Aperture 7
void hex_quantize_ap7(double icosa_triangle_x, double icosa_triangle_y, int resolution,
                      long long& out_i, long long& out_j);
void hex_center_ap7(long long i, long long j, int resolution,
                    double& out_cx, double& out_cy);
void hex_corners_ap7(long long i, long long j, int resolution,
                     double hex_radius, double* out_x, double* out_y);

} // namespace hexify
