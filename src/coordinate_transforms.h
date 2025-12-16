// coordinate_transforms.h - Convert between ISEA DGGS coordinate systems
//
// Coordinate transformations for ISEA DGGS grids.
// Copyright (c) 2024 hexify authors. MIT License.
//
// ============================================================================
// COORDINATE SYSTEM GLOSSARY
// ============================================================================
//
// This module handles transformations between three coordinate systems used
// in ISEA Discrete Global Grid Systems (DGGS):
//
// 1. Icosahedral Triangle Coordinates (from Snyder projection)
//    - Variables: icosa_triangle_face, icosa_triangle_x, icosa_triangle_y
//    - icosa_triangle_face: Triangle/face number (0-19), one of the 20 icosahedral faces
//    - icosa_triangle_x, icosa_triangle_y: Projected coordinates within that face, typically [0, 1]
//    - This is what the Snyder forward projection produces from lon/lat
//
// 2. Quad XY (continuous quad coordinates)
//    - Variables: quad, quad_x, quad_y
//    - quad: Quad number (0-11), pairs of triangles forming diamond shapes
//    - quad_x, quad_y: Continuous floating-point coordinates within the quad
//    - Intermediate representation between icosa triangle coords and quad IJ
//
// 3. Quad IJ (quantized cell indices)
//    - Variables: quad, i, j
//    - quad: Quad number (0-11)
//    - i, j: Integer cell indices within the quad at a given resolution
//    - Cell IDs are derived from this
//
// ICOSA TRIANGLE TO QUAD MAPPING:
// -------------------------------
// The 20 triangular faces are grouped into 12 quads:
//   - Quad 0:     North polar region (vertex)
//   - Quads 1-5:  Upper hemisphere rhombi (each contains 2 triangles)
//   - Quads 6-10: Lower hemisphere rhombi (each contains 2 triangles)
//   - Quad 11:    South polar region (vertex)
//
// ============================================================================

#ifndef HEXIFY_COORDINATE_TRANSFORMS_H
#define HEXIFY_COORDINATE_TRANSFORMS_H

namespace hexify {

// Convert from icosahedral triangle coordinates to quad XY coordinates
// This applies triTable rotation and translation
//
// Parameters:
//   icosa_triangle_face: Triangle/face number (0-19)
//   icosa_triangle_x, icosa_triangle_y: Projected triangle coordinates
//   out_quad: Output quad number (0-11)
//   out_quad_x, out_quad_y: Output continuous quad coordinates
void icosa_tri_to_quad_xy(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int& out_quad, double& out_quad_x, double& out_quad_y);

// Convert from quad XY to quad IJ (cell indices)
// This quantizes continuous coords to integer cell indices
//
// Parameters:
//   quad: Quad number (0-11)
//   quad_x, quad_y: Continuous quad coordinates
//   aperture: Grid aperture (3, 4, or 7)
//   resolution: Grid resolution level
//   out_quad: Output quad (may change due to edge overflow)
//   out_i, out_j: Output integer cell indices
void quad_xy_to_ij(int quad, double quad_x, double quad_y,
                   int aperture, int resolution,
                   int& out_quad, long long& out_i, long long& out_j);

// Full pipeline: icosa triangle coords → quad IJ
// Combines icosa_tri_to_quad_xy and quad_xy_to_ij
void icosa_tri_to_quad_ij(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int aperture, int resolution,
                          int& out_quad, long long& out_i, long long& out_j);

// Inverse: quad IJ → quad XY (for computing cell centers)
void quad_ij_to_xy(int quad, long long i, long long j,
                   int aperture, int resolution,
                   double& out_quad_x, double& out_quad_y);

// Inverse: quad XY → icosa triangle coords (throws on invalid region)
void quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                          int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y);

// Inverse: quad XY → icosa triangle coords (returns false on invalid region)
bool try_quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                              int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y);

// Get maxI/maxJ for a given aperture and resolution
long long get_max_ij(int aperture, int resolution);

// Edge handling: map edge cells to adjacent quads
// Returns true if coord was on edge and was adjusted
bool handle_edge_overflow(int& quad, long long& i, long long& j,
                          int aperture, int resolution);

} // namespace hexify

#endif // HEXIFY_COORDINATE_TRANSFORMS_H
