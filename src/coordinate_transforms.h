// coordinate_transforms.h - Convert between the coordinate systems of a solid's DGGS
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
//    - icosa_triangle_face: Face number of the solid (0-19 on the icosahedron)
//    - icosa_triangle_x, icosa_triangle_y: Projected coordinates within that face, typically [0, 1]
//    - This is what the Snyder forward projection produces from lon/lat
//
// 2. Quad XY (continuous quad coordinates)
//    - Variables: quad, quad_x, quad_y
//    - quad: Quad number, pairs of faces forming diamond shapes (0-11 on the
//      icosahedron, 0-5 on the octahedron)
//    - quad_x, quad_y: Continuous floating-point coordinates within the quad
//    - Intermediate representation between icosa triangle coords and quad IJ
//
// 3. Quad IJ (quantized cell indices)
//    - Variables: quad, i, j
//    - quad: Quad number
//    - i, j: Integer cell indices within the quad at a given resolution
//    - Cell IDs are derived from this
//
// FACE TO QUAD MAPPING:
// ---------------------
// The faces of the solid pair into diamond quads, one per vertex the diamond
// starts at, and the two vertices no diamond starts at are single-cell quads:
// on the icosahedron
//   - Quad 0:     the vertex quad of vertex 0
//   - Quads 1-5:  upper rhombi (each contains 2 triangles)
//   - Quads 6-10: lower rhombi (each contains 2 triangles)
//   - Quad 11:    the vertex quad of vertex 11
//
// ============================================================================

#ifndef HEXIFY_COORDINATE_TRANSFORMS_H
#define HEXIFY_COORDINATE_TRANSFORMS_H

#include <vector>
#include <cstdint>

namespace hexify {

struct HexGridForm;

// Convert from icosahedral triangle coordinates to quad XY coordinates
// This applies triTable rotation and translation
//
// Parameters:
//   icosa_triangle_face: Face number of the solid
//   icosa_triangle_x, icosa_triangle_y: Projected triangle coordinates
//   out_quad: Output quad number
//   out_quad_x, out_quad_y: Output continuous quad coordinates
void icosa_tri_to_quad_xy(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int& out_quad, double& out_quad_x, double& out_quad_y);

// Convert from quad XY to quad IJ (cell indices)
// This quantizes continuous coords to integer cell indices
//
// Parameters:
//   quad: Quad number
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

// Move a quad-plane point lying past a far edge of its quad into the quad that
// owns it, in that quad's frame. False when the point lies past both far edges,
// beyond the far vertex. A point below a near edge stays, as the fan of faces
// around the quad's origin reads it; with `across_near_edges` it is carried
// across that edge into the quad there instead.
bool quad_xy_canonicalize(int& quad, double& quad_x, double& quad_y,
                          bool across_near_edges = false);

// Inverse: quad XY → icosa triangle coords (returns false on invalid region)
bool try_quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                              int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y);

// Substrate steps along a quad edge, which is the coordinate of the quad's
// far corner: exact for every aperture and resolution. A pure grid takes
// 'resolution' refinement steps of 'aperture'; a mixed sequence takes
// ap_seq[1..] (see aperture_sequence.h). For aperture 7 this is also the Class
// I substrate scale 7^ceil(resolution / 2) the surrogate is quantized at.
long long quad_edge_dim(int aperture, int resolution);
long long quad_edge_dim(const std::vector<int>& ap_seq);

// Aperture 7: exact-integer conversion between the Class I substrate IJK (what
// z7 operates on) and the resolution-r surrogate IJK (hexify's stored cell
// coordinate). Even resolutions are the identity; odd resolutions apply one
// exact aperture-7 level (upAp7r / downAp7r).
void ap7_substrate_to_surrogate_ijk(long long sub_i, long long sub_j, int resolution,
                                    long long& sur_i, long long& sur_j);
void ap7_surrogate_to_substrate_ijk(long long sur_i, long long sur_j, int resolution,
                                    long long& sub_i, long long& sub_j);

// The substrate coordinates of the centre of the aperture-7 cell holding the
// substrate-scaled point (px, py), given its nearest substrate point.
void ap7_nearest_centre(double px, double py, long long sub_i, long long sub_j,
                        int resolution, long long& ctr_i, long long& ctr_j);

// Aperture 7: dense cell index within a quad, and its inverse.
//
// A quad owns exactly 7^resolution aperture-7 cells: those whose centre falls in
// the half-open Class I substrate box [0, S)^2, S = quad_edge_dim(7, r). Even
// resolutions store that centre directly, so the index is its row-major position
// in the box. Odd resolutions store the coarsened surrogate, whose centre is a
// point of the aperture-7 sublattice {(u, v) : 2u + v = 0 (mod 7)}; each row
// therefore holds S/7 centres and the index counts those.
//
// The result spans [0, 7^resolution) with no gaps, so cell IDs run 1 ..
// (number of diamonds) * 7^resolution + 2 exactly as they do for apertures 3
// and 4. A quad with a fold axis numbers the centres it holds at their places
// in the box (quad_slot()).
uint64_t ap7_surrogate_to_quad_index(int quad, long long sur_i, long long sur_j,
                                     int resolution);
void ap7_quad_index_to_surrogate(int quad, uint64_t index, int resolution,
                                 long long& sur_i, long long& sur_j);

// Aperture 7: is this surrogate's centre held by the given diamond quad, i.e.
// inside the substrate box [0, S)^2 as quad_holds() reads it?
bool ap7_surrogate_in_quad(int quad, long long sur_i, long long sur_j, int resolution);

// A diamond quad with a fold axis holds the interior of the far edge across
// that axis (quad_holds()). Index strings write such a cell one quad edge
// back along the axis, on the near edge the quad gives up, and mark it. This
// moves a stored (i, j), in the aperture's own cell coordinate, back so when it
// lies on that far edge, and says whether it did.
bool quad_ij_from_far_edge(int quad, long long& i, long long& j,
                           int aperture, int resolution);

// The inverse: (i, j) moved one quad edge forward along the quad's fold axis.
void quad_ij_to_far_edge(int quad, long long& i, long long& j,
                         int aperture, int resolution);

// Re-express an (i, j) that has stepped outside its quad in the quad that owns
// it, via DGGRID's edge table. Coordinates are the aperture's own cell
// coordinate (the aperture-7 surrogate is expanded and coarsened around the
// call). Returns false when the coordinate lands outside every adjacent quad,
// which happens where the solid folds at a vertex.
bool quad_ij_canonicalize(int& quad, long long& i, long long& j,
                          int aperture, int resolution);

// The same walk on a bare substrate coordinate, given the coordinate of the
// quad's far corner. A mixed aperture sequence stores its cells on the
// substrate, so its quad edge is all the edge table needs. With
// `vertex_cells` false, for a lattice with no cell at the solid's vertices
// (Hex9), no point is moved to a vertex.
bool substrate_ij_canonicalize(int& quad, long long& i, long long& j,
                               long long top_edge, bool vertex_cells = true);

// Aperture 7: Inverse - surrogate IJ back to quad XY coordinates.
void surrogate_ij_to_quad_xy_ap7(long long sur_i, long long sur_j, int resolution,
                                  double& out_quad_x, double& out_quad_y);

// Mixed aperture sequence: quad XY -> quad IJ, for the grid of the given form
// (hex_form_sequence()) whose quad edge is 'edge' (quad_edge_dim()). The
// returned (i, j) are the substrate coordinates of the nearest cell centre,
// moved into the quad that owns it. The inverse is center_form().
void quad_xy_to_ij_mixed(int quad, double quad_x, double quad_y,
                         const HexGridForm& form, long long edge,
                         int& out_quad, long long& out_i, long long& out_j);


} // namespace hexify

#endif // HEXIFY_COORDINATE_TRANSFORMS_H
