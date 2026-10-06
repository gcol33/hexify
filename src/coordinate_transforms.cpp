// coordinate_transforms.cpp - Convert between the coordinate systems of a solid's DGGS
//
// ============================================================================
// COORDINATE TRANSFORMATION FLOW
// ============================================================================
//
// This file implements the coordinate transformations between three systems:
//
//     +---------------------+     +---------------+     +-----------+
//     |   Icosa Triangle    | --> |    Quad XY    | --> |  Quad IJ  |
//     | (icosa_triangle_    |     | (quad,        |     | (quad,    |
//     |  face, _x, _y)      |     |  quad_x,      |     |  i, j)    |
//     +---------------------+     |  quad_y)      |     +-----------+
//          |                      +---------------+           |
//          v                            |                     v
//     Icosahedral face            Quad continuous        Quad integer
//     coordinates                 coordinates            cell indices
//
// Icosa Triangle: Output from Snyder forward projection
//   - icosa_triangle_face: Face index of the solid
//   - icosa_triangle_x, icosa_triangle_y: Normalized coords within triangle [0,1]
//
// Quad XY: Quad with continuous (double) coordinates
//   - quad: Quad index (the first and last are the vertex quads, see below)
//   - quad_x, quad_y: Continuous position within quad
//
// Quad IJ: Quad with integer cell indices (used for cell ID computation)
//   - quad: Same as Quad XY
//   - i, j: Integer cell coordinates (resolution-dependent)
//
// ============================================================================
// SOLID GEOMETRY
// ============================================================================
//
// The faces of the solid (polyhedron.h) pair into diamonds across shared
// edges, one diamond per quad. On the icosahedron the 20 faces form 12 quads:
//
//           Quad 0 (vertex 0)
//                  /\
//                 /  \
//           +----+----+----+----+----+
//           | Q1 | Q2 | Q3 | Q4 | Q5 |  <- Upper quads 1-5
//           +----+----+----+----+----+
//           | Q6 | Q7 | Q8 | Q9 |Q10 |  <- Lower quads 6-10
//           +----+----+----+----+----+
//                 \  /
//                  \/
//           Quad 11 (antipode of vertex 0)
//
// and on the octahedron the 8 faces form 6: four diamonds, each a northern and
// a southern face, and the two vertices they meet at. Quad q's origin corner
// is the solid's vertex q, so the first and last quads are the two degenerate
// quads: each holds the single cell at a vertex no diamond starts at. Those
// vertices sit at the poles only under a pole-aligned orientation; under the
// ISEA default (vertex 0 at 11.25, 58.28) they do not, so their lon/lat comes
// from folding (0, 0) through the quad frame like any other cell.
//
// Every table the conversions below read -- which quad a face lies in, which
// face holds each region around a quad's origin, and the map across each quad
// edge -- is derived from the solid's face list (SolidTopology).
//
// ============================================================================
// QUANTIZATION CLASSES
// ============================================================================
//
// Different apertures use different quantization schemes:
//
//   Aperture 3:
//     - Even resolutions: Class I (aligned hexagons)
//     - Odd resolutions: Class II (rotated hexagons)
//
//   Aperture 4:
//     - All resolutions: Class I (aligned hexagons)
//
//   Aperture 7:
//     - Even resolutions: Class III-I
//     - Odd resolutions: Class III-II
//
// Mathematical foundation from Sahr et al. publications on ISEA grids.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "coordinate_transforms.h"
#include "cube_coordinates.h"
#include "grid_math.h"
#include "ijk_coordinates.h"
#include "index_z7.h"
#include "constants.h"
#include "polyhedron.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace hexify {

namespace {

// ============================================================================
// Rotation Helper
// ============================================================================

void rotate_60deg_ccw(double& x, double& y, int n_rotations) {
    // Each 60° counter-clockwise rotation: [cos(60) -sin(60); sin(60) cos(60)]
    // cos(60°) = 0.5, sin(60°) = sqrt(3)/2
    constexpr double c60 = 0.5;
    constexpr double s60 = kSin60;

    n_rotations = ((n_rotations % 6) + 6) % 6;  // Normalize to 0-5

    for (int i = 0; i < n_rotations; ++i) {
        double nx = c60 * x - s60 * y;  // counter-clockwise: x*cos - y*sin
        double ny = s60 * x + c60 * y;  // counter-clockwise: x*sin + y*cos
        x = nx;
        y = ny;
    }
}

// ============================================================================
// Hex Quantization - Precise hexagonal grid rounding
// ============================================================================
// This implementation handles all edge cases at hexagon boundaries correctly
// using a decision-tree approach that carefully handles the fractional parts
// of the continuous coordinates. This is more robust than simple cube-coordinate
// rounding at cell boundaries.

// ============================================================================
// Boundary Classification for Hex Quantization
// ============================================================================
//
// The unit cell is divided into 6 regions based on fractional coordinates (frac_i, frac_j).
// Each region determines the (delta_i, delta_j) offset from the base cell (floor_i, floor_j).
//
// The regions form a hexagonal Voronoi partition:
//   - Region A: frac_i < 1/3, frac_j < (1+frac_i)/2    -> (0, 0)
//   - Region B: frac_i < 1/3, frac_j >= (1+frac_i)/2   -> (0, 1)
//   - Region C: 1/3 <= frac_i < 1/2                    -> complex boundary (see below)
//   - Region D: 1/2 <= frac_i < 2/3                    -> complex boundary (see below)
//   - Region E: frac_i >= 2/3, frac_j < frac_i/2       -> (1, 0)
//   - Region F: frac_i >= 2/3, frac_j >= frac_i/2      -> (1, 1)
//
// For regions C and D, the i-offset depends on whether frac_j falls in the
// "middle band" between two linear thresholds.

// Classify which boundary region based on fractional coords
// Returns: 0=A, 1=B, 2=C_lower, 3=C_upper, 4=C_mid, 5=D_lower, 6=D_upper, 7=D_mid, 8=E, 9=F
inline int classify_hex_boundary(double frac_i, double frac_j) {
    if (frac_i < 1.0/3.0) {
        return (frac_j < (1.0 + frac_i) / 2.0) ? 0 : 1;  // A or B
    }
    if (frac_i < 0.5) {
        // Region C: thresholds at (1-frac_i) and (2*frac_i)
        double lower_threshold = 1.0 - frac_i;
        double upper_threshold = 2.0 * frac_i;
        if (frac_j < lower_threshold) return 2;       // C_lower: j=floor_j
        if (frac_j >= upper_threshold) return 3;      // C_upper: j=floor_j+1
        return 4;                                      // C_mid: i=floor_i+1
    }
    if (frac_i < 2.0/3.0) {
        // Region D: thresholds at (2*frac_i-1) and (1-frac_i)
        double lower_threshold = 2.0 * frac_i - 1.0;
        double upper_threshold = 1.0 - frac_i;
        if (frac_j <= lower_threshold) return 5;      // D_lower: j=floor_j, i=floor_i+1
        if (frac_j >= upper_threshold) return 6;      // D_upper: j=floor_j+1, i=floor_i+1
        return 7;                                      // D_mid: i=floor_i
    }
    return (frac_j < frac_i / 2.0) ? 8 : 9;  // E or F
}

// Lookup table: boundary_region -> (delta_i, delta_j) offset
// Indexed by classify_hex_boundary() return value
static const int kBoundaryOffset[10][2] = {
    {0, 0},  // 0: Region A
    {0, 1},  // 1: Region B
    {0, 0},  // 2: Region C_lower (j=floor_j)
    {0, 1},  // 3: Region C_upper (j=floor_j+1)
    {1, 0},  // 4: Region C_mid - special: j depends on frac_j < (1-frac_i)
    {1, 0},  // 5: Region D_lower (j=floor_j)
    {1, 1},  // 6: Region D_upper (j=floor_j+1)
    {0, 0},  // 7: Region D_mid - special: j depends on frac_j < (1-frac_i)
    {1, 0},  // 8: Region E
    {1, 1},  // 9: Region F
};

// Fold i-coordinate across x-axis when x was negative
inline long long fold_i_negative_x(long long i, long long j) {
    if ((j % 2) == 0) {
        long long axis = j / 2;
        return i - 2 * (i - axis);
    } else {
        long long axis = (j + 1) / 2;
        return i - (2 * (i - axis) + 1);
    }
}

// Class I (flat-top) quantization
void quantize_class1(double x, double y, long long& out_i, long long& out_j) {
    // Guard against NaN/Inf inputs to avoid undefined behavior in integer cast
    if (!std::isfinite(x) || !std::isfinite(y)) {
        out_i = 0;
        out_j = 0;
        return;
    }

    // Work in positive quadrant
    double abs_x = std::fabs(x);
    double abs_y = std::fabs(y);

    // Convert to fractional hex indices
    double idx_j = abs_y / kSin60;
    double idx_i = abs_x + idx_j / 2.0;

    // Integer (floor) and fractional parts
    long long floor_i = static_cast<long long>(idx_i);
    long long floor_j = static_cast<long long>(idx_j);
    double frac_i = idx_i - floor_i;
    double frac_j = idx_j - floor_j;

    // Classify and look up base offset
    int region = classify_hex_boundary(frac_i, frac_j);
    long long delta_i = kBoundaryOffset[region][0];
    long long delta_j = kBoundaryOffset[region][1];

    // Handle special cases where j depends on secondary threshold
    if (region == 4) {  // C_mid
        delta_j = (frac_j < (1.0 - frac_i)) ? 0 : 1;
    } else if (region == 7) {  // D_mid
        delta_j = (frac_j < (1.0 - frac_i)) ? 0 : 1;
    }

    long long i_result = floor_i + delta_i;
    long long j_result = floor_j + delta_j;

    // Fold back to original quadrant
    if (x < 0.0) {
        i_result = fold_i_negative_x(i_result, j_result);
    }
    if (y < 0.0) {
        i_result = i_result - (2 * j_result + 1) / 2;
        j_result = -j_result;
    }

    out_i = i_result;
    out_j = j_result;
}

// Class I inverse: (i,j) to (x,y)
void inv_quantize_class1(long long i, long long j, double& x, double& y) {
    cube_to_cartesian(static_cast<double>(i), static_cast<double>(j), x, y, kSin60);
}

// Class II (pointy-top / 30° rotated) quantization
void quantize_class2(double x, double y, long long& out_i, long long& out_j) {
    constexpr double angle = -kPi / 6.0;  // -30°
    double c = std::cos(angle);
    double s = std::sin(angle);

    // Rotate to surrogate Class I orientation
    double rx = x * c - y * s;
    double ry = x * s + y * c;

    // Quantize in surrogate
    long long sur_i, sur_j;
    quantize_class1(rx, ry, sur_i, sur_j);

    // Get surrogate center and rotate back
    double sur_x, sur_y;
    inv_quantize_class1(sur_i, sur_j, sur_x, sur_y);

    double back_x = sur_x * c + sur_y * s;  // Rotate +30°
    double back_y = -sur_x * s + sur_y * c;

    // Scale to substrate and re-quantize
    quantize_class1(back_x * kSqrt3, back_y * kSqrt3, out_i, out_j);
}


} // anonymous namespace

// ============================================================================
// Aperture 7: exact-integer surrogate machinery (matches DGGRID / H3)
// ============================================================================
// The "surrogate" is hexify's canonical aperture-7 cell coordinate: the exact
// integer IJK of the resolution-r cell. It is obtained by a clean, unrotated
// Class I quantization of the shared quad_xy frame at the Class I substrate
// scale (7^numClassI = sqrt(7)^effectiveRes), DGGRID's edgeTable quad
// canonicalization, and, for odd resolutions, one exact aperture-7 coarsen
// (upAp7r). This replaces the earlier floating-point-rotation surrogate, whose
// re-quantization rounded boundary cells to a neighbour and diverged from the
// exact integer grid.

namespace {

// Folds one term of an edge map into a running sum: the first nonzero term
// starts the sum, later ones are added or subtracted.
template <typename T>
inline void edge_map_term(int coef, T value, bool& started, T& sum) {
    if (coef == 0) return;
    if (!started) {
        sum = (coef == 1) ? value : (coef == -1) ? -value : static_cast<T>(coef) * value;
        started = true;
    } else if (coef == 1) {
        sum = sum + value;
    } else if (coef == -1) {
        sum = sum - value;
    } else {
        sum = sum + static_cast<T>(coef) * value;
    }
}

// One coordinate of an edge map: k[0] * topEdge + k[1] * along + k[2] * d,
// summed in that order.
template <typename T>
inline T edge_map_coord(const int k[3], T topEdge, T along, T d) {
    bool started = false;
    T sum = 0;
    edge_map_term(k[0], topEdge, started, sum);
    edge_map_term(k[1], along, started, sum);
    edge_map_term(k[2], d, started, sum);
    return started ? sum : T(0);
}

// Reassign an out-of-box quad coordinate (i,j) to the quad that owns it, as
// DGGRID's DgQ2DDtoIConverter does through its edge table, with the solid's
// edge maps (QuadEdgeMap). topEdge = maxI + 1 = maxJ + 1. Every edge map is
// affine in (i, j), so the same maps carry a continuous coordinate across a
// quad edge when topEdge is the quad's side. Returns false when the
// coordinate lies beyond both far edges or below both near edges, which a
// cell centre never does: an integer centre beyond both far edges is the far
// vertex, and is moved there. A vertex quad has no box to leave.
template <typename T>
bool canonicalize_q2d(T topEdge, int& quadNum, T& i, T& j) {
    const T maxI = topEdge - 1, maxJ = topEdge - 1;
    const bool integral = std::is_integral<T>::value;

    bool underI = i < 0, underJ = j < 0;
    bool overI = integral ? i > maxI : i >= topEdge;
    bool overJ = integral ? j > maxJ : j >= topEdge;
    int numOver = (int)underI + (int)underJ + (int)overI + (int)overJ;
    if (!numOver) return true;

    const SolidTopology& t = topo();
    if (t.is_pole(quadNum)) return false;

    if (overI && overJ) {
        if (!integral) return false;
        quadNum = t.corner[quadNum][kCornerFar];
        i = 0; j = 0;
        return true;
    }
    if (numOver > 1) return false;

    const int e = underI ? kEdgeLeft : underJ ? kEdgeDown : overI ? kEdgeRight : kEdgeUp;
    const QuadEdgeMap& m = t.edge[quadNum][e];
    const T along = (e == kEdgeLeft || e == kEdgeRight) ? j : i;
    const T d = (e == kEdgeLeft) ? i : (e == kEdgeDown) ? j
              : (e == kEdgeRight) ? i - topEdge : j - topEdge;
    if (m.pole >= 0 && along == 0) {
        quadNum = m.pole;
        i = 0; j = 0;
        return true;
    }
    const T ni = edge_map_coord(m.k[0], topEdge, along, d);
    const T nj = edge_map_coord(m.k[1], topEdge, along, d);
    quadNum = m.quad;
    i = ni;
    j = nj;
    return true;
}

void dggrid_canonicalize_q2di(long long topEdge, int& quadNum,
                              long long& i, long long& j) {
    canonicalize_q2d<long long>(topEdge, quadNum, i, j);
}

} // anonymous namespace

// Substrate steps along a quad edge after n3, n4 and n7 refinement steps of
// apertures 3, 4 and 7, in any order. The square of the substrate scale is the
// product of the apertures times the norm of the grid form's generator. Two
// aperture-3 steps compose to 3 times a unit, two aperture-7 steps to 7 times
// a unit and an aperture-4 step is 2, so that norm is 3 after an odd number of
// aperture-3 steps, 7 after an odd number of aperture-7 steps, 21 after both
// and 1 otherwise. The edge is therefore 2^n4 * 3^ceil(n3/2) * 7^ceil(n7/2).
static long long edge_dim_of_steps(int n3, int n4, int n7) {
    long long d = 1;
    for (int k = 0; k < n4; ++k) d *= 2;
    for (int k = 0; k < (n3 + 1) / 2; ++k) d *= 3;
    for (int k = 0; k < (n7 + 1) / 2; ++k) d *= 7;
    return d;
}

long long quad_edge_dim(int aperture, int resolution) {
    switch (aperture) {
        case 3: return edge_dim_of_steps(resolution, 0, 0);
        case 4: return edge_dim_of_steps(0, resolution, 0);
        case 7: return edge_dim_of_steps(0, 0, resolution);
        default: throw std::runtime_error("quad_edge_dim: aperture must be 3, 4, or 7");
    }
}

long long quad_edge_dim(const std::vector<int>& ap_seq) {
    int n[8] = {0};
    for (size_t k = 1; k < ap_seq.size(); ++k) {
        int a = ap_seq[k];
        if (a != 3 && a != 4 && a != 7) {
            throw std::runtime_error("quad_edge_dim: aperture must be 3, 4, or 7");
        }
        n[a]++;
    }
    return edge_dim_of_steps(n[3], n[4], n[7]);
}

void ap7_substrate_to_surrogate_ijk(long long sub_i, long long sub_j, int resolution,
                                    long long& sur_i, long long& sur_j) {
    if (resolution % 2 == 0) { sur_i = sub_i; sur_j = sub_j; return; }
    z7::IVec3D v(sub_i, sub_j, 0);
    v.upAp7r();
    z7::IVec2D a(v);
    sur_i = a.i();
    sur_j = a.j();
}

void ap7_surrogate_to_substrate_ijk(long long sur_i, long long sur_j, int resolution,
                                    long long& sub_i, long long& sub_j) {
    if (resolution % 2 == 0) { sub_i = sur_i; sub_j = sur_j; return; }
    z7::IVec3D v(sur_i, sur_j, 0);
    v.downAp7r();
    z7::IVec2D a(v);
    sub_i = a.i();
    sub_j = a.j();
}

void ap7_nearest_centre(double px, double py, long long sub_i, long long sub_j,
                        int resolution, long long& ctr_i, long long& ctr_j) {
    if (resolution % 2 == 0) { ctr_i = sub_i; ctr_j = sub_j; return; }
    // At odd resolutions the cells are the Voronoi regions of a rotated
    // sublattice of the substrate, one centre per seven substrate points. The
    // point lies within one substrate circumradius of (sub_i, sub_j) and within
    // sqrt(7) of them of its own centre, so that centre is at most 2.1 lattice
    // spacings from (sub_i, sub_j): inside the 5 x 5 block of indices around it.
    double best = -1.0;
    for (long long di = -2; di <= 2; ++di) {
        for (long long dj = -2; dj <= 2; ++dj) {
            long long ci = sub_i + di, cj = sub_j + dj;
            long long sur_i, sur_j, back_i, back_j;
            ap7_substrate_to_surrogate_ijk(ci, cj, resolution, sur_i, sur_j);
            ap7_surrogate_to_substrate_ijk(sur_i, sur_j, resolution, back_i, back_j);
            if (back_i != ci || back_j != cj) continue;
            double cx, cy;
            inv_quantize_class1(ci, cj, cx, cy);
            double d = (cx - px) * (cx - px) + (cy - py) * (cy - py);
            if (best < 0.0 || d < best) {
                best = d;
                ctr_i = ci;
                ctr_j = cj;
            }
        }
    }
}

uint64_t ap7_surrogate_to_quad_index(long long sur_i, long long sur_j, int resolution) {
    const long long S = quad_edge_dim(7, resolution);
    long long u, v;
    ap7_surrogate_to_substrate_ijk(sur_i, sur_j, resolution, u, v);
    if (resolution % 2 == 0) {
        return static_cast<uint64_t>(u * S + v);
    }
    // v is fixed modulo 7 once u is known (2u + v = 0 mod 7), so v / 7 names the
    // centre within its row on its own.
    return static_cast<uint64_t>(u * (S / 7) + v / 7);
}

void ap7_quad_index_to_surrogate(uint64_t index, int resolution,
                                 long long& sur_i, long long& sur_j) {
    const long long S = quad_edge_dim(7, resolution);
    const long long idx = static_cast<long long>(index);
    long long u, v;
    if (resolution % 2 == 0) {
        u = idx / S;
        v = idx % S;
    } else {
        const long long rows = S / 7;
        u = idx / rows;
        v = 7 * (idx % rows) + (((-2 * u) % 7) + 7) % 7;
    }
    ap7_substrate_to_surrogate_ijk(u, v, resolution, sur_i, sur_j);
}

bool substrate_ij_canonicalize(int& quad, long long& i, long long& j,
                               long long top_edge) {
    int q = quad;
    long long ci = i, cj = j;

    dggrid_canonicalize_q2di(top_edge, q, ci, cj);
    if (ci < 0 || ci >= top_edge || cj < 0 || cj >= top_edge) {
        return false;
    }

    quad = q;
    i = ci;
    j = cj;
    return true;
}

bool quad_ij_canonicalize(int& quad, long long& i, long long& j,
                          int aperture, int resolution) {
    const long long edge = quad_edge_dim(aperture, resolution);
    if (aperture != 7) {
        return substrate_ij_canonicalize(quad, i, j, edge);
    }

    long long ci, cj;
    ap7_surrogate_to_substrate_ijk(i, j, resolution, ci, cj);
    if (!substrate_ij_canonicalize(quad, ci, cj, edge)) {
        return false;
    }

    ap7_substrate_to_surrogate_ijk(ci, cj, resolution, i, j);
    return true;
}

bool ap7_surrogate_in_quad(long long sur_i, long long sur_j, int resolution) {
    const long long S = quad_edge_dim(7, resolution);
    long long u, v;
    ap7_surrogate_to_substrate_ijk(sur_i, sur_j, resolution, u, v);
    return u >= 0 && u < S && v >= 0 && v < S;
}

void surrogate_ij_to_quad_xy_ap7(long long sur_i, long long sur_j, int resolution,
                                  double& out_quad_x, double& out_quad_y) {
    long long S = quad_edge_dim(7, resolution);
    long long sub_i, sub_j;
    ap7_surrogate_to_substrate_ijk(sur_i, sur_j, resolution, sub_i, sub_j);
    double cx, cy;
    inv_quantize_class1(sub_i, sub_j, cx, cy);
    out_quad_x = cx / static_cast<double>(S);
    out_quad_y = cy / static_cast<double>(S);
}

// ============================================================================
// Public API Implementation
// ============================================================================

void icosa_tri_to_quad_xy(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int& out_quad, double& out_quad_x, double& out_quad_y) {
    const SolidTopology& t = topo();
    if (icosa_triangle_face < 0 || icosa_triangle_face >= t.n_faces) {
        throw std::runtime_error("icosa_tri_to_quad_xy: face out of range for the solid");
    }

    const FacePlacement& mapping = t.placement[icosa_triangle_face];

    out_quad = mapping.quad;
    out_quad_x = icosa_triangle_x;
    out_quad_y = icosa_triangle_y;

    // Apply rotation then translation
    rotate_60deg_ccw(out_quad_x, out_quad_y, mapping.rotations);
    out_quad_x -= mapping.offset_x;
    out_quad_y -= mapping.offset_y;
}

void quad_xy_to_ij(int quad, double quad_x, double quad_y,
                   int aperture, int resolution,
                   int& out_quad, long long& out_i, long long& out_j) {

    // Aperture 7: exact-integer route. Clean unrotated Class I quantization at
    // the substrate scale, DGGRID edgeTable quad canonicalization (an out-of-box
    // coordinate belongs to the neighbouring quad), then (odd res) one exact
    // aperture-7 coarsen -- yielding the exact resolution-r cell IJK. This keeps
    // forward/inverse geometry consistent and replaces the float-rotation Class
    // III quantization, which rounded boundary cells.
    if (aperture == 7) {
        if (resolution == 0) {
            // Resolution 0: one cell per vertex. A point nearest a corner of
            // the quad other than its origin goes to that corner's vertex,
            // which the edge maps name; the z7 hierarchy is empty at res 0.
            quantize_class1(quad_x, quad_y, out_i, out_j);
            out_quad = quad;
            dggrid_canonicalize_q2di(1, out_quad, out_i, out_j);
            return;
        }
        long long S = quad_edge_dim(7, resolution);
        double px = quad_x * static_cast<double>(S);
        double py = quad_y * static_cast<double>(S);
        long long sub_i, sub_j;
        quantize_class1(px, py, sub_i, sub_j);
        long long ctr_i, ctr_j;
        ap7_nearest_centre(px, py, sub_i, sub_j, resolution, ctr_i, ctr_j);
        // Canonicalize the cell CENTRE. A cell on a quad edge covers points on
        // both sides of it, so canonicalizing the sampled point would give that
        // one cell an address in either quad. Its centre lies in exactly one
        // quad and so fixes the owner.
        out_quad = quad;
        dggrid_canonicalize_q2di(S, out_quad, ctr_i, ctr_j);
        ap7_substrate_to_surrogate_ijk(ctr_i, ctr_j, resolution, out_i, out_j);
        return;
    }

    // Compute scale factor
    double scale;
    if (aperture == 3) {
        scale = std::pow(kSqrt3, resolution);
    } else if (aperture == 4) {
        scale = std::pow(2.0, resolution);
    } else {
        throw std::runtime_error("quad_xy_to_ij: unsupported aperture");
    }

    double scaled_x = quad_x * scale;
    double scaled_y = quad_y * scale;

    // Select quantization based on aperture and grid class
    if (aperture == 4 || (aperture == 3 && resolution % 2 == 0)) {
        // Class I quantization
        quantize_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        // Class II quantization (aperture 3 odd resolutions)
        quantize_class2(scaled_x, scaled_y, out_i, out_j);
    }

    // A cell centre beyond the quad belongs to the neighbouring quad. A Class II
    // cell edge runs along the icosahedron edge, so a point on or just across
    // that edge can quantize to a centre more than one row outside the quad;
    // DGGRID's edgeTable map reassigns any such centre, not only the first row.
    out_quad = quad;
    dggrid_canonicalize_q2di(quad_edge_dim(aperture, resolution), out_quad, out_i, out_j);
}

void icosa_tri_to_quad_ij(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int aperture, int resolution,
                          int& out_quad, long long& out_i, long long& out_j) {
    int quad;
    double quad_x, quad_y;
    icosa_tri_to_quad_xy(icosa_triangle_face, icosa_triangle_x, icosa_triangle_y, quad, quad_x, quad_y);
    quad_xy_to_ij(quad, quad_x, quad_y, aperture, resolution, out_quad, out_i, out_j);
}

void quad_ij_to_xy(int quad, long long i, long long j,
                   int aperture, int resolution,
                   double& out_quad_x, double& out_quad_y) {

    // Aperture 7: exact-integer route, the inverse of the one quad_xy_to_ij()
    // takes. DGGRID's DgHexGrid2DS toggles Class III on every aperture-7 level
    // (DgHexGrid2DS.cpp), so even resolutions are an unrotated Class I grid with
    // no substrate and odd resolutions carry one aperture-7 level; both come out
    // as the integer divisor 7^numClassI that surrogate_ij_to_quad_xy_ap7()
    // applies.
    if (aperture == 7) {
        surrogate_ij_to_quad_xy_ap7(i, j, resolution, out_quad_x, out_quad_y);
        return;
    }

    double x, y;
    inv_quantize_class1(i, j, x, y);

    // Compute inverse scale accounting for substrate
    double scale;
    if (aperture == 3) {
        bool is_class1 = (resolution % 2 == 0);
        scale = is_class1
            ? std::pow(kSqrt3, resolution)
            : std::pow(kSqrt3, resolution + 1);  // Class II substrate
    } else if (aperture == 4) {
        scale = std::pow(2.0, resolution);
    } else {
        throw std::runtime_error("quad_ij_to_xy: unsupported aperture");
    }

    out_quad_x = x / scale;
    out_quad_y = y / scale;
}

void quad_xy_to_ij_mixed(int quad, double quad_x, double quad_y,
                         const HexGridForm& form, long long edge,
                         int& out_quad, long long& out_i, long long& out_j) {
    quantize_form(form, quad_x, quad_y, out_i, out_j);

    out_quad = quad;

    // (out_i, out_j) is the substrate coordinate of the nearest cell centre.
    // A lattice rotated off the substrate axes by an odd number of aperture-7
    // steps has no mirror symmetry across a quad edge, so the cells along an
    // edge straddle it and the nearest centre of a point inside the quad can
    // lie outside [0, edge]^2. That centre is a cell of the neighbouring
    // quad; the edge table moves it there, and sends a centre on the far edges
    // or vertices to the quad that owns it.
    dggrid_canonicalize_q2di(edge, out_quad, out_i, out_j);
}

// ============================================================================
// Sub-triangle Region Detection
// ============================================================================
//
// Divides the Quad XY coordinate space into 6 wedge-shaped regions emanating
// from the origin. The boundaries are lines at angles 0°, 60°, 120°, 180°,
// 240°, 300° from the positive x-axis. The key boundary is y = ±sqrt(3)*x.
//
//          Region 0 (Upper)
//             /\
//      Reg 5 /  \ Reg 1
//      -----+----+-----
//      Reg 4 \  / Reg 2
//             \/
//          Region 3 (Lower)
//
// Each region maps to a face of the solid around the quad's origin
// (SolidTopology::region).

// Check if point is at origin (within tolerance)
inline bool is_origin(double x, double y, double tol) {
    return std::fabs(x) <= tol && std::fabs(y) <= tol;
}

// Compute which of 6 sub-regions a Quad XY point falls into
// Uses 6-way wedge classification based on y = ±sqrt(3)*x boundaries
static int compute_subtriangle(double x, double y) {
    constexpr double tol = 1e-15;

    // Origin -> Region 1 (center/upper-right by convention)
    if (is_origin(x, y, tol)) return 1;

    // Pre-compute boundary lines: y = ±sqrt(3)*x with tolerance
    const double xs = kSqrt3 * x;
    const double xs_plus  = xs + tol;   // y = sqrt(3)*x + tol
    const double xs_minus = xs - tol;   // y = sqrt(3)*x - tol
    const double neg_xs_plus  = -xs + tol;  // y = -sqrt(3)*x + tol
    const double neg_xs_minus = -xs - tol;  // y = -sqrt(3)*x - tol

    // Region 0: Upper (above both diagonal lines)
    if (y >= neg_xs_minus && y > xs_plus) return 0;

    // Region 1: Upper-right (below y=sqrt(3)*x, above y=0)
    if (y <= xs_plus && y >= -tol) return 1;

    // Region 2: Lower-right (below y=0, above y=-sqrt(3)*x)
    if (y < -tol && y > neg_xs_plus) return 2;

    // Region 3: Lower (below both diagonal lines)
    if (y <= neg_xs_plus && y < xs_minus) return 3;

    // Region 4: Lower-left (above y=sqrt(3)*x, below y=0)
    if (y >= xs_minus && y < -tol) return 4;

    // Region 5: Upper-left (above y=0, below y=-sqrt(3)*x)
    if (y >= -tol && y < neg_xs_minus) return 5;

    // Fallback (should not occur for valid quad coordinates)
    return 1;
}

// Try to convert quad XY to icosa triangle coords. Returns true on success,
// false if the point is in an invalid region (e.g., outside the valid quad bounds).
bool quad_xy_canonicalize(int& quad, double& quad_x, double& quad_y) {
    const SolidTopology& t = topo();
    if (t.is_pole(quad)) return true;
    // The quad box is the unit rhombus of the Class I lattice basis.
    double v = quad_y / kSin60;
    double u = quad_x + v / 2.0;
    // Below a near edge the point is still in the fan of faces around the
    // quad's origin vertex, which try_quad_xy_to_icosa_tri() reads directly,
    // dropped sectors included.
    if (u < 1.0 && v < 1.0) return true;
    int q = quad;
    if (!canonicalize_q2d<double>(1.0, q, u, v)) return false;
    quad = q;
    quad_x = u - v / 2.0;
    quad_y = v * kSin60;
    return true;
}


bool try_quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                              int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y) {
    const SolidTopology& t = topo();
    if (quad < 0 || quad >= t.n_quads()) {
        throw std::invalid_argument("try_quad_xy_to_icosa_tri: quad out of range for the solid");
    }

    // Detect which of 6 sub-regions the point falls into
    int subTri = compute_subtriangle(quad_x, quad_y);

    const QuadRegion& triVal = t.region[quad][subTri];

    if (!triVal.keep) {
        // The region lies in the solid's angular deficit at the vertex
        return false;
    }

    out_icosa_triangle_face = triVal.face;

    // Apply inverse transformation:
    //   coord += trans
    //   coord.rotate(rot60 * -60.0)  // rotate by -60*rot60 degrees CCW
    out_icosa_triangle_x = quad_x + triVal.trans_x;
    out_icosa_triangle_y = quad_y + triVal.trans_y;

    // Rotate: rot60 * -60 degrees = -60 * rot60 degrees CCW
    // Which is the same as 60 * rot60 degrees CW
    // If rot60=4, rotate -240 degrees CCW = 120 degrees CW = -2 rotations of 60deg CCW
    rotate_60deg_ccw(out_icosa_triangle_x, out_icosa_triangle_y, -triVal.rot60);
    return true;
}

void quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                          int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y) {
    if (!try_quad_xy_to_icosa_tri(quad, quad_x, quad_y, out_icosa_triangle_face, out_icosa_triangle_x, out_icosa_triangle_y)) {
        throw std::runtime_error("quad_xy_to_icosa_tri: point in invalid region");
    }
}

} // namespace hexify
