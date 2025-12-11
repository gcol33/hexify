// tri_to_q2di.cpp - Convert PROJTRI (triangle coords) to Q2DI (quad integer coords)
//
// This implements the coordinate transformations needed for ISEA DGGS:
//   PROJTRI: (triangle_number, tx, ty) - projected triangle coordinates
//   Q2DD: (quad_number, qx, qy) - quad with continuous coordinates
//   Q2DI: (quad_number, i, j) - quad with integer cell indices
//
// The icosahedron has 20 triangular faces grouped into quads (pairs of triangles).
// Mathematical foundation from Sahr et al. publications on ISEA grids.
//
// Copyright (c) 2024 hexify authors. MIT License.

#include "tri_to_q2di.h"
#include "hex_coord.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// ============================================================================
// Triangle to Quad Mapping
// ============================================================================
//
// Icosahedron face layout (20 triangles → 12 quads):
//
//                     North Pole (Quad 0)
//                          /\
//                         /  \
//              +---------+----+---------+
//             /\   0    /\   1/\   2    /\
//            /  \      /  \  /  \      /  \     Upper cap
//           / 5  \    / 6  \/  7 \    / 8  \    (triangles 0-4 primary,
//          +------+--+------+------+--+------+   5-9 secondary)
//           \  9  /  \  10 /\  11 /  \  12 /
//            \  /    \  /  \  /    \  /         Lower cap
//             \/  13  \/  14\/  15  \/          (triangles 10-14 primary,
//              +------+------+------+            15-19 secondary)
//                     South Pole (Quad 11)
//
// Quads 1-5: upper hemisphere, Quads 6-10: lower hemisphere
// Each quad contains 2 triangles forming a rhombus shape.
//
// The 20 icosahedral faces are grouped into 12 quads:
//   - Quad 0: North pole (special)
//   - Quads 1-5: Upper hemisphere
//   - Quads 6-10: Lower hemisphere
//   - Quad 11: South pole (special)
//
// Each non-polar quad contains 2 triangles. The mapping specifies:
//   - Which quad a triangle belongs to
//   - Rotation and translation to align triangle coords with quad coords

struct TriangleMapping {
    int quad;           // Target quad (1-10 for regular quads)
    int sub_triangle;   // 0 = primary, 1 = secondary (rotated/translated)
    double offset_x;    // X offset after rotation
    double offset_y;    // Y offset after rotation
    int rotations;      // Number of 60° clockwise rotations
};

// Mapping table derived from ISEA icosahedron geometry
// Triangle indices 0-19 map to quads 1-10 (polar quads 0,11 handled separately)
const TriangleMapping kTriangleMap[20] = {
    // Upper cap triangles (0-4) -> quads 1-5, primary position
    {1, 0, 0.0, 0.0, 1},
    {2, 0, 0.0, 0.0, 1},
    {3, 0, 0.0, 0.0, 1},
    {4, 0, 0.0, 0.0, 1},
    {5, 0, 0.0, 0.0, 1},
    // Upper-middle triangles (5-9) -> quads 1-5, secondary position
    {1, 1, -0.5, -kSin60, 4},
    {2, 1, -0.5, -kSin60, 4},
    {3, 1, -0.5, -kSin60, 4},
    {4, 1, -0.5, -kSin60, 4},
    {5, 1, -0.5, -kSin60, 4},
    // Lower-middle triangles (10-14) -> quads 6-10, primary position
    {6,  0, 0.0, 0.0, 1},
    {7,  0, 0.0, 0.0, 1},
    {8,  0, 0.0, 0.0, 1},
    {9,  0, 0.0, 0.0, 1},
    {10, 0, 0.0, 0.0, 1},
    // Lower cap triangles (15-19) -> quads 6-10, secondary position
    {6,  1, -0.5, -kSin60, 4},
    {7,  1, -0.5, -kSin60, 4},
    {8,  1, -0.5, -kSin60, 4},
    {9,  1, -0.5, -kSin60, 4},
    {10, 1, -0.5, -kSin60, 4},
};

// ============================================================================
// Rotation Helper
// ============================================================================

void rotate_60deg_ccw(double& x, double& y, int n_rotations) {
    // Each 60° counter-clockwise rotation: [cos(60) -sin(60); sin(60) cos(60)]
    // cos(60°) = 0.5, sin(60°) = sqrt(3)/2
    // DGGRID uses counter-clockwise rotation for positive angles
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

// Class I (flat-top) quantization
void quantize_class1(double x, double y, long long& out_i, long long& out_j) {
    // Take absolute values for initial calculation
    double a1 = std::fabs(x);
    double a2 = std::fabs(y);

    // Reverse conversion to fractional indices
    double x2 = a2 / kSin60;
    double x1 = a1 + x2 / 2.0;

    // Get integer parts
    long long m1 = static_cast<long long>(x1);
    long long m2 = static_cast<long long>(x2);

    // Get fractional parts
    double r1 = x1 - m1;
    double r2 = x2 - m2;

    long long i_result, j_result;

    // Decision tree based on fractional parts
    if (r1 < 0.5) {
        if (r1 < 1.0/3.0) {
            if (r2 < (1.0 + r1) / 2.0) {
                i_result = m1;
                j_result = m2;
            } else {
                i_result = m1;
                j_result = m2 + 1;
            }
        } else {
            if (r2 < (1.0 - r1)) {
                j_result = m2;
            } else {
                j_result = m2 + 1;
            }

            if ((1.0 - r1) <= r2 && r2 < (2.0 * r1)) {
                i_result = m1 + 1;
            } else {
                i_result = m1;
            }
        }
    } else {
        if (r1 < 2.0/3.0) {
            if (r2 < (1.0 - r1)) {
                j_result = m2;
            } else {
                j_result = m2 + 1;
            }

            if ((2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)) {
                i_result = m1;
            } else {
                i_result = m1 + 1;
            }
        } else {
            if (r2 < (r1 / 2.0)) {
                i_result = m1 + 1;
                j_result = m2;
            } else {
                i_result = m1 + 1;
                j_result = m2 + 1;
            }
        }
    }

    // Fold across axes if original coordinates were negative
    if (x < 0.0) {
        if ((j_result % 2) == 0) {  // even
            long long axisi = j_result / 2;
            long long diff = i_result - axisi;
            i_result = i_result - 2 * diff;
        } else {
            long long axisi = (j_result + 1) / 2;
            long long diff = i_result - axisi;
            i_result = i_result - (2 * diff + 1);
        }
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

// ============================================================================
// Quad Edge Adjacency
// ============================================================================
//
// When a cell falls on the edge of a quad, it may belong to an adjacent quad.
// This table defines the adjacency relationships.

struct QuadAdjacency {
    bool is_upper;      // Upper hemisphere quad (1-5) vs lower (6-10)
    int up_neighbor;    // Quad above (for top edge overflow)
    int right_neighbor; // Quad to the right (for right edge overflow)
};

const QuadAdjacency kQuadAdjacency[12] = {
    {true,  0,  0},  // Quad 0: north pole (unused)
    {true,  2,  6},  // Quad 1
    {true,  3,  7},  // Quad 2
    {true,  4,  8},  // Quad 3
    {true,  5,  9},  // Quad 4
    {true,  1, 10},  // Quad 5
    {false, 2,  7},  // Quad 6
    {false, 3,  8},  // Quad 7
    {false, 4,  9},  // Quad 8
    {false, 5, 10},  // Quad 9
    {false, 1,  6},  // Quad 10
    {false, 0,  0},  // Quad 11: south pole (unused)
};

} // anonymous namespace

// ============================================================================
// Public API Implementation
// ============================================================================

void projtri_to_q2dd(int tnum, double tx, double ty,
                     int& out_quad, double& out_qx, double& out_qy) {
    if (tnum < 0 || tnum >= 20) {
        throw std::runtime_error("projtri_to_q2dd: triangle number must be 0-19");
    }

    const TriangleMapping& mapping = kTriangleMap[tnum];

    out_quad = mapping.quad;
    out_qx = tx;
    out_qy = ty;

    // Apply rotation then translation (DGGRID order)
    rotate_60deg_ccw(out_qx, out_qy, mapping.rotations);
    out_qx -= mapping.offset_x;
    out_qy -= mapping.offset_y;
}

long long get_max_ij(int aperture, int resolution) {
    if (resolution <= 0) return 0;

    double factor;
    if (aperture == 3) {
        factor = std::pow(kSqrt3, resolution);
        // Class II (odd res) uses finer substrate
        if (resolution % 2 != 0) {
            factor *= kSqrt3;
        }
    } else if (aperture == 4) {
        factor = std::pow(2.0, resolution);
    } else if (aperture == 7) {
        factor = std::pow(std::sqrt(7.0), resolution);
        factor *= std::sqrt(7.0);  // Class III substrate
    } else {
        return 0;
    }

    return static_cast<long long>(factor + 1e-9) - 1;
}

bool handle_edge_overflow(int& quad, long long& i, long long& j,
                          int aperture, int resolution) {
    long long max_coord = get_max_ij(aperture, resolution);
    long long edge_coord = max_coord + 1;

    if (i != edge_coord && j != edge_coord) {
        return false;  // Not on edge
    }

    if (quad < 1 || quad > 10) {
        return false;  // Polar quads don't overflow
    }

    const QuadAdjacency& adj = kQuadAdjacency[quad];

    if (adj.is_upper) {
        // Upper hemisphere quads (1-5)
        if (j == edge_coord) {
            if (i == 0) {
                // Top-left corner -> pole vertex
                quad = 0;  // North pole
                i = j = 0;
            } else {
                // Top edge -> up neighbor
                quad = adj.up_neighbor;
                long long new_j = edge_coord - i;
                i = 0;
                j = new_j;
            }
        } else if (i == edge_coord) {
            // Right edge -> right neighbor
            quad = adj.right_neighbor;
            i = 0;
            // j unchanged
        }
    } else {
        // Lower hemisphere quads (6-10)
        if (i == edge_coord) {
            if (j == 0) {
                // Bottom-right corner -> pole vertex
                quad = 11;  // South pole
                i = j = 0;
            } else {
                // Right edge -> right neighbor
                quad = adj.right_neighbor;
                long long new_i = edge_coord - j;
                i = new_i;
                j = 0;
            }
        } else if (j == edge_coord) {
            // Top edge -> up neighbor
            quad = adj.up_neighbor;
            j = 0;
            // i unchanged
        }
    }

    return true;
}

void q2dd_to_q2di(int quad, double qx, double qy,
                  int aperture, int resolution,
                  int& out_quad, long long& out_i, long long& out_j) {

    // Compute scale factor
    double scale;
    if (aperture == 3) {
        scale = std::pow(kSqrt3, resolution);
    } else if (aperture == 4) {
        scale = std::pow(2.0, resolution);
    } else if (aperture == 7) {
        scale = std::pow(std::sqrt(7.0), resolution);
    } else {
        throw std::runtime_error("q2dd_to_q2di: unsupported aperture");
    }

    double scaled_x = qx * scale;
    double scaled_y = qy * scale;

    // Select quantization based on grid class
    bool is_class1 = (aperture == 4) || (aperture == 3 && resolution % 2 == 0);

    if (is_class1) {
        quantize_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        quantize_class2(scaled_x, scaled_y, out_i, out_j);
    }

    out_quad = quad;
    handle_edge_overflow(out_quad, out_i, out_j, aperture, resolution);
}

void projtri_to_q2di(int tnum, double tx, double ty,
                     int aperture, int resolution,
                     int& out_quad, long long& out_i, long long& out_j) {
    int quad;
    double qx, qy;
    projtri_to_q2dd(tnum, tx, ty, quad, qx, qy);
    q2dd_to_q2di(quad, qx, qy, aperture, resolution, out_quad, out_i, out_j);
}

void q2di_to_q2dd(int quad, long long i, long long j,
                  int aperture, int resolution,
                  double& out_qx, double& out_qy) {

    bool is_class1 = (aperture == 4) || (aperture == 3 && resolution % 2 == 0);

    double x, y;
    inv_quantize_class1(i, j, x, y);

    // Compute inverse scale
    double scale;
    if (aperture == 3) {
        scale = is_class1
            ? std::pow(kSqrt3, resolution)
            : std::pow(kSqrt3, resolution + 1);  // Class II substrate
    } else if (aperture == 4) {
        scale = std::pow(2.0, resolution);
    } else if (aperture == 7) {
        scale = std::pow(std::sqrt(7.0), resolution);
    } else {
        throw std::runtime_error("q2di_to_q2dd: unsupported aperture");
    }

    out_qx = x / scale;
    out_qy = y / scale;
}

// ============================================================================
// DGGRID-compatible vertTable - Derived from First Principles
// ============================================================================
//
// This table maps (quad, subTriRegion) -> (triNum, trans, rot60)
//
// Each quad is divided into 6 regions based on the hex geometry:
//   Region 0: Upper (y > sqrt(3)*x AND y >= -sqrt(3)*x)
//   Region 1: Upper-right (y <= sqrt(3)*x AND y >= 0)
//   Region 2: Lower-right (y < 0 AND y > -sqrt(3)*x)
//   Region 3: Lower (y <= -sqrt(3)*x AND y < sqrt(3)*x)
//   Region 4: Lower-left (y >= sqrt(3)*x AND y < 0)
//   Region 5: Upper-left (y >= 0 AND y < -sqrt(3)*x)
//
// ============================================================================
// DERIVATION FROM FIRST PRINCIPLES
// ============================================================================
//
// The vertTable is the inverse of the triTable. For each (quad, region), we
// need to find which triangle contains that region and what transformation
// brings Q2DD coordinates back to PROJTRI coordinates.
//
// ICOSAHEDRON STRUCTURE:
// ---------------------
// 20 triangular faces are numbered 0-19:
//   - Faces 0-4:   North cap (around vertex 0, touching north pole)
//   - Faces 5-9:   Upper-middle band (connecting north cap to lower band)
//   - Faces 10-14: Lower-middle band (connecting upper band to south cap)
//   - Faces 15-19: South cap (around vertex 11, touching south pole)
//
// QUAD STRUCTURE:
// ---------------
// 12 quads (rhombus shapes), each containing 2 triangles:
//   - Quad 0:     North pole vertex (special - not a rhombus)
//   - Quads 1-5:  Upper hemisphere, each contains triangles (n-1, n+4) for n=1..5
//   - Quads 6-10: Lower hemisphere, each contains triangles (n+4, n+9) for n=6..10
//   - Quad 11:    South pole vertex (special - not a rhombus)
//
// TRIANGLE-TO-QUAD MAPPING (triTable, forward direction):
// -------------------------------------------------------
// From triTable in DGGRID, each triangle maps to a quad with a transformation:
//
//   Triangle | Quad | Rotation | Translation
//   ---------|------|----------|-------------
//   0        |  1   |    1     | (0, 0)        <- primary
//   1        |  2   |    1     | (0, 0)        <- primary
//   2        |  3   |    1     | (0, 0)        <- primary
//   3        |  4   |    1     | (0, 0)        <- primary
//   4        |  5   |    1     | (0, 0)        <- primary
//   5        |  1   |    4     | (-0.5, -sin60) <- secondary
//   6        |  2   |    4     | (-0.5, -sin60) <- secondary
//   7        |  3   |    4     | (-0.5, -sin60) <- secondary
//   8        |  4   |    4     | (-0.5, -sin60) <- secondary
//   9        |  5   |    4     | (-0.5, -sin60) <- secondary
//   10       |  6   |    1     | (0, 0)        <- primary
//   11       |  7   |    1     | (0, 0)        <- primary
//   12       |  8   |    1     | (0, 0)        <- primary
//   13       |  9   |    1     | (0, 0)        <- primary
//   14       | 10   |    1     | (0, 0)        <- primary
//   15       |  6   |    4     | (-0.5, -sin60) <- secondary
//   16       |  7   |    4     | (-0.5, -sin60) <- secondary
//   17       |  8   |    4     | (-0.5, -sin60) <- secondary
//   18       |  9   |    4     | (-0.5, -sin60) <- secondary
//   19       | 10   |    4     | (-0.5, -sin60) <- secondary
//
// Forward transform: rotate(rot * 60°) then subtract(trans)
// Inverse transform: add(trans) then rotate(-rot * 60°)
//
// QUAD-TO-TRIANGLE MAPPING (vertTable, inverse direction):
// --------------------------------------------------------
// For each quad, the 6 regions map to triangles based on adjacency:
//
// Upper quads (1-5) - each contains primary triangle P and secondary S:
//   Region 0: Primary triangle P (rot=-1, trans=negate of primary's)
//   Region 1: Secondary triangle S (rot=-4, trans=negate of secondary's)
//   Region 2: Lower-mid triangle (adjacent via icosahedron edge)
//   Region 3: INVALID (extends beyond icosahedron)
//   Region 4: Adjacent upper-mid secondary triangle
//   Region 5: Previous quad's primary triangle
//
// Lower quads (6-10) - similar structure but mirrored:
//   Region 0: Primary triangle (from lower-mid band)
//   Region 1: Secondary triangle (from south cap)
//   Region 2: Adjacent south cap triangle
//   Region 3: Upper quad's secondary triangle
//   Region 4: INVALID
//   Region 5: Adjacent upper-mid secondary triangle
//
// ADJACENCY DERIVATION:
// ---------------------
// From icosahedron face definition:
//   faces[20][3] = {
//     {0,1,2},{0,2,3},{0,3,4},{0,4,5},{0,5,1},     // 0-4: North cap
//     {6,2,1},{7,3,2},{8,4,3},{9,5,4},{10,1,5},    // 5-9: Upper-mid band
//     {2,6,7},{3,7,8},{4,8,9},{5,9,10},{1,10,6},   // 10-14: Lower-mid band
//     {11,7,6},{11,8,7},{11,9,8},{11,10,9},{11,6,10} // 15-19: South cap
//   }
//
// Two faces are adjacent if they share 2 vertices. For each quad region,
// the adjacent triangle is determined by which face shares the edge
// corresponding to that region's direction.
//
// For quad q (1-5):
//   - Region 0 → triangle (q-1): primary triangle of this quad
//   - Region 1 → triangle (q+4): secondary triangle of this quad
//   - Region 2 → triangle (q+9): lower-mid band (shares edge going southeast)
//   - Region 3 → INVALID (no icosahedron face in this direction)
//   - Region 4 → triangle ((q+3)%5+5): previous quad's secondary
//   - Region 5 → triangle ((q-2+5)%5): next quad's primary
//
// For quad q (6-10):
//   - Region 0 → triangle (q+4): lower-mid band primary
//   - Region 1 → triangle (q+9): south cap secondary
//   - Region 2 → triangle ((q-6+4)%5+15): adjacent south cap
//   - Region 3 → triangle (q-6+10): this quad's lower-mid adjacent
//   - Region 4 → INVALID
//   - Region 5 → triangle ((q-6+4)%5+5): upper-mid secondary
//
// TRANSFORMATION DERIVATION:
// --------------------------
// The inverse transformation parameters are computed as:
//   - rot60: Negate the forward rotation
//   - trans: The translation needed to move from Q2DD back to PROJTRI
//
// For a primary triangle (forward: rot=1, trans=(0,0)):
//   Inverse: rot=-1, but DGGRID convention stores as 1 with sign in usage
//
// For a secondary triangle (forward: rot=4, trans=(-0.5,-sin60)):
//   Inverse: rot=-4, trans is negated after rotation adjustment
//
// Cross-quad adjacencies require additional transformations based on how
// the triangles are oriented relative to each other.
//
// ============================================================================

struct VertTriVals {
    int triNum;       // Output triangle number
    double trans_x;   // Translation x (added to Q2DD before rotation)
    double trans_y;   // Translation y (added to Q2DD before rotation)
    int rot60;        // Number of 60-degree rotations (multiply by -60 for actual rotation)
    bool keep;        // Whether to keep this vertex
};

// vertTable[quad][subTri] - Derived from icosahedron geometry
//
// The derivation uses these key relationships:
//
// 1. Primary triangles of quads 1-5 are faces 0-4 (north cap)
// 2. Secondary triangles of quads 1-5 are faces 5-9 (upper-mid band)
// 3. Primary triangles of quads 6-10 are faces 10-14 (lower-mid band)
// 4. Secondary triangles of quads 6-10 are faces 15-19 (south cap)
//
// 5. Each region maps to an adjacent triangle with a specific transformation:
//    - Region 0: The "upper" direction in Q2DD space
//    - Region 1: The "upper-right" direction (60° clockwise from up)
//    - Region 2: The "lower-right" direction (120° clockwise from up)
//    - Region 3: The "lower" direction (180° from up)
//    - Region 4: The "lower-left" direction (240° clockwise from up)
//    - Region 5: The "upper-left" direction (300° clockwise from up)
//
// 6. The transformations are computed to reverse the forward triTable mapping
//    while accounting for the hexagonal geometry.
//
static const VertTriVals kVertTable[12][6] = {
    // ========================================================================
    // Quad 0 (North pole vertex)
    // ========================================================================
    // The north pole (vertex 0) is surrounded by triangles 0-4.
    // This is a special case where 5 triangles meet at a point.
    // The 6 regions map to these 5 triangles with one invalid region.
    //
    // From vertex 0, going around counter-clockwise:
    //   Triangle 0: shares edge with triangles 4 and 1
    //   Triangle 1: shares edge with triangles 0 and 2
    //   Triangle 2: shares edge with triangles 1 and 3
    //   Triangle 3: shares edge with triangles 2 and 4
    //   Triangle 4: shares edge with triangles 3 and 0
    //
    // Region assignments (empirically verified):
    //   Region 0 → Triangle 1 (rot=3)
    //   Region 1 → Triangle 0 (rot=2)
    //   Region 2 → Triangle 4 (rot=1)
    //   Region 3 → INVALID (pentagon vertex, no 6th triangle)
    //   Region 4 → Triangle 3 (rot=-1)
    //   Region 5 → Triangle 2 (rot=-2)
    {
        { 1, -0.5, -kSin60,  3, true},   // Region 0 → tri 1
        { 0, -1.0,  0.0,     2, true},   // Region 1 → tri 0
        { 4, -0.5,  kSin60,  1, true},   // Region 2 → tri 4
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 3,  1.0,  0.0,    -1, true},   // Region 4 → tri 3
        { 2,  0.5, -kSin60, -2, true}    // Region 5 → tri 2
    },

    // ========================================================================
    // Quads 1-5 (Upper hemisphere)
    // ========================================================================
    // Each quad q contains:
    //   - Primary triangle: (q-1) from north cap (faces 0-4)
    //   - Secondary triangle: (q+4) from upper-mid band (faces 5-9)
    //
    // The primary triangle transformation is: rot=1, trans=(0,0)
    // The secondary triangle transformation is: rot=4, trans=(-0.5,-sin60)
    //
    // Inverse transformations:
    //   - For primary: add (0,0), rotate -1*60° = rotate(-60°)
    //   - For secondary: add (0.5,sin60) rotated, then rotate -4*60°
    //
    // Cross-quad adjacencies (computed from icosahedron edge sharing):
    //   Region 2: Lower-mid band triangle (q+9) with special transform
    //   Region 4: Previous quad's secondary triangle
    //   Region 5: Next quad's primary triangle (wrapping around)
    //
    // Quad 1: primary=tri0, secondary=tri5
    {
        { 0,  0.0,  0.0,     1, true},   // Region 0 → tri 0 (primary)
        { 5, -0.5, -kSin60,  4, true},   // Region 1 → tri 5 (secondary)
        {14, -0.5,  kSin60,  1, true},   // Region 2 → tri 14 (lower-mid, adjacent)
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 9,  0.0,  0.0,     3, true},   // Region 4 → tri 9 (quad5's secondary)
        { 4,  1.0,  0.0,     0, true}    // Region 5 → tri 4 (quad5's primary)
    },
    // Quad 2: primary=tri1, secondary=tri6
    {
        { 1,  0.0,  0.0,     1, true},   // Region 0 → tri 1 (primary)
        { 6, -0.5, -kSin60,  4, true},   // Region 1 → tri 6 (secondary)
        {10, -0.5,  kSin60,  1, true},   // Region 2 → tri 10 (lower-mid, adjacent)
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 5,  0.0,  0.0,     3, true},   // Region 4 → tri 5 (quad1's secondary)
        { 0,  1.0,  0.0,     0, true}    // Region 5 → tri 0 (quad1's primary)
    },
    // Quad 3: primary=tri2, secondary=tri7
    {
        { 2,  0.0,  0.0,     1, true},   // Region 0 → tri 2 (primary)
        { 7, -0.5, -kSin60,  4, true},   // Region 1 → tri 7 (secondary)
        {11, -0.5,  kSin60,  1, true},   // Region 2 → tri 11 (lower-mid, adjacent)
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 6,  0.0,  0.0,     3, true},   // Region 4 → tri 6 (quad2's secondary)
        { 1,  1.0,  0.0,     0, true}    // Region 5 → tri 1 (quad2's primary)
    },
    // Quad 4: primary=tri3, secondary=tri8
    {
        { 3,  0.0,  0.0,     1, true},   // Region 0 → tri 3 (primary)
        { 8, -0.5, -kSin60,  4, true},   // Region 1 → tri 8 (secondary)
        {12, -0.5,  kSin60,  1, true},   // Region 2 → tri 12 (lower-mid, adjacent)
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 7,  0.0,  0.0,     3, true},   // Region 4 → tri 7 (quad3's secondary)
        { 2,  1.0,  0.0,     0, true}    // Region 5 → tri 2 (quad3's primary)
    },
    // Quad 5: primary=tri4, secondary=tri9
    {
        { 4,  0.0,  0.0,     1, true},   // Region 0 → tri 4 (primary)
        { 9, -0.5, -kSin60,  4, true},   // Region 1 → tri 9 (secondary)
        {13, -0.5,  kSin60,  1, true},   // Region 2 → tri 13 (lower-mid, adjacent)
        {-1, -0.5,  kSin60,  1, false},  // Region 3 → INVALID
        { 8,  0.0,  0.0,     3, true},   // Region 4 → tri 8 (quad4's secondary)
        { 3,  1.0,  0.0,     0, true}    // Region 5 → tri 3 (quad4's primary)
    },

    // ========================================================================
    // Quads 6-10 (Lower hemisphere)
    // ========================================================================
    // Each quad q contains:
    //   - Primary triangle: (q+4) from lower-mid band (faces 10-14)
    //   - Secondary triangle: (q+9) from south cap (faces 15-19)
    //
    // Lower hemisphere quads have different adjacency patterns:
    //   Region 0: Primary triangle (lower-mid band)
    //   Region 1: Secondary triangle (south cap)
    //   Region 2: Adjacent south cap triangle (wrapping)
    //   Region 3: Upper quad's lower-mid triangle (cross-hemisphere)
    //   Region 4: INVALID
    //   Region 5: Upper quad's secondary triangle
    //
    // Quad 6: primary=tri10, secondary=tri15
    {
        {10,  0.0,  0.0,     1, true},   // Region 0 → tri 10 (primary)
        {15, -0.5, -kSin60,  4, true},   // Region 1 → tri 15 (secondary)
        {19,  0.0,  0.0,    -1, true},   // Region 2 → tri 19 (adjacent south cap)
        {14, -0.5,  kSin60,  2, true},   // Region 3 → tri 14 (lower-mid, cross)
        {-1, -0.5,  kSin60,  1, false},  // Region 4 → INVALID
        { 5,  0.5, -kSin60,  4, true}    // Region 5 → tri 5 (upper secondary)
    },
    // Quad 7: primary=tri11, secondary=tri16
    {
        {11,  0.0,  0.0,     1, true},   // Region 0 → tri 11 (primary)
        {16, -0.5, -kSin60,  4, true},   // Region 1 → tri 16 (secondary)
        {15,  0.0,  0.0,    -1, true},   // Region 2 → tri 15 (adjacent south cap)
        {10, -0.5,  kSin60,  2, true},   // Region 3 → tri 10 (lower-mid, cross)
        {-1, -0.5,  kSin60,  1, false},  // Region 4 → INVALID
        { 6,  0.5, -kSin60,  4, true}    // Region 5 → tri 6 (upper secondary)
    },
    // Quad 8: primary=tri12, secondary=tri17
    {
        {12,  0.0,  0.0,     1, true},   // Region 0 → tri 12 (primary)
        {17, -0.5, -kSin60,  4, true},   // Region 1 → tri 17 (secondary)
        {16,  0.0,  0.0,    -1, true},   // Region 2 → tri 16 (adjacent south cap)
        {11, -0.5,  kSin60,  2, true},   // Region 3 → tri 11 (lower-mid, cross)
        {-1, -0.5,  kSin60,  1, false},  // Region 4 → INVALID
        { 7,  0.5, -kSin60,  4, true}    // Region 5 → tri 7 (upper secondary)
    },
    // Quad 9: primary=tri13, secondary=tri18
    {
        {13,  0.0,  0.0,     1, true},   // Region 0 → tri 13 (primary)
        {18, -0.5, -kSin60,  4, true},   // Region 1 → tri 18 (secondary)
        {17,  0.0,  0.0,    -1, true},   // Region 2 → tri 17 (adjacent south cap)
        {12, -0.5,  kSin60,  2, true},   // Region 3 → tri 12 (lower-mid, cross)
        {-1, -0.5,  kSin60,  1, false},  // Region 4 → INVALID
        { 8,  0.5, -kSin60,  4, true}    // Region 5 → tri 8 (upper secondary)
    },
    // Quad 10: primary=tri14, secondary=tri19
    {
        {14,  0.0,  0.0,     1, true},   // Region 0 → tri 14 (primary)
        {19, -0.5, -kSin60,  4, true},   // Region 1 → tri 19 (secondary)
        {18,  0.0,  0.0,    -1, true},   // Region 2 → tri 18 (adjacent south cap)
        {13, -0.5,  kSin60,  2, true},   // Region 3 → tri 13 (lower-mid, cross)
        {-1, -0.5,  kSin60,  1, false},  // Region 4 → INVALID
        { 9,  0.5, -kSin60,  4, true}    // Region 5 → tri 9 (upper secondary)
    },

    // ========================================================================
    // Quad 11 (South pole vertex)
    // ========================================================================
    // The south pole (vertex 11) is surrounded by triangles 15-19.
    // This is a special case where 5 triangles meet at a point.
    // The 6 regions map to these 5 triangles with one invalid region.
    //
    // From vertex 11, going around counter-clockwise:
    //   Triangle 15: shares edge with triangles 19 and 16
    //   Triangle 16: shares edge with triangles 15 and 17
    //   Triangle 17: shares edge with triangles 16 and 18
    //   Triangle 18: shares edge with triangles 17 and 19
    //   Triangle 19: shares edge with triangles 18 and 15
    //
    // Region assignments (empirically verified):
    //   Region 0 → Triangle 17 (rot=3)
    //   Region 1 → Triangle 18 (rot=2)
    //   Region 2 → Triangle 19 (rot=1)
    //   Region 3 → Triangle 15 (rot=0)
    //   Region 4 → INVALID (pentagon vertex, no 6th triangle)
    //   Region 5 → Triangle 16 (rot=-2)
    {
        {17, -0.5, -kSin60,  3, true},   // Region 0 → tri 17
        {18, -1.0,  0.0,     2, true},   // Region 1 → tri 18
        {19, -0.5,  kSin60,  1, true},   // Region 2 → tri 19
        {15,  0.5,  kSin60,  0, true},   // Region 3 → tri 15
        {-1,  0.0,  0.0,     0, false},  // Region 4 → INVALID
        {16,  0.5, -kSin60, -2, true}    // Region 5 → tri 16
    }
};

// Compute which of 6 sub-regions a Q2DD point falls into
// Based on DGGRID's DgQ2DDtoVertex2DDConverter::compute_subtriangle
static int compute_subtriangle(double x, double y) {
    const double tol = 1e-15;
    double xs = kSqrt3 * x;

    double xpp = xs + tol;
    double xmp = -xs + tol;
    double xpm = xs - tol;
    double xmm = -xs - tol;

    // Region 0: Upper
    if (y >= xmm && y > xpp) {
        return 0;
    }

    // Region 1: Center or upper-right
    if ((std::fabs(y) <= tol && std::fabs(x) <= tol) ||
        (y <= xpp && y >= (0.0 - tol))) {
        return 1;
    }

    // Region 2: Lower-right
    if (y < (0.0 - tol) && y > xmp) {
        return 2;
    }

    // Region 3: Lower
    if (y <= xmp && y < xpm) {
        return 3;
    }

    // Region 4: Lower-left
    if (y >= xpm && y < (0.0 - tol)) {
        return 4;
    }

    // Region 5: Upper-left
    if (y >= (0.0 - tol) && y < xmm) {
        return 5;
    }

    // Fallback - shouldn't happen for valid quad coords
    return 1;
}

// Try to convert Q2DD to PROJTRI. Returns true on success, false if the point
// is in an invalid region (e.g., outside the valid quad bounds).
bool try_q2dd_to_projtri(int quad, double qx, double qy,
                         int& out_tnum, double& out_tx, double& out_ty) {
    // Detect which of 6 sub-regions the point falls into
    int subTri = compute_subtriangle(qx, qy);

    // Look up transformation from vertTable
    const VertTriVals& triVal = kVertTable[quad][subTri];

    if (!triVal.keep || triVal.triNum < 0) {
        // This region maps to an invalid/dropped vertex
        return false;
    }

    out_tnum = triVal.triNum;

    // Apply DGGRID transformation:
    //   coord += trans
    //   coord.rotate(rot60 * -60.0)  // rotate by -60*rot60 degrees CCW
    out_tx = qx + triVal.trans_x;
    out_ty = qy + triVal.trans_y;

    // Rotate: rot60 * -60 degrees = -60 * rot60 degrees CCW
    // Which is the same as 60 * rot60 degrees CW
    // DGGRID rotates CCW for positive degrees, so rot60 * -60 means:
    //   - If rot60=4, rotate -240 degrees CCW = 120 degrees CW = -2 rotations of 60deg CCW
    rotate_60deg_ccw(out_tx, out_ty, -triVal.rot60);
    return true;
}

void q2dd_to_projtri(int quad, double qx, double qy,
                     int& out_tnum, double& out_tx, double& out_ty) {
    if (!try_q2dd_to_projtri(quad, qx, qy, out_tnum, out_tx, out_ty)) {
        throw std::runtime_error("q2dd_to_projtri: point in invalid region");
    }
}

} // namespace hexify
