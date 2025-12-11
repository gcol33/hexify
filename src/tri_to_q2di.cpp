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
// DGGRID-compatible vertTable
// ============================================================================
//
// This table maps (quad, subTriRegion) -> (triNum, trans, rot60)
// Based on DGGRID's DgVertex2DDRF::vertTable_[12][6]
//
// Each quad is divided into 6 regions based on the hex geometry:
//   Region 0: Upper (y > sqrt(3)*x AND y >= -sqrt(3)*x)
//   Region 1: Upper-right (y <= sqrt(3)*x AND y >= 0)
//   Region 2: Lower-right (y < 0 AND y > -sqrt(3)*x)
//   Region 3: Lower (y <= -sqrt(3)*x AND y < sqrt(3)*x)
//   Region 4: Lower-left (y >= sqrt(3)*x AND y < 0)
//   Region 5: Upper-left (y >= 0 AND y < -sqrt(3)*x)

struct VertTriVals {
    int triNum;       // Output triangle number
    double trans_x;   // Translation x
    double trans_y;   // Translation y
    int rot60;        // Number of 60-degree rotations (multiply by -60 for actual rotation)
    bool keep;        // Whether to keep this vertex
};

// vertTable[quad][subTri] - from DGGRID's DgIDGGutil.cpp
// Note: Only quads 1-10 are used for regular cells (0 and 11 are polar)
static const VertTriVals kVertTable[12][6] = {
    // Quad 0 (North pole - special)
    {
        { 1, -0.5, -kSin60,  3, true},
        { 0, -1.0,  0.0,     2, true},
        { 4, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 3,  1.0,  0.0,    -1, true},
        { 2,  0.5, -kSin60, -2, true}
    },
    // Quad 1
    {
        { 0,  0.0,  0.0,     1, true},
        { 5, -0.5, -kSin60,  4, true},
        {14, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 9,  0.0,  0.0,     3, true},
        { 4,  1.0,  0.0,     0, true}
    },
    // Quad 2
    {
        { 1,  0.0,  0.0,     1, true},
        { 6, -0.5, -kSin60,  4, true},
        {10, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 5,  0.0,  0.0,     3, true},
        { 0,  1.0,  0.0,     0, true}
    },
    // Quad 3
    {
        { 2,  0.0,  0.0,     1, true},
        { 7, -0.5, -kSin60,  4, true},
        {11, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 6,  0.0,  0.0,     3, true},
        { 1,  1.0,  0.0,     0, true}
    },
    // Quad 4
    {
        { 3,  0.0,  0.0,     1, true},
        { 8, -0.5, -kSin60,  4, true},
        {12, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 7,  0.0,  0.0,     3, true},
        { 2,  1.0,  0.0,     0, true}
    },
    // Quad 5
    {
        { 4,  0.0,  0.0,     1, true},
        { 9, -0.5, -kSin60,  4, true},
        {13, -0.5,  kSin60,  1, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 8,  0.0,  0.0,     3, true},
        { 3,  1.0,  0.0,     0, true}
    },
    // Quad 6
    {
        {10,  0.0,  0.0,     1, true},
        {15, -0.5, -kSin60,  4, true},
        {19,  0.0,  0.0,    -1, true},
        {14, -0.5,  kSin60,  2, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 5,  0.5, -kSin60,  4, true}
    },
    // Quad 7
    {
        {11,  0.0,  0.0,     1, true},
        {16, -0.5, -kSin60,  4, true},
        {15,  0.0,  0.0,    -1, true},
        {10, -0.5,  kSin60,  2, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 6,  0.5, -kSin60,  4, true}
    },
    // Quad 8
    {
        {12,  0.0,  0.0,     1, true},
        {17, -0.5, -kSin60,  4, true},
        {16,  0.0,  0.0,    -1, true},
        {11, -0.5,  kSin60,  2, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 7,  0.5, -kSin60,  4, true}
    },
    // Quad 9
    {
        {13,  0.0,  0.0,     1, true},
        {18, -0.5, -kSin60,  4, true},
        {17,  0.0,  0.0,    -1, true},
        {12, -0.5,  kSin60,  2, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 8,  0.5, -kSin60,  4, true}
    },
    // Quad 10
    {
        {14,  0.0,  0.0,     1, true},
        {19, -0.5, -kSin60,  4, true},
        {18,  0.0,  0.0,    -1, true},
        {13, -0.5,  kSin60,  2, true},
        {-1, -0.5,  kSin60,  1, false},  // invalid
        { 9,  0.5, -kSin60,  4, true}
    },
    // Quad 11 (South pole - special)
    {
        {17, -0.5, -kSin60,  3, true},
        {18, -1.0,  0.0,     2, true},
        {19, -0.5,  kSin60,  1, true},
        {15,  0.5,  kSin60,  0, true},
        {-1,  0.0,  0.0,     0, false},  // invalid
        {16,  0.5, -kSin60, -2, true}
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
