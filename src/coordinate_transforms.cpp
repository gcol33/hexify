// coordinate_transforms.cpp - Convert between ISEA DGGS coordinate systems
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
//   - icosa_triangle_face: Triangle index (0-19)
//   - icosa_triangle_x, icosa_triangle_y: Normalized coords within triangle [0,1]
//
// Quad XY: Quad with continuous (double) coordinates
//   - quad: Quad index (0-11, where 0=North pole, 11=South pole)
//   - quad_x, quad_y: Continuous position within quad
//
// Quad IJ: Quad with integer cell indices (used for cell ID computation)
//   - quad: Same as Quad XY
//   - i, j: Integer cell coordinates (resolution-dependent)
//
// ============================================================================
// ICOSAHEDRON GEOMETRY
// ============================================================================
//
// The icosahedron has 20 triangular faces grouped into 12 quads:
//
//           Quad 0 (North Pole)
//                  /\
//                 /  \
//           +----+----+----+----+----+
//           | Q1 | Q2 | Q3 | Q4 | Q5 |  <- Upper hemisphere (quads 1-5)
//           +----+----+----+----+----+
//           | Q6 | Q7 | Q8 | Q9 |Q10 |  <- Lower hemisphere (quads 6-10)
//           +----+----+----+----+----+
//                 \  /
//                  \/
//           Quad 11 (South Pole)
//
// Each non-polar quad contains 2 triangles forming a rhombus.
// Triangles 0-4 and 5-9 map to quads 1-5
// Triangles 10-14 and 15-19 map to quads 6-10
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
#include "constants.h"
#include <cmath>
#include <stdexcept>

namespace hexify {

namespace {

// ============================================================================
// Triangle to Quad Mapping
// ============================================================================
//
// Icosahedron face layout (20 triangles -> 12 quads):
//
// North Pole (Quad 0) at top, South Pole (Quad 11) at bottom.
// Quads 1-5: upper hemisphere, Quads 6-10: lower hemisphere.
// Each quad contains 2 triangles forming a rhombus shape.
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


// ============================================================================
// Class III Quantization (Aperture 7)
// ============================================================================
// Class III hexagons are rotated by arctan(sqrt(3/7)) ~= 19.1deg from Class I.
// This creates a grid where only 1/7 of substrate cells are valid.

// Aperture 7 rotation angle in radians
constexpr double kAp7RotRad = 19.10660535003926406149339781619697490 * kPi / 180.0;

// Class III-I (even resolutions): Class I surrogate rotated by ~19.1deg
void quantize_class3i(double x, double y, long long& out_i, long long& out_j) {
    const double c = std::cos(-kAp7RotRad);
    const double s = std::sin(-kAp7RotRad);

    // Rotate to surrogate
    double rx = x * c - y * s;
    double ry = x * s + y * c;

    // Quantize in Class I surrogate
    long long sur_i, sur_j;
    quantize_class1(rx, ry, sur_i, sur_j);

    // Get surrogate center and rotate back
    double sur_x, sur_y;
    inv_quantize_class1(sur_i, sur_j, sur_x, sur_y);

    double back_x = sur_x * c + sur_y * s;
    double back_y = -sur_x * s + sur_y * c;

    // Scale to substrate (sqrt(7)x finer) and re-quantize
    quantize_class1(back_x * kSqrt7, back_y * kSqrt7, out_i, out_j);
}

// Class III-II (odd resolutions): Class II surrogate rotated by ~19.1deg
// The surrogate is a Class II grid (pointy-top, 30° rotated from Class I).
// We need to:
//   1. Rotate to surrogate frame (-19.1 degrees)
//   2. Quantize in Class II surrogate (which is Class I at -30 degrees)
//   3. Get the Class II surrogate center
//   4. Rotate back to original frame (+19.1 degrees)
//   5. Scale to Class I substrate (sqrt(21)x) and re-quantize
void quantize_class3ii(double x, double y, long long& out_i, long long& out_j) {
    const double c_ap7 = std::cos(-kAp7RotRad);
    const double s_ap7 = std::sin(-kAp7RotRad);

    // Step 1: Rotate to surrogate frame (-19.1 degrees)
    double sur_x = x * c_ap7 - y * s_ap7;
    double sur_y = x * s_ap7 + y * c_ap7;

    // Step 2: Quantize in Class II surrogate
    // Class II = Class I rotated by -30 degrees
    constexpr double c_30 = 0.866025403784438646763723170752936183;  // cos(-30°)
    constexpr double s_30 = -0.5;  // sin(-30°)

    // Rotate to Class I orientation within the surrogate
    double c1_x = sur_x * c_30 - sur_y * s_30;
    double c1_y = sur_x * s_30 + sur_y * c_30;

    // Quantize in Class I
    long long sur1_i, sur1_j;
    quantize_class1(c1_x, c1_y, sur1_i, sur1_j);

    // Get Class I center
    double sur1_cen_x, sur1_cen_y;
    inv_quantize_class1(sur1_i, sur1_j, sur1_cen_x, sur1_cen_y);

    // Rotate back to Class II orientation (+30 degrees)
    double c2_back_x = sur1_cen_x * c_30 + sur1_cen_y * s_30;
    double c2_back_y = -sur1_cen_x * s_30 + sur1_cen_y * c_30;

    // Step 3: Rotate back to original frame (+19.1 degrees)
    double back_x = c2_back_x * c_ap7 + c2_back_y * s_ap7;
    double back_y = -c2_back_x * s_ap7 + c2_back_y * c_ap7;

    // Step 4: Scale to substrate (sqrt(21)x finer for Class III-II) and re-quantize
    quantize_class1(back_x * kSqrt21, back_y * kSqrt21, out_i, out_j);
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

void icosa_tri_to_quad_xy(int icosa_triangle_face, double icosa_triangle_x, double icosa_triangle_y,
                          int& out_quad, double& out_quad_x, double& out_quad_y) {
    if (icosa_triangle_face < 0 || icosa_triangle_face >= 20) {
        throw std::runtime_error("icosa_tri_to_quad_xy: icosa_triangle_face must be 0-19");
    }

    const TriangleMapping& mapping = kTriangleMap[icosa_triangle_face];

    out_quad = mapping.quad;
    out_quad_x = icosa_triangle_x;
    out_quad_y = icosa_triangle_y;

    // Apply rotation then translation
    rotate_60deg_ccw(out_quad_x, out_quad_y, mapping.rotations);
    out_quad_x -= mapping.offset_x;
    out_quad_y -= mapping.offset_y;
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
        // Class III-I (even res) uses sqrt(7) substrate, Class III-II (odd res) uses sqrt(21)
        bool is_class3i = (resolution % 2 == 0);
        factor *= is_class3i ? kSqrt7 : kSqrt21;
    } else {
        return 0;
    }

    return static_cast<long long>(factor + 1e-9) - 1;
}

// Handle edge overflow for upper hemisphere quads (1-5)
// Returns true if overflow was handled
inline bool handle_upper_edge(int& quad, long long& i, long long& j,
                              long long edge_coord, const QuadAdjacency& adj) {
    if (j == edge_coord) {
        // Top edge
        if (i == 0) {
            quad = 0;  // North pole
            i = j = 0;
        } else {
            quad = adj.up_neighbor;
            long long new_j = edge_coord - i;
            i = 0;
            j = new_j;
        }
        return true;
    }
    if (i == edge_coord) {
        // Right edge -> right neighbor
        quad = adj.right_neighbor;
        i = 0;
        return true;
    }
    return false;
}

// Handle edge overflow for lower hemisphere quads (6-10)
// Returns true if overflow was handled
inline bool handle_lower_edge(int& quad, long long& i, long long& j,
                              long long edge_coord, const QuadAdjacency& adj) {
    if (i == edge_coord) {
        // Right edge
        if (j == 0) {
            quad = 11;  // South pole
            i = j = 0;
        } else {
            quad = adj.right_neighbor;
            long long new_i = edge_coord - j;
            i = new_i;
            j = 0;
        }
        return true;
    }
    if (j == edge_coord) {
        // Top edge -> up neighbor
        quad = adj.up_neighbor;
        j = 0;
        return true;
    }
    return false;
}

bool handle_edge_overflow(int& quad, long long& i, long long& j,
                          int aperture, int resolution) {
    long long edge_coord = get_max_ij(aperture, resolution) + 1;

    // Quick exit: not on edge
    if (i != edge_coord && j != edge_coord) return false;

    // Polar quads don't overflow
    if (quad < 1 || quad > 10) return false;

    const QuadAdjacency& adj = kQuadAdjacency[quad];

    return adj.is_upper
        ? handle_upper_edge(quad, i, j, edge_coord, adj)
        : handle_lower_edge(quad, i, j, edge_coord, adj);
}

void quad_xy_to_ij(int quad, double quad_x, double quad_y,
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
        throw std::runtime_error("quad_xy_to_ij: unsupported aperture");
    }

    double scaled_x = quad_x * scale;
    double scaled_y = quad_y * scale;

    // Select quantization based on aperture and grid class
    if (aperture == 7) {
        // Aperture 7: Class III quantization
        bool is_class3i = (resolution % 2 == 0);
        if (is_class3i) {
            quantize_class3i(scaled_x, scaled_y, out_i, out_j);
        } else {
            quantize_class3ii(scaled_x, scaled_y, out_i, out_j);
        }
    } else if (aperture == 4 || (aperture == 3 && resolution % 2 == 0)) {
        // Class I quantization
        quantize_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        // Class II quantization (aperture 3 odd resolutions)
        quantize_class2(scaled_x, scaled_y, out_i, out_j);
    }

    out_quad = quad;
    handle_edge_overflow(out_quad, out_i, out_j, aperture, resolution);
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
    } else if (aperture == 7) {
        // Aperture 7: base scale * substrate multiplier
        double base_scale = std::pow(std::sqrt(7.0), resolution);
        bool is_class3i = (resolution % 2 == 0);
        // Class III-I substrate is sqrt(7)x finer, Class III-II is sqrt(21)x finer
        double substrate_mult = is_class3i ? kSqrt7 : kSqrt21;
        scale = base_scale * substrate_mult;
    } else {
        throw std::runtime_error("quad_ij_to_xy: unsupported aperture");
    }

    out_quad_x = x / scale;
    out_quad_y = y / scale;
}

// ============================================================================
// vertTable - Derived from First Principles
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
// brings Quad XY coordinates back to Icosa Triangle coordinates.
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
// From the triTable, each triangle maps to a quad with a transformation:
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
//   - trans: The translation needed to move from Quad XY back to Icosa Triangle
//
// For a primary triangle (forward: rot=1, trans=(0,0)):
//   Inverse: rot=-1, stored as 1 with sign applied during usage
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
    double trans_x;   // Translation x (added to Quad XY before rotation)
    double trans_y;   // Translation y (added to Quad XY before rotation)
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
//    - Region 0: The "upper" direction in Quad XY space
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
// Each region maps to a different triangle in the icosahedron.

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
bool try_quad_xy_to_icosa_tri(int quad, double quad_x, double quad_y,
                              int& out_icosa_triangle_face, double& out_icosa_triangle_x, double& out_icosa_triangle_y) {
    // Detect which of 6 sub-regions the point falls into
    int subTri = compute_subtriangle(quad_x, quad_y);

    // Look up transformation from vertTable
    const VertTriVals& triVal = kVertTable[quad][subTri];

    if (!triVal.keep || triVal.triNum < 0) {
        // This region maps to an invalid/dropped vertex
        return false;
    }

    out_icosa_triangle_face = triVal.triNum;

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
