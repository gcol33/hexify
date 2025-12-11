// tri_to_q2di.cpp - Convert PROJTRI (triangle coords) to Q2DI (quad integer coords)
// Port of DGGRID's triTable and conversion logic from DgIDGGutil.cpp

#include "tri_to_q2di.h"
#include <cmath>
#include <stdexcept>

namespace {
    const double M_SIN60 = 0.86602540378443864676372317075293618L;
    const double M_SQRT3 = 1.7320508075688772935274463415058723L;
    const double M_HALF = 0.5L;
    const double M_ZERO = 0.0L;
    const double PI = 3.14159265358979323846L;
}

namespace hexify {

// triTable from DGGRID DgIDGGutil.cpp lines 178-200
// Format: { quadNum, triNum, subTri, trans_x, trans_y, rot60 }
// Note: DGGRID uses quadNum 1-10 for the middle quads (0 and 11 are polar)
static const TriTableEntry triTable[20] = {
    {  1,  0, 0,  M_ZERO,  M_ZERO,  1 },   // tri 0 -> quad 1
    {  2,  1, 0,  M_ZERO,  M_ZERO,  1 },   // tri 1 -> quad 2
    {  3,  2, 0,  M_ZERO,  M_ZERO,  1 },   // tri 2 -> quad 3
    {  4,  3, 0,  M_ZERO,  M_ZERO,  1 },   // tri 3 -> quad 4
    {  5,  4, 0,  M_ZERO,  M_ZERO,  1 },   // tri 4 -> quad 5
    {  1,  5, 1, -M_HALF, -M_SIN60, 4 },   // tri 5 -> quad 1 (subtri 1)
    {  2,  6, 1, -M_HALF, -M_SIN60, 4 },   // tri 6 -> quad 2 (subtri 1)
    {  3,  7, 1, -M_HALF, -M_SIN60, 4 },   // tri 7 -> quad 3 (subtri 1)
    {  4,  8, 1, -M_HALF, -M_SIN60, 4 },   // tri 8 -> quad 4 (subtri 1)
    {  5,  9, 1, -M_HALF, -M_SIN60, 4 },   // tri 9 -> quad 5 (subtri 1)
    {  6, 10, 0,  M_ZERO,  M_ZERO,  1 },   // tri 10 -> quad 6
    {  7, 11, 0,  M_ZERO,  M_ZERO,  1 },   // tri 11 -> quad 7
    {  8, 12, 0,  M_ZERO,  M_ZERO,  1 },   // tri 12 -> quad 8
    {  9, 13, 0,  M_ZERO,  M_ZERO,  1 },   // tri 13 -> quad 9
    { 10, 14, 0,  M_ZERO,  M_ZERO,  1 },   // tri 14 -> quad 10
    {  6, 15, 1, -M_HALF, -M_SIN60, 4 },   // tri 15 -> quad 6 (subtri 1)
    {  7, 16, 1, -M_HALF, -M_SIN60, 4 },   // tri 16 -> quad 7 (subtri 1)
    {  8, 17, 1, -M_HALF, -M_SIN60, 4 },   // tri 17 -> quad 8 (subtri 1)
    {  9, 18, 1, -M_HALF, -M_SIN60, 4 },   // tri 18 -> quad 9 (subtri 1)
    { 10, 19, 1, -M_HALF, -M_SIN60, 4 },   // tri 19 -> quad 10 (subtri 1)
};

// Rotate point by angle (in degrees)
static void rotate_point(double& x, double& y, double angle_deg) {
    double angle_rad = angle_deg * PI / 180.0;
    double cos_a = std::cos(angle_rad);
    double sin_a = std::sin(angle_rad);
    double new_x = x * cos_a - y * sin_a;
    double new_y = x * sin_a + y * cos_a;
    x = new_x;
    y = new_y;
}

void projtri_to_q2dd(int tnum, double tx, double ty,
                     int& out_quad, double& out_qx, double& out_qy) {
    if (tnum < 0 || tnum >= 20) {
        throw std::runtime_error("projtri_to_q2dd: invalid triangle number");
    }

    const TriTableEntry& entry = triTable[tnum];

    out_quad = entry.quadNum;

    // Start with input coordinates
    out_qx = tx;
    out_qy = ty;

    // Apply rotation (rot60 * 60 degrees)
    // DGGRID uses deosil (clockwise) rotations
    rotate_point(out_qx, out_qy, entry.rot60 * 60.0);

    // Subtract translation
    out_qx -= entry.trans_x;
    out_qy -= entry.trans_y;
}

// Class I hex grid quantization (from DgHexC1Grid2D::quantify)
static void quantify_class1(double x, double y, long long& out_i, long long& out_j) {
    double a1 = std::fabs(x);
    double a2 = std::fabs(y);

    double x2 = a2 / M_SIN60;
    double x1 = a1 + x2 / 2.0;

    long long m1 = static_cast<long long>(x1);
    long long m2 = static_cast<long long>(x2);

    double r1 = x1 - m1;
    double r2 = x2 - m2;

    long long i_out, j_out;

    if (r1 < 0.5) {
        if (r1 < 1.0/3.0) {
            if (r2 < (1.0 + r1) / 2.0) {
                i_out = m1; j_out = m2;
            } else {
                i_out = m1; j_out = m2 + 1;
            }
        } else {
            j_out = (r2 < (1.0 - r1)) ? m2 : m2 + 1;
            i_out = ((1.0 - r1) <= r2 && r2 < (2.0 * r1)) ? m1 + 1 : m1;
        }
    } else {
        if (r1 < 2.0/3.0) {
            j_out = (r2 < (1.0 - r1)) ? m2 : m2 + 1;
            i_out = ((2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)) ? m1 : m1 + 1;
        } else {
            if (r2 < (r1 / 2.0)) {
                i_out = m1 + 1; j_out = m2;
            } else {
                i_out = m1 + 1; j_out = m2 + 1;
            }
        }
    }

    // Fold across axes for negative coordinates
    if (x < 0.0) {
        if ((j_out % 2) == 0) {
            long long axisi = j_out / 2;
            i_out -= 2 * (i_out - axisi);
        } else {
            long long axisi = (j_out + 1) / 2;
            i_out -= 2 * (i_out - axisi) + 1;
        }
    }

    if (y < 0.0) {
        i_out -= (2 * j_out + 1) / 2;
        j_out = -j_out;
    }

    out_i = i_out;
    out_j = j_out;
}

// Class I hex center from (i,j)
static void ij_to_xy_class1(long long i, long long j, double& x, double& y) {
    x = static_cast<double>(i) - static_cast<double>(j) * 0.5;
    y = static_cast<double>(j) * M_SIN60;
}

// Class II quantization (via rotated surrogate grid)
static void quantify_class2(double x, double y, long long& out_i, long long& out_j) {
    // Rotate by -30° to surrogate grid
    const double angle = -30.0 * PI / 180.0;
    double rot_x = x * std::cos(angle) - y * std::sin(angle);
    double rot_y = x * std::sin(angle) + y * std::cos(angle);

    // Quantify in surrogate (Class I)
    long long sur_i, sur_j;
    quantify_class1(rot_x, rot_y, sur_i, sur_j);

    // Get surrogate center
    double sur_cx, sur_cy;
    ij_to_xy_class1(sur_i, sur_j, sur_cx, sur_cy);

    // Rotate back (+30°)
    double back_x = sur_cx * std::cos(-angle) - sur_cy * std::sin(-angle);
    double back_y = sur_cx * std::sin(-angle) + sur_cy * std::cos(-angle);

    // Express in substrate coordinates (√3 finer)
    double sub_x = back_x * M_SQRT3;
    double sub_y = back_y * M_SQRT3;

    // Quantify in substrate
    quantify_class1(sub_x, sub_y, out_i, out_j);
}

// Edge table from DGGRID DgIDGGBase.cpp
// Format: { quadNum, isType0, loneVert, upQuad, rightQuad }
struct QuadEdgeEntry {
    int quadNum;
    bool isType0;
    int loneVert;
    int upQuad;
    int rightQuad;
};

static const QuadEdgeEntry edgeTable[12] = {
    {  0,  true,  0,  0,  0 },  // quad 0 should never occur
    {  1,  true,  0,  2,  6 },
    {  2,  true,  0,  3,  7 },
    {  3,  true,  0,  4,  8 },
    {  4,  true,  0,  5,  9 },
    {  5,  true,  0,  1, 10 },
    {  6, false, 11,  2,  7 },
    {  7, false, 11,  3,  8 },
    {  8, false, 11,  4,  9 },
    {  9, false, 11,  5, 10 },
    { 10, false, 11,  1,  6 },
    { 11, false, 11,  0,  0 }   // quad 11 should never occur
};

long long get_max_ij(int aperture, int resolution) {
    // maxI = maxJ = maxD, where maxD is computed based on aperture and resolution
    // For Class I grids (even res for ap3, all for ap4): factor = sqrt(aperture)^res
    // For Class II grids (odd res for ap3): factor *= sqrt(3) for substrate grid
    // For Class III (ap7): factor *= sqrt(7) for substrate grid

    if (aperture == 3) {
        double factor = std::pow(M_SQRT3, resolution);
        // Class II adjustment: odd resolution means we use sqrt(3) finer substrate
        bool is_class1 = (resolution % 2 == 0);
        if (!is_class1) {
            factor *= M_SQRT3;
        }
        // Add small epsilon to prevent rounding issues (from DGGRID)
        return static_cast<long long>(factor + 1e-6) - 1;
    } else if (aperture == 4) {
        return (1LL << resolution) - 1;  // 2^res - 1
    } else if (aperture == 7) {
        double factor = std::pow(std::sqrt(7.0), resolution);
        // Class III uses sqrt(7) finer substrate
        factor *= std::sqrt(7.0);
        return static_cast<long long>(factor + 1e-6) - 1;
    }
    return 0;
}

bool handle_edge_overflow(int& quadNum, long long& i, long long& j,
                          int aperture, int resolution) {
    long long maxI = get_max_ij(aperture, resolution);
    long long maxJ = maxI;  // For hex grids, maxI == maxJ
    long long edgeI = maxI + 1;
    long long edgeJ = maxJ + 1;

    // Check if on edge
    if (i != edgeI && j != edgeJ) {
        return false;  // Not on edge
    }

    if (quadNum < 1 || quadNum > 10) {
        return false;  // Invalid quad
    }

    const QuadEdgeEntry& ec = edgeTable[quadNum];

    if (ec.isType0) {
        // Type 0 quads: 1-5 (upper hemisphere)
        if (j == edgeJ) {
            if (i == 0) {
                // Corner vertex -> lone vert (pole)
                quadNum = ec.loneVert;
                i = 0;
                j = 0;
            } else {
                // Top edge -> upQuad
                quadNum = ec.upQuad;
                long long new_i = 0;
                long long new_j = edgeI - i;  // Note: edgeI, not maxI
                i = new_i;
                j = new_j;
            }
        } else if (i == edgeI) {
            // Right edge -> rightQuad
            quadNum = ec.rightQuad;
            i = 0;
            // j stays the same
        }
    } else {
        // Type 1 quads: 6-10 (lower hemisphere)
        if (i == edgeI) {
            if (j == 0) {
                // Corner vertex -> lone vert (pole)
                quadNum = ec.loneVert;
                i = 0;
                j = 0;
            } else {
                // Right edge -> rightQuad
                quadNum = ec.rightQuad;
                long long new_i = edgeJ - j;  // Note: edgeJ, not maxJ
                long long new_j = 0;
                i = new_i;
                j = new_j;
            }
        } else if (j == edgeJ) {
            // Top edge -> upQuad
            quadNum = ec.upQuad;
            j = 0;
            // i stays the same
        }
    }

    return true;
}

void q2dd_to_q2di(int quad, double qx, double qy,
                  int aperture, int resolution,
                  int& out_quad, long long& out_i, long long& out_j) {

    // Scale coordinates based on aperture and resolution
    double scale;
    if (aperture == 3) {
        scale = std::pow(M_SQRT3, resolution);
    } else if (aperture == 4) {
        scale = std::pow(2.0, resolution);
    } else if (aperture == 7) {
        scale = std::pow(std::sqrt(7.0), resolution);
    } else {
        throw std::runtime_error("q2dd_to_q2di: unsupported aperture");
    }

    double scaled_x = qx * scale;
    double scaled_y = qy * scale;

    // Determine grid class based on aperture and resolution
    bool is_class1;
    if (aperture == 3) {
        is_class1 = (resolution % 2 == 0);
    } else if (aperture == 4) {
        is_class1 = true;  // Aperture 4 is always Class I
    } else {
        // Aperture 7 - Class III, but we use Class I quantization
        is_class1 = true;  // Simplified for now
    }

    // Quantize
    if (is_class1) {
        quantify_class1(scaled_x, scaled_y, out_i, out_j);
    } else {
        quantify_class2(scaled_x, scaled_y, out_i, out_j);
    }

    out_quad = quad;

    // Handle edge overflow - map edge cells to adjacent quads
    handle_edge_overflow(out_quad, out_i, out_j, aperture, resolution);
}

void projtri_to_q2di(int tnum, double tx, double ty,
                     int aperture, int resolution,
                     int& out_quad, long long& out_i, long long& out_j) {
    // Step 1: PROJTRI -> Q2DD
    int quad;
    double qx, qy;
    projtri_to_q2dd(tnum, tx, ty, quad, qx, qy);

    // Step 2: Q2DD -> Q2DI
    q2dd_to_q2di(quad, qx, qy, aperture, resolution, out_quad, out_i, out_j);
}

void q2di_to_q2dd(int quad, long long i, long long j,
                  int aperture, int resolution,
                  double& out_qx, double& out_qy) {

    // Determine grid class
    bool is_class1;
    if (aperture == 3) {
        is_class1 = (resolution % 2 == 0);
    } else {
        is_class1 = true;
    }

    // Get center in scaled coordinates
    double x, y;
    ij_to_xy_class1(i, j, x, y);

    // Scale back
    double scale;
    if (aperture == 3) {
        if (is_class1) {
            scale = std::pow(M_SQRT3, resolution);
        } else {
            // Class II substrate is √3 finer
            scale = std::pow(M_SQRT3, resolution + 1);
        }
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

void q2dd_to_projtri(int quad, double qx, double qy,
                     int& out_tnum, double& out_tx, double& out_ty) {
    // Find matching triTable entry for this quad
    // For now, use subTri 0 (primary triangle)
    // TODO: Determine correct subTri based on coordinates

    int tnum = -1;
    for (int i = 0; i < 20; ++i) {
        if (triTable[i].quadNum == quad && triTable[i].subTri == 0) {
            tnum = i;
            break;
        }
    }

    if (tnum < 0) {
        throw std::runtime_error("q2dd_to_projtri: no matching triangle for quad");
    }

    const TriTableEntry& entry = triTable[tnum];

    // Reverse the transformation:
    // Forward was: rotate(rot60*60), then subtract trans
    // Reverse: add trans, then rotate(-rot60*60)

    out_tx = qx + entry.trans_x;
    out_ty = qy + entry.trans_y;

    rotate_point(out_tx, out_ty, -entry.rot60 * 60.0);

    out_tnum = tnum;
}

} // namespace hexify
