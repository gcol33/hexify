// cube_coordinates.h - Cube coordinate system for hexagonal grids
//
// Cube coordinates (q, r, s) satisfy the constraint q + r + s = 0.
// This provides elegant nearest-hex rounding via the "round and fix" algorithm.
//
// Reference: https://www.redblobgames.com/grids/hexagons/#rounding
//
// Copyright (c) 2024 hexify authors. MIT License.

#ifndef HEXIFY_CUBE_COORDINATES_H
#define HEXIFY_CUBE_COORDINATES_H

#include <cmath>

namespace hexify {

// ============================================================================
// Cube Coordinate Quantization
// ============================================================================
//
// Cube coordinates use three axes at 120° angles. The constraint q + r + s = 0
// means we only need two independent coordinates, but using three simplifies
// many operations.
//
// CUBE COORDINATE SYSTEM (q, r, s):
// ---------------------------------
// Standard hexagonal coordinate system where:
//   - q: "column" axis (increases to the right)
//   - r: "row" axis (increases down-right at 60° from q)
//   - s: derived axis (s = -q - r, increases down-left)
//
// These three axes point in directions 120° apart, forming a symmetric
// coordinate system for hexagonal grids. Any point can be uniquely identified
// by (q, r) since s is constrained by q + r + s = 0.
//
// To find the nearest hex center from continuous (q, r, s):
//   1. Round each coordinate independently
//   2. The sum may no longer be zero due to rounding errors
//   3. Fix by adjusting the coordinate with the largest rounding error

// Cube coordinates for hexagonal grid cells
// q, r, s are the three cube coordinate axes (q + r + s = 0)
struct CubeCoord {
    double q;  // "column" axis - increases to the right
    double r;  // "row" axis - increases down-right at 60° from q
    double s;  // derived axis - s = -q - r, increases down-left

    CubeCoord(double q_, double r_, double s_) : q(q_), r(r_), s(s_) {}

    // Round to nearest hex center, maintaining q + r + s = 0
    void round_to_nearest() {
        long long rq = std::llround(q);
        long long rr = std::llround(r);
        long long rs = std::llround(s);

        double dq = std::fabs(static_cast<double>(rq) - q);
        double dr = std::fabs(static_cast<double>(rr) - r);
        double ds = std::fabs(static_cast<double>(rs) - s);

        // Adjust the component with largest rounding error
        if (dq > dr && dq > ds) {
            rq = -rr - rs;
        } else if (dr > ds) {
            rr = -rq - rs;
        } else {
            rs = -rq - rr;
        }

        q = static_cast<double>(rq);
        r = static_cast<double>(rr);
        s = static_cast<double>(rs);
    }
};

// ============================================================================
// Coordinate Conversions
// ============================================================================

// Cartesian to cube for flat-top hex layout
// Hex width = 1 (distance between parallel edges)
inline CubeCoord cartesian_to_cube(double x, double y, double sqrt3) {
    double fj = (2.0 * y) / sqrt3;
    double fi = x + fj * 0.5;
    return CubeCoord(fi, fj, -fi - fj);
}

// Cube to Cartesian for flat-top hex layout
inline void cube_to_cartesian(double q, double r, double& x, double& y, double sin60) {
    x = q - r * 0.5;
    y = r * sin60;
}

} // namespace hexify

#endif // HEXIFY_CUBE_COORDINATES_H
