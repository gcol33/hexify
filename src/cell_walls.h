#ifndef HEXIFY_CELL_WALLS_H
#define HEXIFY_CELL_WALLS_H

// Lengths and areas on the unit sphere along cell walls, shared by every
// backend. A wall runs from one corner to the next as an odd number of unit
// vectors: the even entries lie on the wall, and each odd entry is the wall's
// point halfway between its two neighbours in the wall's own parameter. A
// curved wall is densified before it reaches these functions until each such
// piece is nearly a circular arc; a great-circle wall carries its arcs'
// midpoints (great_circle_wall).

#include <Rcpp.h>
#include <array>
#include <cmath>
#include <vector>

namespace hexify {

using UnitVec = std::array<double, 3>;

inline UnitVec unit_from_lonlat(double lon_deg, double lat_deg) {
    const double d = M_PI / 180.0;
    const double cl = std::cos(lat_deg * d);
    return {cl * std::cos(lon_deg * d), cl * std::sin(lon_deg * d), std::sin(lat_deg * d)};
}

// Angle between two unit vectors, accurate at every separation.
inline double arc_angle(const UnitVec& a, const UnitVec& b) {
    const double cx = a[1] * b[2] - a[2] * b[1];
    const double cy = a[2] * b[0] - a[0] * b[2];
    const double cz = a[0] * b[1] - a[1] * b[0];
    return std::atan2(std::sqrt(cx * cx + cy * cy + cz * cz),
                      a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline UnitVec normalized(const UnitVec& v) {
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return {v[0] / n, v[1] / n, v[2] / n};
}

// The point a fraction t of the way along the great-circle arc a -> b.
inline UnitVec slerp(const UnitVec& a, const UnitVec& b, double t) {
    const double w = arc_angle(a, b);
    if (w < 1e-12) return a;
    const double s = std::sin(w);
    const double fa = std::sin((1.0 - t) * w) / s, fb = std::sin(t * w) / s;
    return normalized({fa * a[0] + fb * b[0], fa * a[1] + fb * b[1], fa * a[2] + fb * b[2]});
}

// A wall through the given points, joined by great-circle arcs, in the form
// above.
inline std::vector<UnitVec> great_circle_wall(const std::vector<UnitVec>& pts) {
    std::vector<UnitVec> out;
    for (size_t i = 0; i + 1 < pts.size(); i++) {
        out.push_back(pts[i]);
        out.push_back(normalized({pts[i][0] + pts[i + 1][0], pts[i][1] + pts[i + 1][1],
                                  pts[i][2] + pts[i + 1][2]}));
    }
    out.push_back(pts.back());
    return out;
}

// Length of the wall piece a -> m -> b. The two chords a m and m b fall short
// of the curve by a quarter of what the chord a b does, to leading order, so
// Richardson's extrapolation (4 (am + mb) - ab) / 3 removes that term. On a
// great circle it is am + mb exactly.
inline double piece_length(const UnitVec& a, const UnitVec& m, const UnitVec& b) {
    return (4.0 * (arc_angle(a, m) + arc_angle(m, b)) - arc_angle(a, b)) / 3.0;
}

// Signed solid angle of the spherical triangle a b c, positive when a b c
// turn counter-clockwise seen from outside (Van Oosterom and Strackee).
inline double triangle_solid_angle(const UnitVec& a, const UnitVec& b, const UnitVec& c) {
    const double triple = a[0] * (b[1] * c[2] - b[2] * c[1]) +
                          a[1] * (b[2] * c[0] - b[0] * c[2]) +
                          a[2] * (b[0] * c[1] - b[1] * c[0]);
    const double denom = 1.0 + (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) +
                         (b[0] * c[0] + b[1] * c[1] + b[2] * c[2]) +
                         (c[0] * a[0] + c[1] * a[1] + c[2] * a[2]);
    return 2.0 * std::atan2(triple, denom);
}

// Solid angle enclosed by a cell's walls, counter-clockwise, each in the form
// above, with every wall ending where the next begins. The polygon through the
// wall points is fanned from their mean, and each piece a -> m -> b adds the
// sliver between its chord and the curve: 4/3 of the triangle a m b, the area
// of a parabolic segment over its inscribed triangle (Archimedes).
inline double enclosed_solid_angle(const std::vector<std::vector<UnitVec>>& walls) {
    UnitVec o = {0.0, 0.0, 0.0};
    for (const auto& w : walls) {
        for (size_t i = 0; i + 1 < w.size(); i += 2) {
            for (int d = 0; d < 3; d++) o[d] += w[i][d];
        }
    }
    o = normalized(o);
    double omega = 0.0;
    for (const auto& w : walls) {
        for (size_t i = 0; i + 2 < w.size(); i += 2) {
            omega += triangle_solid_angle(o, w[i], w[i + 2]) +
                     4.0 / 3.0 * triangle_solid_angle(w[i], w[i + 1], w[i + 2]);
        }
    }
    return std::fabs(omega);
}

// A wall's length and the point halfway along it, both on the unit sphere.
struct WallShape {
    double length;
    UnitVec midpoint;
};

inline WallShape wall_shape(const std::vector<UnitVec>& pts) {
    double length = 0.0;
    for (size_t i = 0; i + 2 < pts.size(); i += 2) {
        length += piece_length(pts[i], pts[i + 1], pts[i + 2]);
    }
    double half = 0.5 * length;
    for (size_t i = 0; i + 2 < pts.size(); i += 2) {
        const double step = piece_length(pts[i], pts[i + 1], pts[i + 2]);
        if (half <= step) {
            // Along the two chords of the piece, in proportion
            const double am = arc_angle(pts[i], pts[i + 1]);
            const double mb = arc_angle(pts[i + 1], pts[i + 2]);
            const double t = half / step * (am + mb);
            if (t <= am) return {length, slerp(pts[i], pts[i + 1], am > 0.0 ? t / am : 0.0)};
            return {length, slerp(pts[i + 1], pts[i + 2], mb > 0.0 ? (t - am) / mb : 0.0)};
        }
        half -= step;
    }
    return {length, pts.back()};
}

// One wall between two cells, on the unit sphere: its length, the distance
// between the two cell centres, and the distance from the wall's midpoint to
// the midpoint of the arc joining the centres. Gregory et al. (2008), after
// Heikes and Randall (1995), divide the last by the first: the cell wall
// midpoint ratio.
struct WallMeasure {
    double wall;
    double centre_distance;
    double midpoint_offset;
};

inline WallMeasure measure_wall(const WallShape& w, const UnitVec& centre,
                                const UnitVec& neighbour) {
    const UnitVec mid = normalized({centre[0] + neighbour[0], centre[1] + neighbour[1],
                                    centre[2] + neighbour[2]});
    return {w.length, arc_angle(centre, neighbour), arc_angle(w.midpoint, mid)};
}

// For each neighbour of a cell, the wall it lies across: the wall whose
// midpoint is nearest the neighbour's centre. On a hexagonal tiling the wall
// facing a neighbour has its midpoint about one inradius from that
// neighbour's centre and every other wall about sqrt(3) inradii, so the
// choice holds under any projection's distortion. Returns false unless every
// wall is chosen by exactly one neighbour.
inline bool walls_of_neighbours(const std::vector<WallShape>& walls,
                                const std::vector<UnitVec>& neighbours,
                                std::vector<int>& wall_of) {
    wall_of.assign(neighbours.size(), -1);
    if (walls.size() != neighbours.size()) return false;
    std::vector<bool> taken(walls.size(), false);
    for (size_t n = 0; n < neighbours.size(); n++) {
        double best = 1e300;
        for (size_t w = 0; w < walls.size(); w++) {
            const double d = arc_angle(walls[w].midpoint, neighbours[n]);
            if (d < best) {
                best = d;
                wall_of[n] = static_cast<int>(w);
            }
        }
        if (taken[wall_of[n]]) return false;
        taken[wall_of[n]] = true;
    }
    return true;
}

// One row per measured wall: the cell (position in the input, from 1) and its
// WallMeasure. The caller names the neighbour column in its own ID type.
struct WallRows {
    std::vector<double> cell, wall, centre_distance, midpoint_offset;

    void add(R_xlen_t k, const WallMeasure& m) {
        cell.push_back(static_cast<double>(k + 1));
        wall.push_back(m.wall);
        centre_distance.push_back(m.centre_distance);
        midpoint_offset.push_back(m.midpoint_offset);
    }

    Rcpp::DataFrame frame(const char* neighbour_name, SEXP neighbour) const {
        using Rcpp::NumericVector;
        Rcpp::List cols = Rcpp::List::create(
            Rcpp::Named("cell") = NumericVector(cell.begin(), cell.end()),
            Rcpp::Named(neighbour_name) = neighbour,
            Rcpp::Named("wall") = NumericVector(wall.begin(), wall.end()),
            Rcpp::Named("centre_distance") =
                NumericVector(centre_distance.begin(), centre_distance.end()),
            Rcpp::Named("midpoint_offset") =
                NumericVector(midpoint_offset.begin(), midpoint_offset.end()));
        return Rcpp::DataFrame(cols);
    }
};

} // namespace hexify

#endif
