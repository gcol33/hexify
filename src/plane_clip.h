#ifndef HEXIFY_PLANE_CLIP_H
#define HEXIFY_PLANE_CLIP_H

// Convex polygons in a plane: the hull of a point set and the intersection of
// two convex polygons, both counter-clockwise. hex_aggregate() measures how a
// cell divides among the coarser cells it overlaps by clipping one against the
// other where both are straight: on a face plane (ISEA-family projections) or
// on the gnomonic plane, where great-circle arcs are straight (H3).

#include <algorithm>
#include <cstddef>
#include <vector>

namespace hexify {

struct PlanePoint {
    double x, y;
};

using PlanePolygon = std::vector<PlanePoint>;

// Twice the signed area of the triangle o a b, positive when it turns
// counter-clockwise.
inline double turn(const PlanePoint& o, const PlanePoint& a, const PlanePoint& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Convex hull, counter-clockwise, points on its edges left out (Andrew's
// monotone chain).
inline PlanePolygon convex_hull(PlanePolygon pts) {
    std::sort(pts.begin(), pts.end(), [](const PlanePoint& a, const PlanePoint& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    if (pts.size() < 3) return pts;
    PlanePolygon hull(2 * pts.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < pts.size(); i++) {
        while (k >= 2 && turn(hull[k - 2], hull[k - 1], pts[i]) <= 0.0) k--;
        hull[k++] = pts[i];
    }
    for (std::size_t i = pts.size() - 1, lower = k + 1; i-- > 0;) {
        while (k >= lower && turn(hull[k - 2], hull[k - 1], pts[i]) <= 0.0) k--;
        hull[k++] = pts[i];
    }
    hull.resize(k - 1);
    return hull;
}

// The part of a polygon inside a convex counter-clockwise polygon, by cutting
// it with the half-plane left of each clip edge in turn (Sutherland and
// Hodgman 1974). A point on a clip edge counts as inside, so edges the two
// polygons share are kept once.
inline PlanePolygon clip_convex(const PlanePolygon& subject, const PlanePolygon& clip) {
    PlanePolygon out = subject, in;
    for (std::size_t e = 0; e < clip.size() && out.size() >= 3; e++) {
        const PlanePoint& a = clip[e];
        const PlanePoint& b = clip[(e + 1) % clip.size()];
        in.swap(out);
        out.clear();
        for (std::size_t k = 0; k < in.size(); k++) {
            const PlanePoint& p = in[k];
            const PlanePoint& q = in[(k + 1) % in.size()];
            const double sp = turn(a, b, p), sq = turn(a, b, q);
            if (sp >= 0.0) out.push_back(p);
            if ((sp >= 0.0) != (sq >= 0.0)) {
                const double t = sp / (sp - sq);
                out.push_back({p.x + t * (q.x - p.x), p.y + t * (q.y - p.y)});
            }
        }
    }
    if (out.size() < 3) out.clear();
    return out;
}

} // namespace hexify

#endif
