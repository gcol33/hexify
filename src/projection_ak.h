#pragma once
#include <utility>

namespace hexify {

// Kaseorg's octahedral projection ("AK"; Griffin 2026, Sec. 11a), the base of
// Hex9's own projection. A point of an octant face with weights (b_0, b_1,
// b_2) on its three vertices, which on the octahedron are orthogonal unit
// vectors V_k, goes to the direction of
//
//   sum_k t_k (t_i^2 + t_j^2 + alpha t_i^2 t_j^2)^(1/4) V_k,
//   t_k = tan(pi b_k / 2),
//
// (i, j) the other two indices and alpha = 3.2278... the coupling Kaseorg
// fitted. The plane-to-sphere direction is in closed form; the sphere to the
// plane is solved by Newton's method. Not equal-area: cell areas vary by
// about 20% across a face.
//
// As for the other face projections, a point on the sphere is its arc z from
// the face centre and its azimuth az from the face's first vertex, and
// face-plane coordinates are those of the face's plane triangle with unit
// edge, origin at its lower-left vertex. Past the face both directions extend
// the same formulas.

// Kaseorg's coupling
constexpr double kKaseorgAlpha = 3.227806237143884260376580;

// With `warped`, Hex9's own projection AKW: the face plane is the lattice
// of the trained warp (hex9_warp.h), a point P = AK^-1(p) sitting at the
// lattice point L with L + d(L) = P.

// (z, az) -> face-plane (x, y), for T = double, or Dual2 to carry the
// derivatives along with the point
template <class T>
std::pair<T,T> ak_face_xy(T z, T az, bool warped);

// Face-plane (x, y) -> (z, az): AK in closed form, after the warp's
// displacement when `warped`
std::pair<double,double> ak_face_polar(double x, double y, bool warped);

} // namespace hexify
