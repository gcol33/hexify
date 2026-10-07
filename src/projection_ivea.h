#pragma once
#include "polyhedron.h"
#include <utility>

namespace hexify {

// van Leeuwen and Strebe's (2006) vertex-oriented great-circle equal-area
// projection of one triangular face, which Jacovella-St-Louis et al. (2025)
// call IVEA on the icosahedron.
//
// Each of the face's six right triangles (A edge midpoint, B vertex, C centre;
// see VertexGcParams) is cut by the great circles through B, which map to
// straight lines through B': a cut at angle rho from BA meets CA at D, and its
// image meets C'A' at D' with C'D' / C'A' the share of the triangle's area
// between the cut and BC (eqs. 20-23). A point P at arc x from B on the cut
// maps to P' on B'D' with (B'P' / B'D')^2 = (1 - cos x) / (1 - cos BD)
// (eq. 28).
//
// A point is given as for Snyder's and Fuller's projections: its arc z from
// the face centre and its azimuth az from the direction of the face's first
// vertex. Face-plane coordinates are those of the ISEA projection: the face's
// plane triangle with unit edge, origin at its lower-left vertex.
//
// Past the face, each direction applies the formulas of the right triangle
// it picks the point's place in, the forward by the azimuth on the sphere and
// the inverse by the azimuth on the plane, so there the two need not undo
// each other; on the face they do.

// (z, az) -> face-plane (x, y), for T = double, or Dual2 to carry the
// derivatives along with the point
template <class T>
std::pair<T,T> ivea_face_xy(const VertexGcParams& p, T z, T az);

// Face-plane (x, y) -> (z, az) in closed form
std::pair<double,double> ivea_face_polar(const VertexGcParams& p, double x, double y);

// The same by Newton's method on the forward projection, kept to check the
// closed form
std::pair<double,double> ivea_face_polar_newton(const VertexGcParams& p, double x, double y);

} // namespace hexify
