#pragma once
#include "polyhedron.h"
#include <utility>

namespace hexify {

// The projection between each spherical face and its plane triangle: Snyder's
// equal-area ISEA projection, Fuller's, or van Leeuwen and Strebe's
// vertex-oriented equal-area one (IVEA). Every forward and inverse face
// projection reads the active one.
enum class FaceProjection { ISEA = 0, Fuller = 1, IVEA = 2 };
void use_projection(FaceProjection p);
FaceProjection active_projection();

// Project a point onto a specific face of the solid (low-level), with the
// active face projection or with `proj`.
// Returns (icosa_triangle_x, icosa_triangle_y) face-plane coordinates
std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face);
std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face,
                                         FaceProjection proj);

// The derivative of the face projection at a point of a face: j[r][c] is the
// rate of change of face-plane coordinate r (x, y) along unit direction c on
// the sphere, c = 0 away from the face centre and c = 1 a quarter turn
// clockwise from it seen from outside (the direction of increasing azimuth).
// Plane lengths are in units of the unit sphere, so the plane triangle has
// the face's area; the singular values of j are Tissot's scale factors.
struct FaceScale { double j[2][2]; };

// The face projection's derivative at `geo` on `face` of the active solid,
// exact (forward-mode automatic differentiation of the projection).
FaceScale face_scale(const Geo& geo, int face);

// High-level projection result
struct ProjectionResult { int face; double icosa_triangle_x; double icosa_triangle_y; };

// Forward projection: (lon, lat) -> (face, icosa_triangle_x, icosa_triangle_y),
// latitude geodetic on the active ellipsoid (see authalic.h)
ProjectionResult snyder_forward(double lon_deg, double lat_deg);

// The same from a point given in latitude on the sphere
ProjectionResult snyder_forward_sphere(double lon_deg, double lat_deg);

// Forward projection to a known face, latitude geodetic on the active
// ellipsoid
std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg);

// Per-face azimuth offset (radians)
double snyder_get_face_azimuth_offset(int face);

} // namespace hexify
