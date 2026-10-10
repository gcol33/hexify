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

// The steps of Lambert's construction of Snyder's projection on a face, with
// T the face centre on the sphere. Lambert: the point S swings about T, in the
// plane of S and the normal at T, down onto the tangent plane at T, keeping
// its distance |TS| = 2 sin(z / 2) (Lambert's azimuthal equal-area
// projection). Nudge: within that plane, Snyder's adjustment moves it to the
// azimuth Az' and the radius 2 f sin(z / 2). Face: the face projection; under
// ISEA, the scaling by R' about the sphere's centre that takes the tangent
// plane onto the face plane, at distance R' from the centre.
enum class ConstructionStage { Lambert = 0, Nudge = 1, Face = 2 };

// The derivative of the map from the sphere onto the plane of `stage` at a
// point of a face: j[r][c] is the rate of change of plane coordinate r (x, y)
// along unit direction c on the sphere, c = 0 away from the face centre and
// c = 1 a quarter turn clockwise from it seen from outside (the direction of
// increasing azimuth). Plane lengths are in units of the unit sphere, so the
// face projection's plane triangle has the face's area; the singular values
// of j are Tissot's scale factors.
struct FaceScale { double j[2][2]; };

// The derivative at `geo` on `face` of the active solid, exact (forward-mode
// automatic differentiation of the projection). The nudge needs the ISEA
// projection.
FaceScale face_scale(const Geo& geo, int face,
                     ConstructionStage stage = ConstructionStage::Face);

// A point of the sphere through Lambert's construction of Snyder's projection
// onto `face` (ConstructionStage), all in radians and in units of the unit
// sphere: the arc z from the face centre and the azimuth az from the face's
// first vertex, in [0, 2 pi); the Lambert point (u, v) on the tangent plane;
// Snyder's azimuth az_prime, from the same vertex, and factor f, so that the
// nudged point is (plane_u, plane_v) / R'; Snyder's point (plane_u, plane_v)
// on the face plane, and its triangle coordinates (tx, ty), which are those
// of the forward projection. (u, v) run along the plane triangle's x and y
// axes from the face centre, and an azimuth a lies along (sin a, cos a).
// Needs the ISEA projection.
struct ConstructionPoint {
  double z, az;
  double lambert_u, lambert_v;
  double az_prime, f;
  double plane_u, plane_v;
  double tx, ty;
};

ConstructionPoint construction_point(const Geo& geo, int face);

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
