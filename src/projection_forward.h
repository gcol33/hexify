#pragma once
#include "polyhedron.h"
#include <utility>

namespace hexify {

// The projection between each spherical face and its plane triangle: Snyder's
// equal-area ISEA projection or Fuller's. Every forward and inverse face
// projection reads the active one.
enum class FaceProjection { ISEA = 0, Fuller = 1 };
void use_projection(FaceProjection p);
FaceProjection active_projection();

// Project a point onto a specific face of the solid (low-level), with the
// active face projection or with `proj`.
// Returns (icosa_triangle_x, icosa_triangle_y) face-plane coordinates
std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face);
std::pair<double,double> project_to_face(const Geo& geo, const PolyData& ico_data, int face,
                                         FaceProjection proj);

// High-level projection result
struct ProjectionResult { int face; double icosa_triangle_x; double icosa_triangle_y; };

// Forward projection: (lon, lat) -> (face, icosa_triangle_x, icosa_triangle_y)
ProjectionResult snyder_forward(double lon_deg, double lat_deg);

// Forward projection to a known face
std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg);

// Per-face azimuth offset (radians)
double snyder_get_face_azimuth_offset(int face);

} // namespace hexify
