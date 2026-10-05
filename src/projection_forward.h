#pragma once
#include "icosahedron.h"
#include <utility>

namespace hexify {

// Project a point onto a specific icosahedral face (low-level)
// Returns (icosa_triangle_x, icosa_triangle_y) face-plane coordinates
std::pair<double,double> project_to_face(const Geo& geo, const IcosaData& ico_data, int face);

// High-level projection result
struct ProjectionResult { int face; double icosa_triangle_x; double icosa_triangle_y; };

// Forward projection: (lon, lat) -> (face, icosa_triangle_x, icosa_triangle_y)
ProjectionResult snyder_forward(double lon_deg, double lat_deg);

// Forward projection to a known face
std::pair<double,double> snyder_forward_to_face(int face, double lon_deg, double lat_deg);

// Per-face azimuth offset (radians), face = 0..19
double snyder_get_face_azimuth_offset(int face);

// The derived constants the forward projection works with, for a renderer
// that runs the same projection
struct SnyderConstants { double tan_el, cos_el, cot_30, sin_g, cos_g; };
SnyderConstants snyder_constants();

} // namespace hexify
