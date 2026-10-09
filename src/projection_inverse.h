#pragma once
#include <utility>

namespace hexify {

// How the inverse is solved: in closed form, or by Newton's method (on
// Snyder's auxiliary azimuth, or on Gray's equation (39) for Fuller), kept to
// check the closed form.
enum class InverseSolver { Closed = 0, Newton = 1 };

// Face-plane (x, y) of `face` -> (lon_deg, lat_deg), with the active solid
// and face projection, latitude geodetic on the active ellipsoid (see
// authalic.h).
std::pair<double,double> face_xy_to_ll(double x, double y, int face,
                                       InverseSolver solver = InverseSolver::Closed);

// The same with latitude on the sphere.
std::pair<double,double> face_xy_to_sphere_ll(double x, double y, int face,
                                              InverseSolver solver = InverseSolver::Closed);

} // namespace hexify
