#pragma once
#include <utility>

namespace hexify {

// How the ISEA inverse is solved: Snyder's projection in closed form, or by
// Newton's method on Snyder's auxiliary azimuth, kept to check the closed
// form. Fuller's inverse ignores it.
enum class InverseSolver { Closed = 0, Newton = 1 };

// Face-plane (x, y) of `face` -> (lon_deg, lat_deg), with the active solid
// and face projection.
std::pair<double,double> face_xy_to_ll(double x, double y, int face,
                                       InverseSolver solver = InverseSolver::Closed);

} // namespace hexify
