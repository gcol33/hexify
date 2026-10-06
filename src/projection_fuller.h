#pragma once
#include <utility>

namespace hexify {

// Fuller's projection of one icosahedron face (Gray 1995, Crider 2008).
//
// A point on a face is given in polar form about the face centre: its
// great-circle distance z from the centre and its azimuth az, measured from
// the direction of the face's first vertex as the ISEA projection measures
// it. Face-plane coordinates are those of the ISEA projection: the face's
// plane triangle with unit edge, origin at its lower-left vertex.

// (z, az) -> face-plane (x, y)
std::pair<double,double> fuller_face_xy(double z, double az);

// Face-plane (x, y) -> (z, az), by Newton's method on Gray's equation (39)
std::pair<double,double> fuller_face_polar(double x, double y);

} // namespace hexify
