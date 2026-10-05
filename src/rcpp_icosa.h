#pragma once
// The icosahedron an R call asks for: where it sits on the sphere and how its
// faces are projected onto their plane triangles.
//
// Every entry point that reads the icosahedron takes `icosa` as its first
// argument:
//   - an empty vector: the default orientation (see hexify::build_icosa_full())
//     with the ISEA projection;
//   - c(projection): the default orientation with that projection;
//   - c(vert0_lon, vert0_lat, azimuth, projection): a grid's own orientation in
//     degrees, and its projection.
// The projection is 0 for ISEA and 1 for Fuller (hexify::FaceProjection).
// Setting both on entry means no call reads a state a previous call left
// active.

#include <Rcpp.h>
#include "icosahedron.h"
#include "projection_forward.h"

inline hexify::FaceProjection icosa_projection(double code) {
  if (code == 0.0) return hexify::FaceProjection::ISEA;
  if (code == 1.0) return hexify::FaceProjection::Fuller;
  Rcpp::stop("icosa projection must be 0 (ISEA) or 1 (Fuller)");
}

inline void activate_default_icosa() {
  hexify::use_default_orientation();
  hexify::use_projection(hexify::FaceProjection::ISEA);
}

inline void activate_icosa(const Rcpp::NumericVector& icosa) {
  if (icosa.size() == 0) {
    activate_default_icosa();
    return;
  }
  if (icosa.size() == 1) {
    hexify::use_default_orientation();
    hexify::use_projection(icosa_projection(icosa[0]));
    return;
  }
  if (icosa.size() != 4) {
    Rcpp::stop("icosa must be empty, c(projection), or "
               "c(vert0_lon, vert0_lat, azimuth, projection)");
  }
  hexify::Orientation o;
  o.vert0_lon_deg = icosa[0];
  o.vert0_lat_deg = icosa[1];
  o.azimuth_deg = icosa[2];
  hexify::use_orientation(o);
  hexify::use_projection(icosa_projection(icosa[3]));
}
