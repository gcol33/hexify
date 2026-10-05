#pragma once
// The icosahedron orientation an R call asks for.
//
// Every entry point that reads the icosahedron takes `orient` as its first
// argument: c(vert0_lon, vert0_lat, azimuth) in degrees for a grid's own
// orientation, or an empty vector for the default one (see
// hexify::build_icosa_full()). Setting it on entry means no call reads an
// orientation a previous call left active.

#include <Rcpp.h>
#include "icosahedron.h"

inline void activate_orientation(const Rcpp::NumericVector& orient) {
  if (orient.size() == 0) {
    hexify::use_default_orientation();
    return;
  }
  if (orient.size() != 3) {
    Rcpp::stop("orient must be empty or hold vert0_lon, vert0_lat and azimuth");
  }
  hexify::Orientation o;
  o.vert0_lon_deg = orient[0];
  o.vert0_lat_deg = orient[1];
  o.azimuth_deg = orient[2];
  hexify::use_orientation(o);
}
