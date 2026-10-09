#pragma once
// The solid an R call asks for: which solid, where it sits on the sphere and
// how its faces are projected onto their plane triangles.
//
// Every entry point that reads the solid takes `icosa` as its first argument:
//   - an empty vector: the default icosahedron (see hexify::build_icosa_full())
//     with the ISEA projection;
//   - c(projection) or c(projection, solid): that solid in its default
//     orientation, with that projection;
//   - c(vert0_lon, vert0_lat, azimuth, projection),
//     c(vert0_lon, vert0_lat, azimuth, projection, solid) or
//     c(vert0_lon, vert0_lat, azimuth, projection, solid, flattening): a
//     grid's own orientation in degrees, its projection, its solid and the
//     flattening of its ellipsoid.
// The projection is 0 for ISEA, 1 for Fuller, 2 for IVEA, 3 for AK and 4
// for AKW (hexify::FaceProjection), and
// the solid 0 for the icosahedron, 1 for the octahedron and 2 for the
// tetrahedron (hexify::Solid); a missing solid is the icosahedron. A missing
// flattening, or 0, is the sphere: latitudes are read on it as they are given.
// Otherwise lon/lat are geodetic on that ellipsoid (hexify::use_ellipsoid()).
// Setting all of them on entry means no call reads a state a previous call
// left active.

#include <Rcpp.h>
#include "polyhedron.h"
#include "projection_forward.h"
#include "authalic.h"
#include "hex9_warp.h"

inline hexify::FaceProjection icosa_projection(double code) {
  if (code == 0.0) return hexify::FaceProjection::ISEA;
  if (code == 1.0) return hexify::FaceProjection::Fuller;
  if (code == 2.0) return hexify::FaceProjection::IVEA;
  if (code == 3.0) return hexify::FaceProjection::AK;
  if (code == 4.0) return hexify::FaceProjection::AKW;
  Rcpp::stop("icosa projection must be 0 (ISEA), 1 (Fuller), 2 (IVEA), 3 (AK) or 4 (AKW)");
}

inline hexify::Solid icosa_solid(double code) {
  if (code == 0.0) return hexify::Solid::Icosahedron;
  if (code == 1.0) return hexify::Solid::Octahedron;
  if (code == 2.0) return hexify::Solid::Tetrahedron;
  Rcpp::stop("icosa solid must be 0 (icosahedron), 1 (octahedron) or 2 (tetrahedron)");
}

inline void activate_default_icosa() {
  hexify::use_default_orientation();
  hexify::use_projection(hexify::FaceProjection::ISEA);
  hexify::use_ellipsoid(0.0);
}

inline void activate_icosa(const Rcpp::NumericVector& icosa) {
  const R_xlen_t n = icosa.size();
  if (n == 0) {
    activate_default_icosa();
    return;
  }
  if (n != 1 && n != 2 && n != 4 && n != 5 && n != 6) {
    Rcpp::stop("icosa must be empty, c(projection), c(projection, solid), "
               "c(vert0_lon, vert0_lat, azimuth, projection), "
               "c(vert0_lon, vert0_lat, azimuth, projection, solid) or "
               "c(vert0_lon, vert0_lat, azimuth, projection, solid, flattening)");
  }
  const bool oriented = (n >= 4);
  const hexify::FaceProjection projection = icosa_projection(icosa[oriented ? 3 : 0]);
  const hexify::Solid solid = (n == 2 || n == 5 || n == 6) ? icosa_solid(icosa[oriented ? 4 : 1])
                                                           : hexify::Solid::Icosahedron;
  const double flattening = n == 6 ? icosa[5] : 0.0;
  if (!(flattening >= 0.0 && flattening < 1.0)) {
    Rcpp::stop("icosa flattening must lie in [0, 1)");
  }
  if (projection == hexify::FaceProjection::Fuller && solid != hexify::Solid::Icosahedron) {
    Rcpp::stop("Fuller's projection is defined on the icosahedron only");
  }
  const bool kaseorg = projection == hexify::FaceProjection::AK ||
                       projection == hexify::FaceProjection::AKW;
  if (kaseorg && solid != hexify::Solid::Octahedron) {
    Rcpp::stop("Kaseorg's projection is defined on the octahedron only");
  }
  if (projection == hexify::FaceProjection::AKW && !hexify::hex9::warp_ready()) {
    Rcpp::stop("the Hex9 warp field is not loaded (hex9_warp_download())");
  }
  if (oriented) {
    hexify::Orientation o;
    o.vert0_lon_deg = icosa[0];
    o.vert0_lat_deg = icosa[1];
    o.azimuth_deg = icosa[2];
    hexify::use_solid(solid, o);
  } else {
    hexify::use_default_solid(solid);
  }
  hexify::use_projection(projection);
  hexify::use_ellipsoid(flattening);
}

// Activates the solid an entry point that reads cells or quads is given, which
// must carry a hexagonal grid: its faces pair into the diamonds cell IDs are
// numbered in.
inline void activate_grid(const Rcpp::NumericVector& icosa) {
  activate_icosa(icosa);
  const hexify::SolidTopology& t = hexify::topo();
  if (!t.has_quads) {
    Rcpp::stop("the %s carries no hexagonal grid: its faces do not pair into "
               "the diamonds cell IDs are numbered in", t.name);
  }
}
