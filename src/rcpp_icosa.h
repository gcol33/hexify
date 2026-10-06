#pragma once
// The solid an R call asks for: which solid, where it sits on the sphere and
// how its faces are projected onto their plane triangles.
//
// Every entry point that reads the solid takes `icosa` as its first argument:
//   - an empty vector: the default icosahedron (see hexify::build_icosa_full())
//     with the ISEA projection;
//   - c(projection) or c(projection, solid): that solid in its default
//     orientation, with that projection;
//   - c(vert0_lon, vert0_lat, azimuth, projection) or
//     c(vert0_lon, vert0_lat, azimuth, projection, solid): a grid's own
//     orientation in degrees, its projection and its solid.
// The projection is 0 for ISEA and 1 for Fuller (hexify::FaceProjection), and
// the solid 0 for the icosahedron, 1 for the octahedron and 2 for the
// tetrahedron (hexify::Solid); a missing solid is the icosahedron. Setting all
// of them on entry means no call reads a state a previous call left active.

#include <Rcpp.h>
#include "polyhedron.h"
#include "projection_forward.h"

inline hexify::FaceProjection icosa_projection(double code) {
  if (code == 0.0) return hexify::FaceProjection::ISEA;
  if (code == 1.0) return hexify::FaceProjection::Fuller;
  Rcpp::stop("icosa projection must be 0 (ISEA) or 1 (Fuller)");
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
}

inline void activate_icosa(const Rcpp::NumericVector& icosa) {
  const R_xlen_t n = icosa.size();
  if (n == 0) {
    activate_default_icosa();
    return;
  }
  if (n != 1 && n != 2 && n != 4 && n != 5) {
    Rcpp::stop("icosa must be empty, c(projection), c(projection, solid), "
               "c(vert0_lon, vert0_lat, azimuth, projection) or "
               "c(vert0_lon, vert0_lat, azimuth, projection, solid)");
  }
  const bool oriented = (n >= 4);
  const hexify::FaceProjection projection = icosa_projection(icosa[oriented ? 3 : 0]);
  const hexify::Solid solid = (n == 2 || n == 5) ? icosa_solid(icosa[n - 1])
                                                 : hexify::Solid::Icosahedron;
  if (projection == hexify::FaceProjection::Fuller && solid != hexify::Solid::Icosahedron) {
    Rcpp::stop("Fuller's projection is defined on the icosahedron only");
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
