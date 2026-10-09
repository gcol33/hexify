# orientation.R
# The solid an ISEA-family grid is built on, and where it sits on the sphere.
#
# An orientation is vertex 0 of the solid (vert0_lon, vert0_lat) and the
# azimuth of vertex 1 seen from vertex 0, in degrees: DGGRID's dggs_vert0_lon,
# dggs_vert0_lat and dggs_vert0_azimuth. A rotation of the sphere carries one
# orientation to another and the grid with it, so cell IDs, the hierarchy and
# neighbours read the same under every orientation; only where cells sit on the
# sphere changes.

#' Solids a grid or face projection can be built on
#'
#' The codes the C++ layer reads (\code{hexify::Solid}). The tetrahedron
#' carries the face projection only: its faces do not pair into the diamond
#' quads cell IDs are numbered in.
#' @noRd
POLYHEDRA <- c(icosahedron = 0, octahedron = 1, tetrahedron = 2)

#' Solids that carry a hexagonal grid
#' @noRd
GRID_POLYHEDRA <- c("icosahedron", "octahedron")

#' The standard orientation of each solid: the ISEA orientation for the
#' icosahedron; vertex 0 at the north pole and vertex 1 on the prime meridian
#' for the octahedron and the tetrahedron
#' @noRd
POLYHEDRON_ORIENTATION <- list(
  icosahedron = ISEA_ORIENTATION,
  octahedron = c(vert0_lon = 0, vert0_lat = 90, azimuth = 180),
  tetrahedron = c(vert0_lon = 0, vert0_lat = 90, azimuth = 180)
)

#' Fuller's Dymaxion orientation of the icosahedron (Sahr et al. 2003, p. 125):
#' vertex 0 at 5.2454W, 2.3009N, an adjacent vertex at azimuth 7.46658, which
#' puts all twelve vertices in the ocean
#' @noRd
DYMAXION_ORIENTATION <- c(vert0_lon = -5.2454, vert0_lat = 2.3009,
                          azimuth = 7.46658)

#' The octahedron with the poles at midpoints of its edges, as Van de Sande's
#' Gosper World places them: vertex 0 at 21.25W, 45N and vertex 1 across the
#' north pole from it. Of the longitudes in steps of 0.25 degrees this one
#' puts the six vertices furthest from the land of hexify_world, all at least
#' 619 km offshore.
#' @noRd
GOSPER_ORIENTATION <- c(vert0_lon = -21.25, vert0_lat = 45, azimuth = 0)

#' The orientation of the ISEA3H and ISEA7H DGGRS definitions registered with
#' OGC, and of DGGAL's grids: vertex 0 at 11.20E and at latitude
#' arctan(golden ratio) on the authalic sphere, azimuth 0. On WGS84 that
#' vertex lies at geodetic latitude 58.397145907431 degrees.
#' @noRd
OGC_ORIENTATION <- c(vert0_lon = 11.20, vert0_lat = ISEA_VERT0_LAT_DEG,
                     azimuth = 0)

#' Solid of a grid
#'
#' A grid saved before grids carried a solid, and a legacy \code{hexify_grid}
#' list, is built on the icosahedron. H3 grids have none.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @return \code{"icosahedron"} or \code{"octahedron"}, or \code{NA} for an
#'   H3 grid
#' @noRd
grid_polyhedron <- function(x) {
  if (isS4(x)) {
    if (is_h3_grid(x)) return(NA_character_)
    if (.hasSlot(x, "polyhedron") && length(x@polyhedron) == 1L) return(x@polyhedron)
    return("icosahedron")
  }
  if (is.null(x$polyhedron)) "icosahedron" else tolower(x$polyhedron)
}

#' What the C++ layer knows of a solid
#'
#' Its numbers of faces, vertices and diamond quads, whether it carries a
#' grid, each vertex's valence, Snyder's g and G, and the arc of an edge.
#' @param polyhedron Name of the solid
#' @noRd
solid_info <- function(polyhedron = "icosahedron") {
  cpp_solid_info(c(FACE_PROJECTIONS[["isea"]], POLYHEDRA[[polyhedron]]))
}

#' Number of diamond quads of a solid: the faces paired, 10 on the icosahedron
#' @noRd
polyhedron_diamonds <- function(polyhedron = "icosahedron") {
  solid_info(polyhedron)$n_diamonds
}

#' Orientation of an ISEA grid's solid
#'
#' A grid saved before grids carried an orientation, and a legacy
#' \code{hexify_grid} list without the DGGRID fields, has the standard ISEA
#' orientation. H3 grids have none.
#'
#' @param x HexGridInfo object or legacy hexify_grid list
#' @return Named numeric \code{c(vert0_lon, vert0_lat, azimuth)} in degrees,
#'   or \code{numeric(0)} for an H3 grid
#' @noRd
grid_orientation <- function(x) {
  if (isS4(x)) {
    if (is_h3_grid(x)) return(numeric(0))
    if (.hasSlot(x, "orientation") && length(x@orientation) == 3L) {
      return(x@orientation)
    }
    return(ISEA_ORIENTATION)
  }
  o <- c(x$pole_lon_deg, x$pole_lat_deg, x$azimuth_deg)
  if (length(o) != 3L) return(ISEA_ORIENTATION)
  stats::setNames(as.numeric(o), names(ISEA_ORIENTATION))
}

#' Face projections an ISEA-family grid can be built on
#'
#' Snyder's equal-area ISEA projection, Fuller's projection, van Leeuwen
#' and Strebe's vertex-oriented equal-area projection (IVEA) and Kaseorg's
#' octahedral projection (AK), with the codes the C++ layer reads
#' (\code{hexify::FaceProjection}).
#' @noRd
FACE_PROJECTIONS <- c(isea = 0, fuller = 1, ivea = 2, ak = 3, akw = 4)

#' Face projections that keep every cell of a resolution the same area
#' @noRd
EQUAL_AREA_PROJECTIONS <- c("isea", "ivea")

#' Face projections defined on the icosahedron only; the others take every
#' solid with triangular faces
#' @noRd
ICOSAHEDRON_PROJECTIONS <- "fuller"

#' Face projections defined on the octahedron only: Kaseorg's, plain or with
#' Hex9's warp, needs the vertices of a face to be orthogonal
#' @noRd
OCTAHEDRON_PROJECTIONS <- c("ak", "akw")

#' DGGRID's names for the face projections it has (\code{dggs_proj}); DGGRID
#' has no IVEA
#' @noRd
DGGS_PROJECTIONS <- c("ISEA", "FULLER")

#' Face projections built on Lambert's azimuthal equal-area projection about
#' each face centre (Snyder's), whose steps projection_stages() returns
#' @noRd
LAMBERT_PROJECTIONS <- "isea"

#' Whether a face projection keeps cells equal-area
#' @param projection A name in \code{FACE_PROJECTIONS}
#' @noRd
is_equal_area_projection <- function(projection) {
  projection %in% EQUAL_AREA_PROJECTIONS
}

#' Whether a face projection is Lambert's construction, adjusted
#' @param projection A name in \code{FACE_PROJECTIONS}
#' @noRd
is_lambert_projection <- function(projection) {
  projection %in% LAMBERT_PROJECTIONS
}

#' Stops when a face projection is not defined on a solid
#' @noRd
check_projection_solid <- function(projection, polyhedron) {
  if (projection %in% ICOSAHEDRON_PROJECTIONS && polyhedron != "icosahedron") {
    stop("Fuller's projection is defined on the icosahedron only", call. = FALSE)
  }
  if (projection %in% OCTAHEDRON_PROJECTIONS && polyhedron != "octahedron") {
    stop("Kaseorg's projection (\"", projection, "\") is defined on the ",
         "octahedron only: it needs the vertices of a face to be orthogonal",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Face projection of a grid
#'
#' A grid saved before grids carried a projection, and a legacy
#' \code{hexify_grid} list, uses ISEA. H3 grids have none.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @return \code{"isea"}, \code{"fuller"}, \code{"ivea"} or \code{"ak"}, or
#'   \code{NA} for an H3 grid
#' @noRd
grid_projection <- function(x) {
  if (isS4(x)) {
    if (is_h3_grid(x)) return(NA_character_)
    if (.hasSlot(x, "projection") && length(x@projection) == 1L) return(x@projection)
    return("isea")
  }
  if (is.null(x$projection)) "isea" else tolower(x$projection)
}

#' The icosa argument for a projection and solid on the solid's default
#' orientation
#' @param projection A name in \code{FACE_PROJECTIONS}, or the choices vector
#'   of a function argument
#' @param polyhedron Name of the solid, or the choices vector of a function
#'   argument
#' @noRd
projection_icosa <- function(projection, polyhedron = "icosahedron") {
  projection <- match.arg(projection, names(FACE_PROJECTIONS))
  polyhedron <- match.arg(polyhedron, names(POLYHEDRA))
  check_projection_solid(projection, polyhedron)
  if (projection == "akw") ensure_hex9_warp()
  c(unname(FACE_PROJECTIONS[projection]), unname(POLYHEDRA[polyhedron]))
}

#' The icosa argument of a solid in its standard orientation on the ISEA
#' projection
#'
#' Cell IDs, the hierarchy and neighbours do not depend on where the solid
#' sits, so code that reads only those runs in this frame.
#' @param polyhedron Name of the solid
#' @noRd
standard_icosa <- function(polyhedron = "icosahedron") {
  c(unname(POLYHEDRON_ORIENTATION[[polyhedron]]), FACE_PROJECTIONS[["isea"]],
    POLYHEDRA[[polyhedron]])
}

#' The icosa argument the C++ layer takes
#'
#' A grid's own orientation, face projection and solid,
#' \code{c(vert0_lon, vert0_lat, azimuth, projection, solid)}, followed by
#' the flattening of its ellipsoid when it reads geodetic latitude on one;
#' \code{numeric(0)} for no grid or an H3 grid, which the C++ layer reads as
#' the icosahedron in the default orientation set by
#' \code{hexify_build_icosa()} with the ISEA projection.
#' @param g HexGridInfo object, legacy hexify_grid list, or NULL
#' @noRd
icosa_arg <- function(g) {
  if (is.null(g)) return(numeric(0))
  o <- grid_orientation(g)
  if (length(o) == 0L) return(numeric(0))
  f <- grid_flattening(g)
  if (identical(grid_projection(g), "akw")) ensure_hex9_warp()
  c(unname(o), unname(FACE_PROJECTIONS[grid_projection(g)]),
    unname(POLYHEDRA[grid_polyhedron(g)]), if (f > 0) f)
}

#' An icosa argument that carries an ellipsoid of flattening f and nothing
#' else a latitude conversion reads; numeric(0) for the sphere
#' @noRd
ellipsoid_icosa <- function(flattening) {
  if (flattening > 0) c(standard_icosa("icosahedron"), flattening) else numeric(0)
}

#' Latitudes taken between geodetic and the sphere
#'
#' A grid on an ellipsoid reads geodetic latitude, and its cells live on the
#' ellipsoid's authalic sphere, where the latitude is the authalic one.
#' `sphere_lat()` gives the sphere's latitude of a geodetic one,
#' `geodetic_lat()` the reverse; on a grid without an ellipsoid both return
#' their input.
#' @param lat Latitudes in degrees
#' @param icosa The grid's icosa argument (icosa_arg())
#' @noRd
sphere_lat <- function(lat, icosa) {
  if (length(icosa) < 6L) return(lat)
  cpp_sphere_latitude(icosa, as.numeric(lat), inverse = FALSE)
}

#' @rdname sphere_lat
#' @noRd
geodetic_lat <- function(lat, icosa) {
  if (length(icosa) < 6L) return(lat)
  cpp_sphere_latitude(icosa, as.numeric(lat), inverse = TRUE)
}

#' Is this the standard orientation of the solid?
#' @noRd
is_standard_orientation <- function(o, polyhedron = "icosahedron") {
  length(o) == 3L &&
    isTRUE(all.equal(unname(o), unname(POLYHEDRON_ORIENTATION[[polyhedron]]),
                     tolerance = 0, check.attributes = FALSE))
}

#' Resolve hex_grid()'s orientation argument
#'
#' @param orientation "standard", "ogc", "dymaxion", "gosper", "random",
#'   "region", "face", or a numeric \code{c(vert0_lon, vert0_lat, azimuth)}
#' @param region The area "region" and "face" centre the grid on
#' @param polyhedron The solid the orientation places
#' @param flattening Flattening of the grid's ellipsoid, 0 for the sphere: a
#'   region's centre is read in geodetic latitude and placed on the sphere
#' @return Named numeric \code{c(vert0_lon, vert0_lat, azimuth)}, longitude in
#'   [-180, 180) and azimuth in [0, 360)
#' @noRd
resolve_orientation <- function(orientation, region = NULL,
                                polyhedron = "icosahedron", flattening = 0) {
  placed <- c("region", "face")
  region_misuse <- "region applies to orientation = \"region\" or \"face\""
  if (is.numeric(orientation)) {
    if (!is.null(region)) stop(region_misuse, call. = FALSE)
    return(check_orientation(orientation))
  }
  if (!is.character(orientation) || length(orientation) != 1L ||
      !orientation %in% c("standard", "ogc", "dymaxion", "gosper", "random", placed)) {
    stop("orientation must be \"standard\", \"ogc\", \"dymaxion\", \"gosper\", ",
         "\"random\", \"region\", \"face\", or c(vert0_lon, vert0_lat, azimuth) ",
         "in degrees", call. = FALSE)
  }
  if (orientation == "ogc" && polyhedron != "icosahedron") {
    stop("orientation = \"ogc\" places an icosahedron", call. = FALSE)
  }
  if (!orientation %in% placed && !is.null(region)) {
    stop(region_misuse, call. = FALSE)
  }
  if (orientation == "dymaxion" && polyhedron != "icosahedron") {
    stop("orientation = \"dymaxion\" places an icosahedron", call. = FALSE)
  }
  if (orientation == "gosper" && polyhedron != "octahedron") {
    stop("orientation = \"gosper\" places an octahedron", call. = FALSE)
  }
  if (orientation %in% placed) {
    if (is.null(region)) {
      stop("orientation = \"", orientation, "\" needs a region: c(lon, lat) ",
           "or an sf object", call. = FALSE)
    }
    centre <- region_centre(region)
    centre[2] <- sphere_lat(centre[2], ellipsoid_icosa(flattening))
  }
  switch(orientation,
    standard = POLYHEDRON_ORIENTATION[[polyhedron]],
    ogc = OGC_ORIENTATION,
    dymaxion = DYMAXION_ORIENTATION,
    gosper = GOSPER_ORIENTATION,
    random = check_orientation(c(
      stats::runif(1, -180, 180),
      asin(stats::runif(1, -1, 1)) * 180 / pi,
      stats::runif(1, 0, 360)
    )),
    region = region_orientation(centre[1], centre[2], polyhedron),
    face = face_orientation(centre[1], centre[2], polyhedron)
  )
}

#' Check and normalise a numeric orientation
#' @noRd
check_orientation <- function(o) {
  if (length(o) != 3L || anyNA(o) || !all(is.finite(o))) {
    stop("orientation must be three finite numbers: ",
         "c(vert0_lon, vert0_lat, azimuth) in degrees", call. = FALSE)
  }
  if (!is.null(names(o))) {
    if (!setequal(names(o), names(ISEA_ORIENTATION))) {
      stop("orientation names must be vert0_lon, vert0_lat and azimuth",
           call. = FALSE)
    }
    o <- o[names(ISEA_ORIENTATION)]
  }
  if (o[2] < -90 || o[2] > 90) {
    stop("vert0_lat must lie between -90 and 90", call. = FALSE)
  }
  stats::setNames(c(wrap_lon_deg(o[1]), unname(o[2]), unname(o[3]) %% 360),
                  names(ISEA_ORIENTATION))
}

#' Longitude in [-180, 180)
#' @noRd
wrap_lon_deg <- function(lon) {
  unname(((lon + 180) %% 360) - 180)
}

#' The centre a region orientation is placed on
#'
#' A longitude/latitude pair, or the spherical centroid of an sf, sfc, sfg or
#' bbox object's area, read in longitude/latitude.
#' @noRd
region_centre <- function(region) {
  if (is.numeric(region) && !inherits(region, "bbox")) {
    if (length(region) != 2L || !all(is.finite(region)) ||
        region[2] < -90 || region[2] > 90) {
      stop("region must be c(lon, lat) in degrees or an sf object", call. = FALSE)
    }
    return(unname(region))
  }
  geom <- if (inherits(region, "bbox")) {
    sf::st_as_sfc(region)
  } else if (inherits(region, c("sf", "sfc", "sfg"))) {
    sf::st_geometry(if (inherits(region, "sfg")) sf::st_sfc(region, crs = 4326) else region)
  } else {
    stop("region must be c(lon, lat) in degrees or an sf object", call. = FALSE)
  }
  if (is.na(sf::st_crs(geom))) sf::st_crs(geom) <- 4326
  if (!isTRUE(sf::st_is_longlat(geom))) geom <- sf::st_transform(geom, 4326)
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  centre <- sf::st_coordinates(sf::st_centroid(sf::st_union(geom)))
  unname(centre[1, 1:2])
}

#' Orientation that places a region centre on the midpoint of an edge
#'
#' On the icosahedron this is DGGRID's REGION_CENTER: DGGRID places vertex 0
#' and reads the azimuth from two fixed points of the gnomonic projection about
#' the centre. The centre then lies at the midpoint of an icosahedron edge, the
#' middle of the two faces that share it. On another solid vertex 0 lies half
#' an edge arc due north of the centre and vertex 1 the same distance due
#' south, so the centre is the midpoint of the edge they span; the north
#' direction at a pole is read from the centre's longitude.
#' @param lon,lat Region centre in degrees
#' @param polyhedron The solid the orientation places
#' @noRd
region_orientation <- function(lon, lat, polyhedron = "icosahedron") {
  if (polyhedron == "icosahedron") {
    p0 <- gnomonic_inverse(lon, lat, DGGRID_REGION_VERT0_M / DGGRID_AUTHALIC_RADIUS_M)
    p1 <- gnomonic_inverse(lon, lat, DGGRID_REGION_AZ_POINT_M / DGGRID_AUTHALIC_RADIUS_M)
    return(check_orientation(c(p0[1], p0[2], gc_azimuth_deg(p0, p1))))
  }
  r <- pi / 180
  half <- solid_info(polyhedron)$edge_arc_deg / 2 * r
  centre <- drop(unit_vec(lon, lat))
  north <- c(-sin(lat * r) * cos(lon * r), -sin(lat * r) * sin(lon * r), cos(lat * r))
  v0 <- cos(half) * centre + sin(half) * north
  v1 <- cos(half) * centre - sin(half) * north
  p <- vec_lonlat(rbind(v0, v1))
  check_orientation(c(p[1, 1], p[1, 2], gc_azimuth_deg(p[1, ], p[2, ])))
}

#' Orientation that places a region centre on the centre of a face
#'
#' Vertex 0 lies due north of the centre, one face circumradius away, and
#' vertex 1 is vertex 0 turned a third of a turn about the centre, so the face
#' they span with the third such vertex is centred on the region. The north
#' direction at a pole is read from the centre's longitude.
#' @param lon,lat Region centre in degrees
#' @param polyhedron The solid the orientation places
#' @noRd
face_orientation <- function(lon, lat, polyhedron = "icosahedron") {
  r <- pi / 180
  # Angle between a face centre and its vertices
  circum <- if (polyhedron == "icosahedron") {
    phi <- (1 + sqrt(5)) / 2
    acos(sqrt((3 * phi + 2) / (3 * (phi + 2))))
  } else {
    solid_info(polyhedron)$g_deg * r
  }
  lo <- lon * r
  la <- lat * r
  centre <- drop(unit_vec(lon, lat))
  north <- c(-sin(la) * cos(lo), -sin(la) * sin(lo), cos(la))
  v0 <- cos(circum) * centre + sin(circum) * north
  # Rodrigues rotation of v0 by 120 degrees about the centre
  turn <- 2 * pi / 3
  v1 <- cos(turn) * v0 + sin(turn) * cross3(centre, v0) +
    (1 - cos(turn)) * sum(centre * v0) * centre
  p <- vec_lonlat(rbind(v0, v1))
  check_orientation(c(p[1, 1], p[1, 2], gc_azimuth_deg(p[1, ], p[2, ])))
}

#' Inverse gnomonic projection about (lon0, lat0) of a point (x, y) on the
#' unit sphere's tangent plane, x east and y north, as DGGRID's
#' DgProjGnomonicRF::projInverse() reads it
#' @return c(lon, lat) in degrees
#' @noRd
gnomonic_inverse <- function(lon0, lat0, xy) {
  phi0 <- lat0 * pi / 180
  x <- xy[1]
  y <- xy[2]
  rh <- sqrt(x^2 + y^2)
  z <- atan(rh)
  if (abs(abs(phi0) - pi / 2) < 1e-10) {
    lat <- if (phi0 > 0) pi / 2 - z else z - pi / 2
    dlon <- if (phi0 > 0) atan2(x, -y) else atan2(x, y)
  } else {
    lat <- asin(max(-1, min(1, cos(z) * sin(phi0) + y * sin(z) * cos(phi0) / rh)))
    dlon <- atan2(x * sin(z) * cos(phi0), (cos(z) - sin(phi0) * sin(lat)) * rh)
  }
  c(wrap_lon_deg(lon0 + dlon * 180 / pi), lat * 180 / pi)
}

#' Azimuth of p2 seen from p1 along the great circle, in degrees
#' @noRd
gc_azimuth_deg <- function(p1, p2) {
  r <- pi / 180
  atan2(cos(p2[2] * r) * sin((p2[1] - p1[1]) * r),
        cos(p1[2] * r) * sin(p2[2] * r) -
          sin(p1[2] * r) * cos(p2[2] * r) * cos((p2[1] - p1[1]) * r)) / r
}
