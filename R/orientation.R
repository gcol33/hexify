# orientation.R
# Where an ISEA grid's icosahedron sits on the sphere.
#
# An orientation is vertex 0 of the icosahedron (vert0_lon, vert0_lat) and the
# azimuth of vertex 1 seen from vertex 0, in degrees: DGGRID's dggs_vert0_lon,
# dggs_vert0_lat and dggs_vert0_azimuth. A rotation of the sphere carries one
# orientation to another and the grid with it, so cell IDs, the hierarchy and
# neighbours read the same under every orientation; only where cells sit on the
# sphere changes.

#' Orientation of an ISEA grid's icosahedron
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
#' Snyder's equal-area ISEA projection and Fuller's projection, with the codes
#' the C++ layer reads (\code{hexify::FaceProjection}).
#' @noRd
FACE_PROJECTIONS <- c(isea = 0, fuller = 1)

#' DGGRID's names for the face projections (\code{dggs_proj})
#' @noRd
DGGS_PROJECTIONS <- toupper(names(FACE_PROJECTIONS))

#' Face projection of a grid
#'
#' A grid saved before grids carried a projection, and a legacy
#' \code{hexify_grid} list, uses ISEA. H3 grids have none.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @return \code{"isea"} or \code{"fuller"}, or \code{NA} for an H3 grid
#' @noRd
grid_projection <- function(x) {
  if (isS4(x)) {
    if (is_h3_grid(x)) return(NA_character_)
    if (.hasSlot(x, "projection") && length(x@projection) == 1L) return(x@projection)
    return("isea")
  }
  if (is.null(x$projection)) "isea" else tolower(x$projection)
}

#' The icosa argument for a projection on the default orientation
#' @param projection \code{"isea"} or \code{"fuller"}, or the choices vector
#'   of a function argument
#' @noRd
projection_icosa <- function(projection) {
  projection <- match.arg(projection, names(FACE_PROJECTIONS))
  unname(FACE_PROJECTIONS[projection])
}

#' The icosa argument of the standard ISEA orientation on the ISEA projection
#' @noRd
standard_icosa <- function() {
  c(unname(ISEA_ORIENTATION), FACE_PROJECTIONS[["isea"]])
}

#' The icosa argument the C++ layer takes
#'
#' A grid's own orientation and face projection,
#' \code{c(vert0_lon, vert0_lat, azimuth, projection)};
#' \code{numeric(0)} for no grid or an H3 grid, which the C++ layer reads as
#' the default orientation set by \code{hexify_build_icosa()} with the ISEA
#' projection.
#' @param g HexGridInfo object, legacy hexify_grid list, or NULL
#' @noRd
icosa_arg <- function(g) {
  if (is.null(g)) return(numeric(0))
  o <- grid_orientation(g)
  if (length(o) == 0L) return(numeric(0))
  c(unname(o), unname(FACE_PROJECTIONS[grid_projection(g)]))
}

#' Is this the standard ISEA orientation?
#' @noRd
is_standard_orientation <- function(o) {
  length(o) == 3L && isTRUE(all.equal(unname(o), unname(ISEA_ORIENTATION),
                                      tolerance = 0, check.attributes = FALSE))
}

#' Resolve hex_grid()'s orientation argument
#'
#' @param orientation "standard", "random", "region", "face", or a numeric
#'   \code{c(vert0_lon, vert0_lat, azimuth)}
#' @param region The area "region" and "face" centre the grid on
#' @return Named numeric \code{c(vert0_lon, vert0_lat, azimuth)}, longitude in
#'   [-180, 180) and azimuth in [0, 360)
#' @noRd
resolve_orientation <- function(orientation, region = NULL) {
  placed <- c("region", "face")
  region_misuse <- "region applies to orientation = \"region\" or \"face\""
  if (is.numeric(orientation)) {
    if (!is.null(region)) stop(region_misuse, call. = FALSE)
    return(check_orientation(orientation))
  }
  if (!is.character(orientation) || length(orientation) != 1L ||
      !orientation %in% c("standard", "random", placed)) {
    stop("orientation must be \"standard\", \"random\", \"region\", \"face\", ",
         "or c(vert0_lon, vert0_lat, azimuth) in degrees", call. = FALSE)
  }
  if (!orientation %in% placed && !is.null(region)) {
    stop(region_misuse, call. = FALSE)
  }
  if (orientation %in% placed) {
    if (is.null(region)) {
      stop("orientation = \"", orientation, "\" needs a region: c(lon, lat) ",
           "or an sf object", call. = FALSE)
    }
    centre <- region_centre(region)
  }
  switch(orientation,
    standard = ISEA_ORIENTATION,
    random = check_orientation(c(
      stats::runif(1, -180, 180),
      asin(stats::runif(1, -1, 1)) * 180 / pi,
      stats::runif(1, 0, 360)
    )),
    region = region_orientation(centre[1], centre[2]),
    face = face_orientation(centre[1], centre[2])
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

#' Orientation that places a region centre where DGGRID's REGION_CENTER does
#'
#' DGGRID places vertex 0 and reads the azimuth from two fixed points of the
#' gnomonic projection about the centre. The centre then lies at the midpoint of
#' an icosahedron edge, the middle of the two faces that share it.
#' @param lon,lat Region centre in degrees
#' @noRd
region_orientation <- function(lon, lat) {
  p0 <- gnomonic_inverse(lon, lat, DGGRID_REGION_VERT0_M / DGGRID_AUTHALIC_RADIUS_M)
  p1 <- gnomonic_inverse(lon, lat, DGGRID_REGION_AZ_POINT_M / DGGRID_AUTHALIC_RADIUS_M)
  check_orientation(c(p0[1], p0[2], gc_azimuth_deg(p0, p1)))
}

#' Orientation that places a region centre on the centre of a face
#'
#' Vertex 0 lies due north of the centre, one face circumradius away, and
#' vertex 1 is vertex 0 turned a third of a turn about the centre, so the face
#' they span with the third such vertex is centred on the region. The north
#' direction at a pole is read from the centre's longitude.
#' @param lon,lat Region centre in degrees
#' @noRd
face_orientation <- function(lon, lat) {
  r <- pi / 180
  phi <- (1 + sqrt(5)) / 2
  # Angle between a face centre and its vertices
  circum <- acos(sqrt((3 * phi + 2) / (3 * (phi + 2))))
  lo <- lon * r
  la <- lat * r
  centre <- c(cos(la) * cos(lo), cos(la) * sin(lo), sin(la))
  north <- c(-sin(la) * cos(lo), -sin(la) * sin(lo), cos(la))
  v0 <- cos(circum) * centre + sin(circum) * north
  # Rodrigues rotation of v0 by 120 degrees about the centre
  turn <- 2 * pi / 3
  cross <- c(centre[2] * v0[3] - centre[3] * v0[2],
             centre[3] * v0[1] - centre[1] * v0[3],
             centre[1] * v0[2] - centre[2] * v0[1])
  v1 <- cos(turn) * v0 + sin(turn) * cross +
    (1 - cos(turn)) * sum(centre * v0) * centre
  ll <- function(v) c(atan2(v[2], v[1]), asin(max(-1, min(1, v[3])))) / r
  p0 <- ll(v0)
  check_orientation(c(p0[1], p0[2], gc_azimuth_deg(p0, ll(v1))))
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
