# R/constants.R - Package-wide constants and helpers
#
# Centralizes magic numbers to avoid duplication and improve maintainability.

# =============================================================================
# Body Geometry Constants
# =============================================================================

#' Total Earth surface area in square kilometers ('WGS84' ellipsoid)
#' @noRd
EARTH_SURFACE_KM2 <- 510065621.724078904704516

#' Mean Earth radius in kilometers ('WGS84' mean radius)
#' Calculated as (2*a + b) / 3 where a = 6378.137 km (equatorial) and b = 6356.7523142 km (polar)
#' @noRd
EARTH_RADIUS_KM <- 6371.0088

#' How far an ISEA cell edge drawn in lon/lat may stray from the true edge
#'
#' A cell edge is straight in the projection plane and curved in lon/lat.
#' Polygon builders split each edge until every piece follows the edge to
#' within this fraction of its length both as a straight lon/lat chord and as
#' a great-circle arc, the edge s2 reads, which keeps a cell's drawn area
#' within about 5e-4 of the cell that `lonlat_to_cell()` assigns, read either way.
#' @noRd
CELL_EDGE_TOLERANCE <- 1e-3

#' How closely cell walls are followed on the sphere to measure them
#'
#' Each wall is halved until the true wall's midpoint between two consecutive
#' points lies within this fraction of their distance of the great-circle
#' plane through them; lengths and areas are then corrected for the curve
#' between them. Cell areas and perimeters come out within about 1e-8 of their
#' converged values, and areas of a whole grid add up to the body's to 1e-14.
#' @noRd
CELL_WALL_TOLERANCE <- 1e-4

#' How far a corner is moved towards its cell's centre to find the coarser
#' cells it lies in
#'
#' As a fraction of the corner's distance from the centre. A corner on the
#' boundary of coarser cells then lies inside one of them, and the parts of a
#' cell in each coarser cell, a third, a half or a twelfth of it, are far
#' wider than the step.
#' @noRd
CORNER_STEP <- 1e-3

#' Longest edge of the triangles hex_globe() draws surfaces with
#'
#' In triangle coordinates, where a face edge is 1. A triangle of the sphere
#' this wide sags below the sphere by under 3e-4 radii, less than the lift
#' the widget gives each layer above the one below.
#' @noRd
GLOBE_MESH_SPACING <- 0.04

#' Fractions of the radius each globe layer is lifted off the surface
#'
#' A later layer lies on top: each step is more than a surface triangle's
#' sag below the sphere at GLOBE_MESH_SPACING. The widget and
#' hexglobe::render_scene() both read these.
#' @noRd
GLOBE_LIFT <- list(ocean = 0, land = 6e-4, cells = 1.2e-3, grid = 2e-3,
                   coast = 2.2e-3, edges = 2.4e-3)

#' Approximate km per degree of latitude (at equator)
#' @noRd
KM_PER_DEGREE <- 111.0

#' Mean radii of solar system bodies, in kilometers
#'
#' 'IAU' mean radii (Archinal et al. 2018, Report of the 'IAU' Working Group on
#' Cartographic Coordinates and Rotational Elements: 2015) as tabulated by 'JPL'
#' Solar System Dynamics. "earth" carries the 'WGS84' mean radius instead, so a
#' grid built by name matches one built from the package default.
#' @noRd
BODY_RADII_KM <- c(
  mercury   = 2439.4,
  venus     = 6051.8,
  earth     = EARTH_RADIUS_KM,
  moon      = 1737.4,
  mars      = 3389.50,
  ceres     = 469.7,
  jupiter   = 69911,
  io        = 1821.49,
  europa    = 1560.80,
  ganymede  = 2631.20,
  callisto  = 2410.30,
  saturn    = 58232,
  enceladus = 252.10,
  titan     = 2574.76,
  uranus    = 25362,
  neptune   = 24622,
  pluto     = 1188.3
)

#' Reference ellipsoids of revolution: semi-major axis a in km and flattening f
#'
#' 'WGS84' and 'GRS80' as 'PROJ' lists them (\code{sf::sf_proj_info("ellps")}:
#' a = 6378137 m, 1/f = 298.257223563 and 298.257222101). "earth" is 'WGS84'.
#' The planets are the 'IAU' size and shape parameters of Archinal et al.
#' (2018), Table 4, equatorial and polar radius: Mars with the average polar
#' radius, which the report recommends for a best-fitting ellipsoid; Jupiter,
#' Saturn, Uranus and Neptune at their one-bar surface. The report's
#' recommended shapes for Mercury, Venus and the Moon are spheres, and its
#' satellites with a measured shape are triaxial.
#' @noRd
ELLIPSOIDS <- list(
  wgs84   = c(a_km = 6378.137, f = 1 / 298.257223563),
  grs80   = c(a_km = 6378.137, f = 1 / 298.257222101),
  mars    = c(a_km = 3396.19, f = 1 - 3376.20 / 3396.19),
  jupiter = c(a_km = 71492, f = 1 - 66854 / 71492),
  saturn  = c(a_km = 60268, f = 1 - 54364 / 60268),
  uranus  = c(a_km = 25559, f = 1 - 24973 / 25559),
  neptune = c(a_km = 24764, f = 1 - 24341 / 24764)
)
ELLIPSOIDS$earth <- ELLIPSOIDS$wgs84

#' Names an ellipsoid prints with
#' @noRd
ELLIPSOID_LABELS <- c(wgs84 = "WGS84", grs80 = "GRS80", mars = "Mars",
                      jupiter = "Jupiter", saturn = "Saturn",
                      uranus = "Uranus", neptune = "Neptune")

# =============================================================================
# Body Geometry Helpers
# =============================================================================

#' Resolve hex_grid()'s ellipsoid argument
#'
#' @param ellipsoid NULL, a name in ELLIPSOIDS, or c(a, f): the semi-major
#'   axis in km and the flattening
#' @return Named numeric c(a_km, f), or numeric(0) for NULL
#' @noRd
resolve_ellipsoid <- function(ellipsoid) {
  if (is.null(ellipsoid)) return(numeric(0))
  if (is.character(ellipsoid)) {
    if (length(ellipsoid) != 1L || is.na(ellipsoid)) {
      stop("ellipsoid must be a single name or c(a, f)", call. = FALSE)
    }
    key <- tolower(trimws(ellipsoid))
    if (!key %in% names(ELLIPSOIDS)) {
      stop(sprintf(paste0(
        "Unknown ellipsoid \"%s\". Named ellipsoids are: %s. Any other ",
        "ellipsoid takes c(a, f): its semi-major axis in km and its flattening."),
        ellipsoid, paste(names(ELLIPSOIDS), collapse = ", ")), call. = FALSE)
    }
    return(ELLIPSOIDS[[key]])
  }
  if (!is.numeric(ellipsoid) || length(ellipsoid) != 2L || anyNA(ellipsoid) ||
      !all(is.finite(ellipsoid)) || ellipsoid[1] <= 0 ||
      ellipsoid[2] < 0 || ellipsoid[2] >= 1) {
    stop("ellipsoid must be a name, or c(a, f) with a > 0 km and 0 <= f < 1",
         call. = FALSE)
  }
  c(a_km = unname(ellipsoid[1]), f = unname(ellipsoid[2]))
}

#' Ellipsoid of a grid
#'
#' A grid saved before grids carried an ellipsoid, an H3 grid and a grid built
#' without one read latitude on the sphere.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @return Named numeric c(a_km, f), or numeric(0)
#' @noRd
grid_ellipsoid <- function(x) {
  e <- if (isS4(x)) {
    if (.hasSlot(x, "ellipsoid")) x@ellipsoid else NULL
  } else {
    x$ellipsoid
  }
  if (length(e) != 2L) return(numeric(0))
  stats::setNames(as.numeric(e), c("a_km", "f"))
}

#' Flattening of a grid's ellipsoid, 0 for the sphere
#' @noRd
grid_flattening <- function(x) {
  e <- grid_ellipsoid(x)
  if (length(e) == 2L) e[["f"]] else 0
}

#' Radius of the sphere of an ellipsoid's area, in km
#' @param e Named numeric c(a_km, f)
#' @noRd
authalic_radius_km <- function(e) {
  e[["a_km"]] * cpp_authalic_radius_ratio(e[["f"]])
}

#' Name of a grid's ellipsoid, or its axes
#' @noRd
ellipsoid_label <- function(e) {
  for (key in names(ELLIPSOID_LABELS)) {
    if (identical(unname(e), unname(ELLIPSOIDS[[key]]))) return(ELLIPSOID_LABELS[[key]])
  }
  sprintf("a = %s km, f = %s", format(e[["a_km"]], digits = 10),
          format(e[["f"]], digits = 10))
}

#' Is this one of Earth's reference ellipsoids?
#' @noRd
is_earth_ellipsoid <- function(e) {
  length(e) == 2L && any(vapply(ELLIPSOIDS[c("wgs84", "grs80")], function(x)
    identical(unname(e), unname(x)), logical(1)))
}

#' Longlat CRS on an ellipsoid, as a 'PROJ' string
#'
#' 'WGS84' is EPSG:4326; any other ellipsoid has no EPSG code to name it, so
#' it carries a longlat CRS on its own axes.
#' @param e Named numeric c(a_km, f)
#' @noRd
ellipsoid_crs <- function(e) {
  if (identical(unname(e), unname(ELLIPSOIDS$wgs84))) return(4326L)
  if (identical(unname(e), unname(ELLIPSOIDS$grs80))) {
    return("+proj=longlat +ellps=GRS80 +no_defs")
  }
  metres <- function(x) sub("[.]?0+$", "", sprintf("%.4f", x * 1000))
  if (e[["f"]] == 0) return(body_crs_string(e[["a_km"]]))
  sprintf("+proj=longlat +a=%s +b=%s +no_defs", metres(e[["a_km"]]),
          metres(e[["a_km"]] * (1 - e[["f"]])))
}

#' Radius of a body in km, from a number or a name
#'
#' @param radius_km Positive number, or a name from BODY_RADII_KM
#' @return Radius in kilometers
#' @noRd
resolve_radius_km <- function(radius_km) {
  if (is.character(radius_km)) {
    if (length(radius_km) != 1L || is.na(radius_km)) {
      stop("radius_km must be a single body name or a positive number of kilometers")
    }
    key <- tolower(trimws(radius_km))
    if (!key %in% names(BODY_RADII_KM)) {
      stop(sprintf(
        "Unknown body \"%s\". Named bodies are: %s. Any other body takes its mean radius in km.",
        radius_km, paste(names(BODY_RADII_KM), collapse = ", ")
      ))
    }
    return(unname(BODY_RADII_KM[[key]]))
  }

  if (!is.numeric(radius_km) || length(radius_km) != 1L || is.na(radius_km) ||
      !is.finite(radius_km) || radius_km <= 0) {
    stop("radius_km must be a single positive number of kilometers, or a body name such as \"mars\"")
  }

  as.numeric(radius_km)
}

#' Surface area of a body in km^2
#'
#' The area of a sphere of the given radius. Earth's radius returns the 'WGS84'
#' ellipsoid area, which is what an Earth grid is sized against; an ellipsoid
#' has a slightly different area from the sphere of its mean radius, so the two
#' part in the seventh significant figure.
#' @param radius_km Radius in kilometers
#' @noRd
body_surface_km2 <- function(radius_km) {
  if (radius_km == EARTH_RADIUS_KM) return(EARTH_SURFACE_KM2)
  4 * pi * radius_km^2
}

#' Kilometers per degree of arc at a body's radius
#' @param radius_km Radius in kilometers
#' @noRd
km_per_degree <- function(radius_km) {
  KM_PER_DEGREE * radius_km / EARTH_RADIUS_KM
}

#' Centre spacing of a regular hexagon of a given area
#'
#' A regular hexagon of area A has flat-to-flat width sqrt(2 A / sqrt(3)),
#' which is also the distance between the centres of neighbouring cells.
#' @param area_km2 Hexagon area in km^2
#' @noRd
hex_spacing_km <- function(area_km2) {
  sqrt(2 * area_km2 / sqrt(3))
}

#' Area of a regular hexagon of a given centre spacing (inverse of
#' hex_spacing_km())
#' @param spacing_km Flat-to-flat width in km
#' @noRd
hex_area_from_spacing <- function(spacing_km) {
  spacing_km^2 * sqrt(3) / 2
}

#' Characteristic length scale: the diameter of a circle of a cell's area
#' @param area_km2 Cell area in km^2
#' @noRd
cls_km <- function(area_km2) {
  2 * sqrt(area_km2 / pi)
}

#' Number of cells of an H3 grid at a resolution
#' @param resolution H3 resolution
#' @noRd
h3_n_cells <- function(resolution) {
  2 + 120 * 7^resolution
}

#' Mean cell area of an ISEA grid in km^2: the body's area over the cell count
#' @param aperture Aperture spelling (pure or mixed)
#' @param resolution Resolution
#' @param radius_km Radius of the body in km
#' @param polyhedron The solid the grid is built on
#' @noRd
mean_cell_area_km2 <- function(aperture, resolution, radius_km = EARTH_RADIUS_KM,
                               polyhedron = "icosahedron") {
  body_surface_km2(radius_km) / aperture_n_cells(aperture, resolution, polyhedron)
}

#' Radius a grid is sized against, in km
#'
#' A grid carrying neither the slot nor the field is an Earth grid.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @noRd
grid_radius_km <- function(x) {
  r <- if (isS4(x)) {
    if (.hasSlot(x, "radius_km")) x@radius_km else EARTH_RADIUS_KM
  } else {
    x$radius_km
  }
  if (is.null(r) || length(r) != 1L || is.na(r)) EARTH_RADIUS_KM else as.numeric(r)
}

#' Is this grid sized against Earth?
#' @param x HexGridInfo object or legacy hexify_grid list
#' @noRd
is_earth_grid <- function(x) {
  grid_radius_km(x) == EARTH_RADIUS_KM || is_earth_ellipsoid(grid_ellipsoid(x))
}

#' Do two grids cover the same body?
#'
#' Earth grids do, whether sized on the sphere or on an Earth ellipsoid; any
#' other two when they share a radius.
#' @noRd
same_body_grids <- function(a, b) {
  (is_earth_grid(a) && is_earth_grid(b)) || grid_radius_km(a) == grid_radius_km(b)
}

# =============================================================================
# Body CRS Helpers
# =============================================================================

#' Longlat CRS on the sphere of a body's radius, as a 'PROJ' string
#'
#' 'EPSG' codes name Earth reference systems, so a grid on another body carries
#' a longlat CRS on the sphere of its own radius instead.
#' @param radius_km Radius in kilometers
#' @noRd
body_crs_string <- function(radius_km) {
  metres <- sub("[.]?0+$", "", sprintf("%.4f", radius_km * 1000))
  sprintf("+proj=longlat +R=%s +no_defs", metres)
}

#' A CRS specification read by sf, or NA_crs_ if sf cannot read it
#' @param crs An 'EPSG' code or a 'PROJ' or 'WKT' string
#' @noRd
parse_crs <- function(crs) {
  tryCatch(sf::st_crs(crs), error = function(e) sf::NA_crs_)
}

#' The CRS a new grid stores
#'
#' NULL takes 'WGS84' on Earth, a longlat CRS on the grid's ellipsoid where it
#' has one, and the body's own sphere elsewhere. An 'EPSG' code stores as an
#' integer, any other CRS as the string sf reads it from.
#' @param crs NULL, an 'EPSG' code, or a 'PROJ' or 'WKT' string
#' @param radius_km Radius the grid is sized against, in kilometers
#' @param ellipsoid The grid's ellipsoid, c(a_km, f), or numeric(0)
#' @noRd
resolve_crs <- function(crs, radius_km, ellipsoid = numeric(0)) {
  if (is.null(crs)) {
    if (length(ellipsoid) == 2L) return(ellipsoid_crs(ellipsoid))
    if (radius_km == EARTH_RADIUS_KM) return(4326L)
    return(body_crs_string(radius_km))
  }

  if (is.character(crs)) {
    if (length(crs) != 1L || is.na(crs) || !nzchar(crs)) {
      stop("crs must be a single non-empty CRS string, or an EPSG code")
    }
    if (is.na(parse_crs(crs))) {
      stop(sprintf(
        "crs \"%s\" is not a coordinate reference system sf can read", crs
      ))
    }
    return(crs)
  }

  if (!is.numeric(crs) || length(crs) != 1L || is.na(crs) ||
      !is.finite(crs) || crs <= 0) {
    stop("crs must be a single positive EPSG code, or a PROJ or WKT string")
  }

  as.integer(crs)
}

#' CRS a grid's coordinates are read in
#'
#' A grid carrying no CRS is read on its own body: 'WGS84' on Earth, and a
#' longlat CRS on the sphere of the radius elsewhere.
#' @param x HexGridInfo object or legacy hexify_grid list
#' @noRd
grid_crs <- function(x) {
  crs <- if (isS4(x)) {
    if (.hasSlot(x, "crs")) x@crs else NULL
  } else {
    x$crs
  }
  if (is.null(crs) || length(crs) != 1L || is.na(crs)) {
    crs <- resolve_crs(NULL, grid_radius_km(x), grid_ellipsoid(x))
  }
  parse_crs(crs)
}

#' A CRS as it prints in a grid specification
#' @param crs An 'EPSG' code or a 'PROJ' or 'WKT' string
#' @noRd
format_crs <- function(crs) {
  if (is.character(crs)) crs else sprintf("EPSG:%d", as.integer(crs))
}

#' A CRS as it reads in a message
#'
#' The SRID where there is one, and whatever sf was given otherwise.
#' @param crs An object of class crs
#' @noRd
crs_label <- function(crs) {
  for (field in c("srid", "input", "proj4string")) {
    label <- crs[[field]]
    if (!is.null(label) && length(label) == 1L && !is.na(label) &&
        nzchar(label)) {
      return(label)
    }
  }
  "an unnamed CRS"
}

#' Are two CRSs on the same body?
#'
#' The semi-major axis is the body: Earth's reference systems all carry the
#' same one, and a grid on another body carries that body's radius. A CRS that
#' does not report one is taken as being on the same body, leaving the question
#' to 'PROJ'.
#' @param a,b Objects of class crs
#' @noRd
same_body_crs <- function(a, b) {
  radii <- suppressWarnings(c(as.numeric(a$SemiMajor), as.numeric(b$SemiMajor)))
  if (length(radii) != 2L || !all(is.finite(radii))) {
    return(TRUE)
  }
  isTRUE(all.equal(radii[1], radii[2], tolerance = 1e-6))
}

#' Geometry read in another CRS
#'
#' A grid on another body carries a longlat CRS on that body's sphere, so the
#' boundary a caller supplies, or the raster a caller extracts from, sits under
#' a CRS 'PROJ' will not reach: it holds no operation between two celestial
#' bodies. Longitude and latitude are the same angles on either sphere, though,
#' and name the same region there, so lon/lat geometry is read in the target
#' CRS and the caller is told. Anything else carries lengths, which are not the
#' same, and is refused with both CRSs named. Two CRSs on one body reproject as
#' they always have.
#'
#' @param x An sf or sfc object
#' @param target The CRS to read it in, as sf reads a CRS
#' @param what The calling function, for the messages
#' @param subject What `x` is, for the messages
#' @param counterpart What `target` belongs to, for the messages
#' @return `x`, in `target`
#' @noRd
geometry_in_crs <- function(x, target, what, subject, counterpart) {
  source <- sf::st_crs(x)

  if (is.na(source) || is.na(target) || source == target) {
    return(x)
  }

  if (same_body_crs(source, target)) {
    return(sf::st_transform(x, target))
  }

  if (!isTRUE(sf::st_is_longlat(source)) || !isTRUE(sf::st_is_longlat(target))) {
    stop(what, "(): ", subject, " is in ", crs_label(source), " and the ",
         counterpart, " in ", crs_label(target),
         ", which are on different bodies, and a projected CRS carries lengths ",
         "that do not carry over. Give both in longitude and latitude, which ",
         "names the same region on either body.", call. = FALSE)
  }

  message(what, "(): ", subject, " is in ", crs_label(source), " and the ",
          counterpart, " in ", crs_label(target), ". Longitude and latitude ",
          "name the same region on both, so ", subject, " is read in the ",
          counterpart, "'s CRS.")
  # sf warns that replacing a CRS does not reproject, which is the point
  suppressWarnings(sf::st_set_crs(x, target))
}

#' CRS of a terra raster, as sf reads a CRS
#' @param raster A terra SpatRaster
#' @noRd
raster_crs <- function(raster) {
  parse_crs(terra::crs(raster))
}

#' 'H3' areas read on another body
#'
#' A cell covers the same solid angle on any sphere, and 'H3' reports an area as
#' that solid angle times H3_EARTH_RADIUS_KM squared, so dividing that radius
#' back out and multiplying by the body's leaves the solid angle times the
#' body's radius squared. Earth returns H3's own figure untouched: hexify's
#' Earth radius is the mean radius and H3's the authalic one, and rescaling
#' between the two would move an Earth area no user asked to move.
#' @param area_km2 Numeric vector of 'H3' areas, in km^2
#' @param radius_km Radius in kilometers
#' @noRd
scale_area_to_body <- function(area_km2, radius_km) {
  if (radius_km == EARTH_RADIUS_KM) return(area_km2)
  area_km2 * (radius_km / H3_EARTH_RADIUS_KM)^2
}

# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Get resolution from grid object (handles both field names)
#'
#' Extracts resolution from a grid object, supporting both 'resolution'
#' (hexify style) and 'res' ('dggridR' style) field names.
#'
#' @param dggs Grid specification object
#' @param require Logical; if TRUE, stops with error if resolution not found
#' @return Integer resolution value, or NULL if not found and require=FALSE
#' @noRd
get_grid_resolution <- function(dggs, require = FALSE) {
  if ("resolution" %in% names(dggs)) {
    dggs$resolution
  } else if ("res" %in% names(dggs)) {
    dggs$res
  } else if (require) {
    stop("Grid object missing resolution field (neither 'resolution' nor 'res')")
  } else {
    NULL
  }
}

# =============================================================================
# Unit Conversion Constants
# =============================================================================

#' Square miles to square kilometers conversion factor
#' @noRd
MI2_TO_KM2 <- 2.58999

#' Miles to kilometers conversion factor
#' @noRd
MI_TO_KM <- 1.60934

# =============================================================================
# ISEA Aperture-3 Calibration Constants
# =============================================================================

#' Cell area at effective resolution 10 (aperture 3) in km^2
#' Used for area-to-resolution conversions.
#' @noRd
ISEA3H_RES10_AREA_KM2 <- 863.8006

# =============================================================================
# ISEA Default Orientation Constants
# =============================================================================
# Standard ISEA orientation with vertex 0 positioned at these coordinates.
# Reference: Snyder (1992) "An Equal-Area Map Projection For Polyhedral Globes"

#' Default longitude for ISEA vertex 0 (degrees)
#' @noRd
ISEA_VERT0_LON_DEG <- 11.25

#' Default latitude for ISEA vertex 0 (degrees)
#'
#' atan(phi) in degrees, where phi = (1 + sqrt(5)) / 2.
#' @noRd
ISEA_VERT0_LAT_DEG <- 58.282525588538995

#' Default azimuth for ISEA orientation (degrees)
#' @noRd
ISEA_AZIMUTH_DEG <- 0.0

#' The standard ISEA orientation, as a grid's orientation slot holds it
#' @noRd
ISEA_ORIENTATION <- c(vert0_lon = ISEA_VERT0_LON_DEG,
                      vert0_lat = ISEA_VERT0_LAT_DEG,
                      azimuth = ISEA_AZIMUTH_DEG)

#' DGGRID's REGION_CENTER placement (SubOpDGG::orientGrid()): vertex 0 and
#' the point its azimuth is taken towards, as gnomonic coordinates in metres
#' about the region centre on the sphere of the WGS84 authalic radius
#' @noRd
DGGRID_REGION_VERT0_M <- c(-7289214.618283, 7289214.618283)
DGGRID_REGION_AZ_POINT_M <- c(2784232.232959, 2784232.232959)
DGGRID_AUTHALIC_RADIUS_M <- 6371007.180918475

# =============================================================================
# Grid Parameter Limits
# =============================================================================

#' Valid aperture values
#' @noRd
VALID_APERTURES <- c(3L, 4L, 7L)

#' Hex9's aperture: Griffin's (2026) shifted aperture 9, on the octahedron
#' only, and never mixed with the others
#' @noRd
HEX9_APERTURE <- 9L

#' Cells of a Hex9 grid at resolution 0, which multiply by 9 per resolution
#' @noRd
HEX9_BASE_CELLS <- 12L

#' Whether an aperture spelling names Hex9
#' @param aperture Character or numeric aperture spelling
#' @noRd
is_hex9_aperture <- function(aperture) {
  length(aperture) == 1L && identical(as.character(aperture), "9")
}

#' Whether a grid is a Hex9 grid
#' @param g HexGridInfo object
#' @noRd
is_hex9_grid <- function(g) {
  !is_h3_grid(g) && is_hex9_aperture(g@aperture)
}

#' Maximum supported resolution
#' @noRd
MAX_RESOLUTION <- 30L

#' Minimum supported resolution
#' @noRd
MIN_RESOLUTION <- 0L

# =============================================================================
# H3 Grid Constants
# =============================================================================

#' Earth radius the vendored H3 library measures areas and lengths against
#'
#' The 'WGS84' authalic radius, from EARTH_RADIUS_KM in src/h3/constants.h.
#' Distinct from EARTH_RADIUS_KM here, which is the 'WGS84' mean radius.
#' @noRd
H3_EARTH_RADIUS_KM <- 6371.007180918475

#' Maximum H3 resolution
#' @noRd
H3_MAX_RESOLUTION <- 15L

#' Minimum H3 resolution
#' @noRd
H3_MIN_RESOLUTION <- 0L

#' Average cell areas from H3 documentation (km^2), indexed by resolution + 1
#' Source: https://h3geo.org/docs/core-library/restable/
#' @noRd
H3_AVG_AREA_KM2 <- c(
  4357449.416,  # res 0
  609788.442,   # res 1
  86801.780,    # res 2
  12393.435,    # res 3
  1770.348,     # res 4
  252.904,      # res 5
  36.129,       # res 6
  5.161,        # res 7
  0.737,        # res 8
  0.105,        # res 9
  0.015,        # res 10
  0.00215,      # res 11
  0.000307,     # res 12
  0.0000439,    # res 13
  0.00000627,   # res 14
  0.000000895   # res 15
)

# =============================================================================
# H3 Resolution Helpers
# =============================================================================

#' Find closest H3 resolution for a target area
#'
#' Shared by hex_grid() and h3_crosswalk() to avoid duplicating the
#' resolution-matching logic.
#'
#' @param area_km2 Target cell area in km^2
#' @param radius_km Radius of the body, in kilometers
#' @return Integer H3 resolution (0-15)
#' @noRd
closest_h3_resolution <- function(area_km2, radius_km = EARTH_RADIUS_KM) {
  diffs <- abs(scale_area_to_body(H3_AVG_AREA_KM2, radius_km) - area_km2)
  which.min(diffs) - 1L
}

#' Average area of an H3 cell at a resolution, on a given body
#' @param resolution Integer H3 resolution (0-15)
#' @param radius_km Radius of the body, in kilometers
#' @noRd
h3_avg_area_km2 <- function(resolution, radius_km = EARTH_RADIUS_KM) {
  scale_area_to_body(H3_AVG_AREA_KM2[resolution + 1L], radius_km)
}

#' Check whether a grid is H3 type
#'
#' Safe check that handles old serialized objects without grid_type slot.
#'
#' @param grid A HexGridInfo object
#' @return Logical
#' @noRd
is_h3_grid <- function(grid) {
  tryCatch(
    identical(grid@grid_type, "h3"),
    error = function(e) FALSE
  )
}

# =============================================================================
# Coordinate Validation Helpers
# =============================================================================

#' Validate longitude values
#' @param lon Numeric vector of longitudes
#' @param warn Whether to warn on out-of-range values (default TRUE)
#' @return Logical vector indicating valid values
#' @noRd
validate_lon <- function(lon, warn = TRUE) {
  if (!is.numeric(lon)) {
    stop("Longitude must be numeric")
  }
  valid <- is.na(lon) | (lon >= -180 & lon <= 180)
  if (warn && any(!valid, na.rm = TRUE)) {
    warning("Some longitude values are outside valid range [-180, 180]")
  }
  valid
}

#' Validate latitude values
#' @param lat Numeric vector of latitudes
#' @param warn Whether to warn on out-of-range values (default TRUE)
#' @return Logical vector indicating valid values
#' @noRd
validate_lat <- function(lat, warn = TRUE) {
  if (!is.numeric(lat)) {
    stop("Latitude must be numeric")
  }
  valid <- is.na(lat) | (lat >= -90 & lat <= 90)
  if (warn && any(!valid, na.rm = TRUE)) {
    warning("Some latitude values are outside valid range [-90, 90]")
  }
  valid
}

#' Validate that vectors read in step carry the same length
#'
#' Names the calling function and the arguments, so a caller passing three
#' longitudes and one latitude reads which of its own arguments is short.
#'
#' @param values Named list of the vectors read in step
#' @param what The calling function, for the message
#' @return TRUE if valid, otherwise throws error
#' @noRd
validate_same_length <- function(values, what) {
  n <- lengths(values)
  if (length(unique(n)) > 1L) {
    stop(sprintf(
      "%s(): %s are read in step and must be the same length; got %s",
      what,
      paste(sprintf("`%s`", names(values)), collapse = " and "),
      paste(sprintf("%s = %d", names(values), n), collapse = ", ")
    ), call. = FALSE)
  }
  TRUE
}

#' Validate resolution value
#' @param resolution Integer resolution value
#' @return TRUE if valid, otherwise throws error
#' @noRd
validate_resolution <- function(resolution) {
  if (!is.numeric(resolution) || length(resolution) != 1) {
    stop("Resolution must be a single numeric value")
  }
  resolution <- as.integer(resolution)
  if (is.na(resolution) ||
      resolution < MIN_RESOLUTION ||
      resolution > MAX_RESOLUTION) {
    stop(sprintf(
      "Resolution must be between %d and %d", MIN_RESOLUTION, MAX_RESOLUTION
    ))
  }
  TRUE
}

#' Validate aperture value
#' @param aperture Integer aperture value
#' @return TRUE if valid, otherwise throws error
#' @noRd
validate_aperture <- function(aperture) {
  if (!is.numeric(aperture) || length(aperture) != 1) {
    stop("Aperture must be a single numeric value")
  }
  aperture <- as.integer(aperture)
  if (!aperture %in% VALID_APERTURES) {
    stop(sprintf(
      "Aperture must be one of: %s", paste(VALID_APERTURES, collapse = ", ")
    ))
  }
  TRUE
}

#' Hierarchical index type for a grid aperture
#' @param aperture Character aperture ("3", "4", "7", "9", or "4/3")
#' @return "z3", "z7", "hex9", or "zorder"
#' @noRd
index_type_for_aperture <- function(aperture) {
  if (aperture == "3") "z3"
  else if (aperture == "7") "z7"
  else if (is_hex9_aperture(aperture)) "hex9"
  else "zorder"
}

#' Integer aperture for the C++ functions that take a single one
#'
#' A mixed sequence has no single aperture, and the 3 it reports is not a
#' refinement factor of the grid. Anything that reads the aperture as one --
#' cell geometry, hierarchy, cell counts -- takes the grid's levels from
#' isea_levels() instead, which carry one aperture per level.
#' @param aperture Character or numeric aperture spelling
#' @return Integer aperture (3L, 4L, or 7L)
#' @noRd
aperture_to_int <- function(aperture) {
  if (is_mixed_aperture(aperture)) 3L else as.integer(aperture)
}

#' Validate cell ID values
#' @param cell_id Cell IDs, in any form as_cell_id() reads
#' @param resolution Integer resolution value
#' @param aperture Integer aperture value
#' @param warn Whether to warn on out-of-range values (default TRUE)
#' @param polyhedron The solid the grid is built on
#' @return Logical vector indicating valid values
#' @noRd
validate_cell_id <- function(cell_id, resolution, aperture, warn = TRUE,
                             polyhedron = "icosahedron") {
  cell_id <- as_cell_id(cell_id)
  max_id <- isea_cell_count(aperture, resolution, polyhedron)
  valid <- is.na(cell_id) | (cell_id >= 1L & cell_id <= max_id)
  if (warn && any(!valid, na.rm = TRUE)) {
    warning(sprintf(
      "Some cell IDs are outside valid range [1, %s] for res %d, ap %d",
      as.character(max_id), resolution, aperture
    ))
  }
  valid
}
