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

# =============================================================================
# Body Geometry Helpers
# =============================================================================

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
  grid_radius_km(x) == EARTH_RADIUS_KM
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
#' @noRd
ISEA_VERT0_LAT_DEG <- 58.28252559

#' Default azimuth for ISEA orientation (degrees)
#' @noRd
ISEA_AZIMUTH_DEG <- 0.0

# =============================================================================
# Grid Parameter Limits
# =============================================================================

#' Valid aperture values
#' @noRd
VALID_APERTURES <- c(3L, 4L, 7L)

#' Maximum supported resolution
#' @noRd
MAX_RESOLUTION <- 30L

#' Minimum supported resolution
#' @noRd
MIN_RESOLUTION <- 0L

# =============================================================================
# H3 Grid Constants
# =============================================================================

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
# Session-Scoped Cache
# =============================================================================

# =============================================================================
# H3 Resolution Helpers
# =============================================================================

#' Find closest H3 resolution for a target area
#'
#' Shared by hex_grid() and h3_crosswalk() to avoid duplicating the
#' resolution-matching logic.
#'
#' @param area_km2 Target cell area in km^2
#' @return Integer H3 resolution (0-15)
#' @noRd
closest_h3_resolution <- function(area_km2) {
  diffs <- abs(H3_AVG_AREA_KM2 - area_km2)
  which.min(diffs) - 1L
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
#' @param aperture Character aperture ("3", "4", "7", or "4/3")
#' @return "z3", "z7", or "zorder"
#' @noRd
index_type_for_aperture <- function(aperture) {
  if (aperture == "3") "z3"
  else if (aperture == "7") "z7"
  else "zorder"
}

#' Integer aperture for the C++ functions that take a single one
#'
#' A mixed sequence has no single aperture; its cell IDs come from the sequence
#' API instead, and the callers that still ask for an integer use it only to
#' select an unrotated substrate branch, which aperture 3 gives.
#' @param aperture Character or numeric aperture spelling
#' @return Integer aperture (3L, 4L, or 7L)
#' @noRd
aperture_to_int <- function(aperture) {
  if (is_mixed_aperture(aperture)) 3L else as.integer(aperture)
}

#' Calculate maximum cell ID for given resolution and aperture
#' @param resolution Integer resolution value
#' @param aperture Integer aperture value
#' @return Maximum valid cell ID (numeric)
#' @noRd
max_cell_id <- function(resolution, aperture) {
  # Cell count formula: N = 10 * aperture^res + 2 (cell IDs are 1..N,
  # matching calc_grid_params() in rcpp_cell.cpp)
  10 * (aperture^resolution) + 2
}

#' Validate cell ID values
#' @param cell_id Numeric vector of cell IDs
#' @param resolution Integer resolution value
#' @param aperture Integer aperture value
#' @param warn Whether to warn on out-of-range values (default TRUE)
#' @return Logical vector indicating valid values
#' @noRd
validate_cell_id <- function(cell_id, resolution, aperture, warn = TRUE) {
  if (!is.numeric(cell_id)) {
    stop("Cell ID must be numeric")
  }
  max_id <- max_cell_id(resolution, aperture)
  valid <- is.na(cell_id) | (cell_id >= 1 & cell_id <= max_id)
  if (warn && any(!valid, na.rm = TRUE)) {
    warning(sprintf(
      "Some cell IDs are outside valid range [1, %.0f] for res %d, ap %d",
      max_id, resolution, aperture
    ))
  }
  valid
}
