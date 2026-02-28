# R/constants.R - Package-wide constants and helpers
#
# Centralizes magic numbers to avoid duplication and improve maintainability.

# =============================================================================
# Earth Geometry Constants
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

#' Calculate maximum cell ID for given resolution and aperture
#' @param resolution Integer resolution value
#' @param aperture Integer aperture value
#' @return Maximum valid cell ID (numeric)
#' @noRd
max_cell_id <- function(resolution, aperture) {
  # Cell count formula: N = 10 * aperture^res + 2
  # But for cell numbering, we use 20 faces with aperture^res cells each
  # Max cell ID = 20 * (max_coord + 1)^2 for the simple linear indexing

  # This matches the C++ cell_numbering.cpp logic
  if (resolution == 0) return(20)
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
