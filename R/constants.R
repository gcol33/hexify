# R/constants.R - Package-wide constants and helpers
#
# Centralizes magic numbers to avoid duplication and improve maintainability.

# =============================================================================
# Earth Geometry Constants
# =============================================================================

#' Total Earth surface area in square kilometers
#' @noRd
EARTH_SURFACE_KM2 <- 510072000

# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Get resolution from grid object (handles both field names)
#'
#' Extracts resolution from a grid object, supporting both 'resolution'
#' (hexify style) and 'res' (dggridR style) field names.
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
