# R/dggrid_compat.R
# dggridR compatibility layer
#
# Functions for converting between hexify and dggridR grid objects.
# These enable interoperability with dggridR without copying its API.

# =============================================================================
# GRID OBJECT CONVERSION
# =============================================================================

#' Convert hexify grid to 'dggridR'-compatible grid object
#'
#' Creates a 'dggridR'-compatible grid specification from a hexify_grid object.
#' The resulting object can be used with 'dggridR' functions that accept a dggs
#' object.
#'
#' @param grid A hexify_grid object from hexify_grid()
#'
#' @return A list with 'dggridR'-compatible fields:
#'   \item{pole_lon_deg}{Longitude of grid pole (default 11.25)}
#'   \item{pole_lat_deg}{Latitude of grid pole (default 58.28252559)}
#'   \item{azimuth_deg}{Grid azimuth rotation (default 0)}
#'   \item{aperture}{Grid aperture (3, 4, or 7)}
#'   \item{res}{Resolution level}
#'   \item{topology}{Grid topology ("HEXAGON")}
#'   \item{projection}{Map projection ('ISEA')}
#'   \item{precision}{Output decimal precision (default 7)}
#'
#' @family 'dggridR' compatibility
#' @export
as_dggrid <- function(grid) {

  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }

  dggs <- list(
    pole_lon_deg = ISEA_VERT0_LON_DEG,
    pole_lat_deg = ISEA_VERT0_LAT_DEG,
    azimuth_deg = ISEA_AZIMUTH_DEG,
    aperture = grid$aperture,
    res = grid$resolution,
    topology = "HEXAGON",
    projection = "ISEA",
    precision = 7L
  )

  class(dggs) <- "list"
  dggs
}

#' Convert 'dggridR' grid object to hexify_grid
#'
#' Creates a hexify_grid object from a 'dggridR' dggs object. This allows
#' using hexify functions with grids created by 'dggridR' dgconstruct().
#'
#' @param dggs A 'dggridR' grid object from dgconstruct()
#'
#' @return A hexify_grid object
#'
#' @details
#' Only 'ISEA' projection with HEXAGON topology is fully supported.
#' Other configurations will generate warnings.
#'
#' The function validates that the 'dggridR' grid uses compatible settings:
#' - Projection must be 'ISEA' (FULLER not supported)
#' - Topology must be "HEXAGON" (DIAMOND, TRIANGLE not supported)
#' - Aperture must be 3, 4, or 7
#'
#' @family 'dggridR' compatibility
#' @export
from_dggrid <- function(dggs) {
  # Validate dggridR object

  if (!is.list(dggs)) {
    stop("dggs must be a list (dggridR grid object)")
  }

  required <- c("res", "aperture", "topology", "projection")
  missing <- setdiff(required, names(dggs))
  if (length(missing) > 0) {
    stop(sprintf("dggs missing required fields: %s", paste(missing, collapse = ", ")))
  }

  # Check compatibility

  if (!is.null(dggs$projection) && dggs$projection != "ISEA") {
    warning("Only ISEA projection is supported. Results may differ from dggridR.")
  }

  if (!is.null(dggs$topology) && dggs$topology != "HEXAGON") {
    warning("Only HEXAGON topology is supported. Results may differ from dggridR.")
  }

  if (!dggs$aperture %in% c(3L, 4L, 7L)) {
    stop(sprintf("Aperture %d not supported. Must be 3, 4, or 7.", dggs$aperture))
  }

  # Check for non-default orientation (not supported)
  if (!is.null(dggs$pole_lon_deg) && abs(dggs$pole_lon_deg - ISEA_VERT0_LON_DEG) > 1e-6) {
    warning("Non-default pole_lon_deg not supported. Using standard ISEA orientation.")
  }
  if (!is.null(dggs$pole_lat_deg) && abs(dggs$pole_lat_deg - ISEA_VERT0_LAT_DEG) > 1e-6) {
    warning("Non-default pole_lat_deg not supported. Using standard ISEA orientation.")
  }
  if (!is.null(dggs$azimuth_deg) && abs(dggs$azimuth_deg - ISEA_AZIMUTH_DEG) > 1e-6) {
    warning("Non-default azimuth_deg not supported. Using standard ISEA orientation.")
  }

  # Create hexify_grid
  hexify_grid(
    area = NA,  # Will be calculated from resolution if needed
    topology = "HEXAGON",
    metric = TRUE,
    resround = "nearest",
    aperture = as.integer(dggs$aperture),
    projection = "ISEA"
  ) -> grid

  # Override resolution to match dggridR exactly
  grid$resolution <- as.integer(dggs$res)
  grid$res <- as.integer(dggs$res)

  # Calculate actual area for this resolution
  n_cells <- 10 * (grid$aperture ^ grid$resolution) + 2
  grid$area <- EARTH_SURFACE_KM2 / n_cells

  grid
}

#' Validate 'dggridR' grid compatibility with hexify
#'
#' Checks whether a 'dggridR' grid object is compatible with hexify functions.
#' Returns TRUE if compatible, or throws an error describing incompatibilities.
#'
#' @param dggs A 'dggridR' grid object
#' @param strict If TRUE (default), throw errors for incompatibilities.
#'   If FALSE, return FALSE instead of throwing errors.
#'
#' @return TRUE if compatible, FALSE if not compatible (when strict=FALSE)
#'
#' @family 'dggridR' compatibility
#' @export
dggrid_is_compatible <- function(dggs, strict = TRUE) {
  issues <- character()

  if (!is.list(dggs)) {
    issues <- c(issues, "dggs must be a list")
  } else {
    if (is.null(dggs$projection) || dggs$projection != "ISEA") {
      issues <- c(issues, "Only ISEA projection supported (not FULLER)")
    }

    if (is.null(dggs$topology) || dggs$topology != "HEXAGON") {
      issues <- c(issues, "Only HEXAGON topology supported (not DIAMOND/TRIANGLE)")
    }

    if (is.null(dggs$aperture) || !dggs$aperture %in% c(3L, 4L, 7L)) {
      issues <- c(issues, "Aperture must be 3, 4, or 7")
    }

    # Check orientation
    if (!is.null(dggs$pole_lon_deg) && abs(dggs$pole_lon_deg - ISEA_VERT0_LON_DEG) > 1e-6) {
      issues <- c(issues, "Non-default pole_lon_deg not supported")
    }
    if (!is.null(dggs$pole_lat_deg) && abs(dggs$pole_lat_deg - ISEA_VERT0_LAT_DEG) > 1e-6) {
      issues <- c(issues, "Non-default pole_lat_deg not supported")
    }
    if (!is.null(dggs$azimuth_deg) && abs(dggs$azimuth_deg - ISEA_AZIMUTH_DEG) > 1e-6) {
      issues <- c(issues, "Non-default azimuth_deg not supported")
    }
  }

  if (length(issues) > 0) {
    if (strict) {
      stop(sprintf("dggridR grid not compatible with hexify:\n- %s",
                   paste(issues, collapse = "\n- ")))
    }
    return(FALSE)
  }

  TRUE
}
