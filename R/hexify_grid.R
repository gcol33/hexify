# hexify_grid.R
# Core grid construction and validation functions
#
# This file contains the fundamental grid construction functions that form
# the foundation of the hexify package.

#' @title Core Grid Construction
#' @description Core functions for hexify grid construction and validation
#' @name hexify-grid
NULL

#' Calculate resolution for target area
#'
#' Uses the 'ISEA3H' cell count formula: N = 10 * aperture^res + 2
#' This matches 'dggridR' resolution numbering exactly.
#'
#' @param target_area_km2 Target area in square kilometers
#' @param aperture Aperture (3, 4, or 7)
#' @return Resolution level
#' @keywords internal
calculate_resolution_for_area <- function(target_area_km2, aperture = 3) {
  # ISEA3H cell count formula (matches dggridR exactly):
  # N = 10 * aperture^res + 2
  #
  # Solving for res given target area:
  # area = EARTH_SURFACE / N
  # N = EARTH_SURFACE / area
  # 10 * aperture^res + 2 = EARTH_SURFACE / area
  # aperture^res = (EARTH_SURFACE / area - 2) / 10
  # res = log((EARTH_SURFACE / area - 2) / 10) / log(aperture)

  n_cells <- EARTH_SURFACE_KM2 / target_area_km2
  resolution <- log((n_cells - 2) / 10) / log(aperture)

  return(resolution)  # Return unrounded for caller to handle rounding mode
}

#' Create a hexagonal grid specification
#'
#' Creates a discrete global grid system (DGGS) object with hexagonal cells
#' at a specified resolution. This is the main constructor for hexify grids.
#'
#' @param area Target cell area in km^2 (if metric=TRUE) or area code
#' @param topology Grid topology (only "HEXAGON" supported)
#' @param metric Whether area is in metric units (km^2)
#' @param resround How to round resolution ("nearest", "up", "down")
#' @param aperture Aperture sequence (3, 4, or 7)
#' @param projection Projection type (only 'ISEA' supported currently)
#'
#' @return A hexify_grid object containing:
#'   \item{area}{Target cell area}
#'   \item{resolution}{Calculated resolution level}
#'   \item{aperture}{Grid aperture (3, 4, or 7)}
#'   \item{topology}{Grid topology ("HEXAGON")}
#'   \item{projection}{Map projection ("ISEA")}
#'   \item{index_type}{Index encoding type ("z3", "z7", or "zorder")}
#'
#' @family hexify main
#' @seealso \code{\link{hexify}} for the main user function,
#'   \code{\link{hexify_grid_to_cell}} for coordinate conversion
#' @export
#' @examples
#' # Create a grid with ~1000 km^2 cells
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' print(grid)
#'
#' # Create a finer resolution grid (~100 km^2 cells)
#' fine_grid <- hexify_grid(area = 100, aperture = 3, resround = "up")
hexify_grid <- function(area, 
                             topology = "HEXAGON", 
                             metric = TRUE,
                             resround = "nearest",
                             aperture = 3,
                             projection = "ISEA") {
  
  # Input validation
  if (topology != "HEXAGON") {
    stop("Only HEXAGON topology is supported")
  }
  
  if (projection != "ISEA") {
    stop("Only ISEA projection is supported")
  }
  
  validate_aperture(aperture)
  
  if (!resround %in% c("nearest", "up", "down")) {
    stop("resround must be 'nearest', 'up', or 'down'")
  }
  
  # Calculate resolution for target area
  resolution <- calculate_resolution_for_area(area, aperture)
  
  # Apply rounding
  if (resround == "up") {
    resolution <- ceiling(resolution)
  } else if (resround == "down") {
    resolution <- floor(resolution)
  } else {
    resolution <- round(resolution)
  }
  
  # Ensure resolution is valid
  resolution <- max(MIN_RESOLUTION, min(MAX_RESOLUTION, resolution))
  
  # Initialize icosahedron geometry (C++ function)
  cpp_build_icosa()
  
  # Determine index type based on aperture
  index_type <- if (aperture == 3) {
    "z3"
  } else if (aperture == 7) {
    "z7"
  } else {
    "zorder"
  }
  
  # Create grid specification with both hexify and dggridR-compatible fields
  grid <- list(
    # Hexify fields
    area = area,
    resolution = resolution,
    aperture = aperture,
    topology = topology,
    projection = projection,
    metric = metric,
    index_type = index_type,
    
    # dggridR-compatible fields (for backwards compatibility)
    res = resolution,
    topology_family = topology,
    metric_radius = if (metric) sqrt(area / pi) else NULL,
    pole_lon_deg = ISEA_VERT0_LON_DEG,
    pole_lat_deg = ISEA_VERT0_LAT_DEG,
    azimuth_deg = ISEA_AZIMUTH_DEG,
    aperture_type = "SEQUENCE",
    res_spec = resolution,
    precision = 7
  )
  
  # Set class for method dispatch
  class(grid) <- c("hexify_grid", "dggs", "list")
  
  return(grid)
}


#' Verify grid object
#'
#' Validates that a grid object has all required fields and valid values.
#' This function is called internally by most hexify functions to ensure
#' grid integrity.
#'
#' @param dggs Grid object to verify (from hexify_grid)
#' @return TRUE (invisibly) if valid, otherwise throws an error
#'
#' @export
#' @examples
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' dgverify(grid)  # Should pass silently
#'
#' # Invalid grid will throw error
#' bad_grid <- list(aperture = 5)
#' try(dgverify(bad_grid))  # Will error
dgverify <- function(dggs) {
  # Check object type
  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a grid object from hexify_grid()")
  }
  
  # Check required fields exist
  required_fields <- c("aperture", "topology", "projection")
  missing_fields <- setdiff(required_fields, names(dggs))
  if (length(missing_fields) > 0) {
    stop(sprintf("Grid object missing required field(s): %s",
                 paste(missing_fields, collapse = ", ")))
  }
  
  resolution <- get_grid_resolution(dggs, require = TRUE)
  
  # Validate resolution and aperture
  validate_resolution(resolution)
  validate_aperture(dggs$aperture)
  
  # Validate topology
  if (dggs$topology != "HEXAGON") {
    warning("Only HEXAGON topology is fully supported")
  }
  
  # Validate projection
  if (dggs$projection != "ISEA") {
    warning("Only ISEA projection is fully supported")
  }
  
  invisible(TRUE)
}

