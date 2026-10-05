# hexify_grid.R
# Core grid construction and validation functions
#
# This file contains the fundamental grid construction functions that form
# the foundation of the hexify package.

#' @title Core Grid Construction
#' @description Core functions for hexify grid construction and validation
#' @name hexify-grid
NULL

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
#' @param projection Face projection, DGGRID's \code{dggs_proj}: 'ISEA'
#'   (default) or 'FULLER'
#' @param radius_km Radius of the body the grid covers, in kilometers, or a body
#'   name such as "mars" (default Earth). See \code{\link{hex_grid}}.
#'
#' @return A hexify_grid object containing:
#'   \item{area}{Target cell area}
#'   \item{resolution}{Calculated resolution level}
#'   \item{aperture}{Grid aperture (3, 4, or 7)}
#'   \item{topology}{Grid topology ("HEXAGON")}
#'   \item{projection}{Face projection ("ISEA" or "FULLER")}
#'   \item{radius_km}{Radius of the body, in kilometers}
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
                             projection = "ISEA",
                             radius_km = EARTH_RADIUS_KM) {

  # Input validation
  if (topology != "HEXAGON") {
    stop("Only HEXAGON topology is supported")
  }

  if (!projection %in% DGGS_PROJECTIONS) {
    stop("projection must be 'ISEA' or 'FULLER'")
  }

  validate_aperture(aperture)

  radius_km <- resolve_radius_km(radius_km)

  if (!is.numeric(area) || length(area) != 1 || is.na(area) || area <= 0) {
    stop("area must be a positive number")
  }

  resolution <- resolve_resolution_from_area(area, aperture, radius_km, resround)
  index_type <- index_type_for_aperture(as.character(aperture))

  # Create grid specification with both hexify and dggridR-compatible fields
  grid <- list(
    # Hexify fields
    area = area,
    resolution = resolution,
    aperture = aperture,
    topology = topology,
    projection = projection,
    metric = metric,
    radius_km = radius_km,
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
#' @param dggs Grid object to verify (from hexify_grid() or hex_grid())
#' @return TRUE (invisibly) if valid, otherwise throws an error
#'
#' @export
#' @examples
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' dgverify(grid)  # Should pass silently
#'
#' # Modern HexGridInfo objects are accepted too
#' dgverify(hex_grid(area_km2 = 1000))
#'
#' # Invalid grid will throw error
#' bad_grid <- list(aperture = 5)
#' try(dgverify(bad_grid))  # Will error
dgverify <- function(dggs) {
  # Accept the modern HexGridInfo (S4) object by converting to the legacy
  # list representation this function checks.
  if (is_hex_grid(dggs)) {
    if (dggs@grid_type == "h3") {
      validate_resolution(dggs@resolution)
      return(invisible(TRUE))
    }
    dggs <- HexGridInfo_to_hexify_grid(dggs)
  }

  # Check object type
  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a grid object from hexify_grid() or hex_grid()")
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
  if (!dggs$projection %in% DGGS_PROJECTIONS) {
    warning("Only the ISEA and FULLER projections are supported")
  }
  
  invisible(TRUE)
}

