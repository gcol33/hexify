# hexify_compat.R
# dggridR compatibility layer
#
# This file provides wrapper functions for compatibility with the dggridR
# package. All functions starting with 'dg' are dggridR-compatible wrappers
# around the native hexify functions.

#' @title dggridR Compatibility Layer
#' @description Wrapper functions for dggridR package compatibility
#' @name hexify-compat
NULL

#' Construct a discrete global grid system (dggridR-compatible)
#' 
#' Compatibility wrapper for hexify_construct() that follows dggridR naming
#' conventions. Creates a discrete global grid system (DGGS) object.
#' 
#' @param area Target cell area in km² (if metric=TRUE)
#' @param topology Grid topology (only "HEXAGON" supported)
#' @param metric Whether area is in metric units (km²)
#' @param resround How to round resolution ("nearest", "up", "down")
#' @param aperture Aperture sequence (3, 4, or 7)
#' @param projection Projection type (only "ISEA" supported currently)
#' 
#' @return A dggs object (same as hexify_grid)
#' 
#' @export
#' @seealso \code{\link{hexify_construct}}
#' @examples
#' \dontrun{
#' # dggridR-style usage
#' dggs <- dgconstruct(area = 1000, metric = TRUE, aperture = 3)
#' print(dggs)
#' 
#' # Equivalent hexify usage
#' grid <- hexify_construct(area = 1000, aperture = 3)
#' }
dgconstruct <- function(area, 
                       topology = "HEXAGON",
                       metric = TRUE, 
                       resround = "nearest",
                       aperture = 3,
                       projection = "ISEA") {
  return(hexify_construct(area, topology, metric, resround, aperture, projection))
}

#' Convert geographic coordinates to cell sequence numbers (dggridR-compatible)
#' 
#' Compatibility wrapper for hexify_geo_to_cell() that follows dggridR naming
#' conventions and return format. Converts lon/lat coordinates to cell indices.
#' 
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#' @param in_lon_deg Longitude vector in degrees
#' @param in_lat_deg Latitude vector in degrees
#' 
#' @return List with component:
#'   \item{seqnum}{Cell indices (character vector)}
#'   
#' @export
#' @seealso \code{\link{hexify_geo_to_cell}}
#' @examples
#' \dontrun{
#' dggs <- dgconstruct(area = 1000, aperture = 3)
#' 
#' # Single point
#' result <- dgGEO_to_SEQNUM(dggs, in_lon_deg = 0, in_lat_deg = 0)
#' print(result$seqnum)
#' 
#' # Multiple points
#' lons <- c(0, 10, 20)
#' lats <- c(0, 45, -30)
#' result <- dgGEO_to_SEQNUM(dggs, lons, lats)
#' print(result$seqnum)
#' }
dgGEO_to_SEQNUM <- function(dggs, in_lon_deg, in_lat_deg) {
  result <- hexify_geo_to_cell(dggs, in_lon_deg, in_lat_deg)
  return(list(seqnum = result$seqnum))
}

#' Convert cell sequence numbers to geographic coordinates (dggridR-compatible)
#' 
#' Compatibility wrapper for hexify_cell_to_geo() that follows dggridR naming
#' conventions and return format. Converts cell indices back to lon/lat.
#' 
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#' @param in_seqnum Cell indices (character vector)
#' 
#' @return List with components:
#'   \item{lon_deg}{Longitude in degrees}
#'   \item{lat_deg}{Latitude in degrees}
#'   
#' @export
#' @seealso \code{\link{hexify_cell_to_geo}}
#' @examples
#' \dontrun{
#' dggs <- dgconstruct(area = 1000, aperture = 3)
#' 
#' # Get some cell indices
#' cells <- dgGEO_to_SEQNUM(dggs, in_lon_deg = c(0, 10), 
#'                                in_lat_deg = c(0, 45))
#' 
#' # Convert back to coordinates
#' coords <- dgSEQNUM_to_GEO(dggs, cells$seqnum)
#' print(coords)
#' }
dgSEQNUM_to_GEO <- function(dggs, in_seqnum) {
  result <- hexify_cell_to_geo(dggs, in_seqnum)
  return(list(lon_deg = result$lon_deg, lat_deg = result$lat_deg))
}

#' Get cell boundaries (dggridR-compatible)
#'
#' Returns boundary coordinates for hexagonal cells as an sf object
#' or data frame of vertex coordinates.
#'
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#' @param in_seqnum Cell indices (integer SEQNUM values)
#' @param return_sf Logical. If TRUE, returns sf object; if FALSE, returns
#'   data frame with vertex coordinates
#'
#' @return If return_sf = TRUE: sf object with polygon geometries
#'   If return_sf = FALSE: data frame with columns seqnum, lon, lat, order
#'
#' @export
#' @seealso \code{\link{hex_polygons}}
dgSEQNUM_to_BOUNDARY <- function(dggs, in_seqnum, return_sf = TRUE) {
  hex_polygons(
    hex_id = as.integer(in_seqnum),
    resolution = dggs$resolution,
    aperture = dggs$aperture,
    return_sf = return_sf
  )
}

#' Set grid resolution (dggridR-compatible)
#' 
#' Changes the resolution of an existing grid object. This is useful for
#' creating multiple grids at different resolutions from the same base
#' configuration.
#' 
#' @param dggs Grid specification from dgconstruct()
#' @param res New resolution level (0-30)
#' 
#' @return Updated grid object with new resolution
#' 
#' @export
#' @examples
#' \dontrun{
#' # Create a grid
#' dggs <- dgconstruct(area = 1000, aperture = 3)
#' print(dggs$resolution)  # e.g., resolution 6
#' 
#' # Change to higher resolution
#' dggs_fine <- dgsetres(dggs, res = 10)
#' print(dggs_fine$resolution)  # resolution 10
#' }
dgsetres <- function(dggs, res) {
  # Validate resolution
  if (!is.numeric(res) || res < 0 || res > 30) {
    stop("Resolution must be a number between 0 and 30")
  }
  
  # Update resolution in grid object
  dggs$resolution <- as.integer(res)
  dggs$res <- as.integer(res)  # Also update 'res' field for compatibility
  
  # Recalculate dependent fields using ISEA3H formula: N = 10 * aperture^res + 2
  n_cells <- 10 * (dggs$aperture ^ res) + 2
  dggs$area <- EARTH_SURFACE_KM2 / n_cells
  
  # Verify the updated grid
  dgverify(dggs)
  
  return(dggs)
}

#' Get information about grid resolution (dggridR-compatible)
#' 
#' Returns information about grid properties at the current resolution,
#' including cell area, spacing, and total number of cells.
#' 
#' @param dggs Grid specification
#' @param res Optional resolution level (uses grid's resolution if not specified)
#' 
#' @return List with grid information
#' 
#' @export
#' @seealso \code{\link{dgearthstat}}
dggetres <- function(dggs, res = NULL) {
  if (!is.null(res)) {
    dggs <- dgsetres(dggs, res)
  }
  
  return(dgearthstat(dggs))
}

#' Generate cells for a rectangular region (dggridR-compatible)
#'
#' Generates hexagon polygons covering a rectangular geographic region.
#'
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#' @param minlat Minimum latitude
#' @param minlon Minimum longitude
#' @param maxlat Maximum latitude
#' @param maxlon Maximum longitude
#'
#' @return sf object with hexagon polygons covering the specified region
#'
#' @export
#' @seealso \code{\link{hex_grid_rect}}
#' @examples
#' \dontrun{
#' dggs <- dgconstruct(area = 10000, aperture = 3)
#'
#' # Get cells for a region (e.g., central Europe)
#' grid <- dgrectgrid(dggs,
#'                   minlat = 45, minlon = 5,
#'                   maxlat = 55, maxlon = 15)
#'
#' # Plot the result
#' library(sf)
#' plot(st_geometry(grid))
#' }
dgrectgrid <- function(dggs, minlat, minlon, maxlat, maxlon) {
  dgverify(dggs)

  hex_grid_rect(
    minlon = minlon, maxlon = maxlon,
    minlat = minlat, maxlat = maxlat,
    area = dggs$area,
    aperture = dggs$aperture
  )
}

#' Generate cells for entire Earth (dggridR-compatible)
#'
#' Generates hexagon polygons covering the entire Earth.
#' WARNING: This can generate millions of cells at high resolution!
#'
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#'
#' @return sf object with hexagon polygons covering the globe
#'
#' @export
#' @seealso \code{\link{hex_grid_global}}
dgearthgrid <- function(dggs) {
  dgverify(dggs)

  hex_grid_global(
    area = dggs$area,
    aperture = dggs$aperture
  )
}

#' Generate cells for shapefile region (dggridR-compatible)
#'
#' Generates hexagon polygons covering the bounding box of a shapefile
#' or sf object.
#'
#' @param dggs Grid specification from dgconstruct() or hexify_construct()
#' @param shpfname Path to shapefile or sf object
#'
#' @return sf object with hexagon polygons covering the region
#'
#' @export
dgshptogrid <- function(dggs, shpfname) {
  dgverify(dggs)

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Read shapefile if it's a path
  if (is.character(shpfname)) {
    shp <- sf::st_read(shpfname, quiet = TRUE)
  } else if (inherits(shpfname, "sf")) {
    shp <- shpfname
  } else {
    stop("shpfname must be a file path or sf object")
  }

  # Get bounding box

  bbox <- sf::st_bbox(shp)

  # Generate cells for bounding box
  dgrectgrid(dggs,
             minlat = bbox["ymin"], minlon = bbox["xmin"],
             maxlat = bbox["ymax"], maxlon = bbox["xmax"])
}
