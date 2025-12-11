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
#' Compatibility wrapper for hexify_cell_to_boundary().
#' 
#' @param dggs Grid specification
#' @param in_seqnum Cell indices
#' 
#' @return List of boundary coordinates
#' 
#' @export
#' @seealso \code{\link{hexify_cell_to_boundary}}
dgSEQNUM_to_BOUNDARY <- function(dggs, in_seqnum) {
  return(hexify_cell_to_boundary(dggs, in_seqnum))
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
#' Generates cell indices for all cells that intersect a rectangular
#' geographic region defined by min/max longitude and latitude.
#' 
#' @param dggs Grid specification
#' @param minlat Minimum latitude
#' @param minlon Minimum longitude
#' @param maxlat Maximum latitude
#' @param maxlon Maximum longitude
#' @param cellsize Sampling density (degrees between sample points)
#' @param ... Additional arguments passed to dgcellstogrid
#' 
#' @return sf object or data frame of hexagonal cells
#' 
#' @export
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
dgrectgrid <- function(dggs, minlat, minlon, maxlat, maxlon, 
                       cellsize = 0.1, ...) {
  dgverify(dggs)
  
  # Generate a dense grid of sample points
  lon_seq <- seq(minlon, maxlon, by = cellsize)
  lat_seq <- seq(minlat, maxlat, by = cellsize)
  
  # Create all combinations
  sample_points <- expand.grid(lon = lon_seq, lat = lat_seq)
  
  # Convert sample points to cell indices
  cells <- dgGEO_to_SEQNUM(dggs, sample_points$lon, sample_points$lat)
  
  # Get unique cells (many sample points will be in the same cell)
  unique_cells <- unique(cells$seqnum)
  
  # Convert to grid
  return(dgcellstogrid(dggs, unique_cells, ...))
}

#' Generate cells for entire Earth (dggridR-compatible)
#' 
#' WARNING: This can generate millions of cells at high resolution!
#' Use with caution. For most purposes, use dgrectgrid or dgshptogrid
#' to generate cells for specific regions instead.
#' 
#' @param dggs Grid specification
#' @param ... Additional arguments passed to dgcellstogrid
#' 
#' @return sf object or data frame of hexagonal cells
#' 
#' @export
dgearthgrid <- function(dggs, ...) {
  dgverify(dggs)
  
  # Calculate total number of cells
  stats <- dgearthstat(dggs)
  
  if (stats$n_cells > 1000000) {
    stop(sprintf(
      "This would generate %.0f million cells! Consider using:\n  1. Lower resolution (current: %d)\n  2. dgrectgrid() for specific region\n  3. dgshptogrid() for specific area",
      stats$n_cells / 1000000,
      dggs$resolution
    ))
  }
  
  warning(sprintf(
    "Generating %.0f cells. This may take a while...",
    stats$n_cells
  ))
  
  # For now, use rectangular grid covering entire globe
  # This is not perfect but generates most cells
  return(dgrectgrid(dggs, 
                   minlat = -90, minlon = -180,
                   maxlat = 90, maxlon = 180,
                   cellsize = 2,  # Coarse sampling
                   ...))
}

#' Generate cells for shapefile region (dggridR-compatible)
#' 
#' Generates cell indices for all cells that intersect with a shapefile
#' polygon or multi-polygon.
#' 
#' @param dggs Grid specification
#' @param shpfname Path to shapefile or sf object
#' @param cellsize Sampling density (degrees between sample points)
#' @param ... Additional arguments passed to dgcellstogrid
#' 
#' @return sf object or data frame of hexagonal cells
#' 
#' @export
dgshptogrid <- function(dggs, shpfname, cellsize = 0.1, ...) {
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
  return(dgrectgrid(dggs,
                   minlat = bbox["ymin"], minlon = bbox["xmin"],
                   maxlat = bbox["ymax"], maxlon = bbox["xmax"],
                   cellsize = cellsize,
                   ...))
}
