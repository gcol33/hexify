# hexify_cell_to_sf.R
# Cell ID to sf polygon conversion and grid generation
#
# This file provides efficient polygon generation from cell IDs,
# with sf integration for modern spatial workflows.

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Get resolution from area (internal)
#' @noRd
.resolution_from_area <- function(area, aperture) {
  grid <- hexify_grid(area = area, aperture = as.integer(aperture),
                      resround = "nearest")
  grid$resolution
}

# =============================================================================
# CORE POLYGON GENERATION
# =============================================================================

#' Convert cell IDs to sf polygons
#'
#' Creates polygon geometries for hexagonal grid cells from their cell IDs.
#' Returns an sf object by default, or a data frame for lightweight workflows.
#'
#' @param cell_id Integer vector of cell identifiers
#' @param resolution Grid resolution level. Can be omitted if grid is provided.
#' @param aperture Grid aperture: 3, 4, or 7. Can be omitted if grid is provided.
#' @param return_sf Logical. If TRUE (default), returns sf object with polygon
#'   geometries. If FALSE, returns data frame with vertex coordinates.
#' @param grid Optional HexGridInfo object. If provided, resolution and aperture
#'   are extracted from it.
#' @param wrap_dateline Logical. If TRUE (default), calls
#'   \code{sf::st_wrap_dateline()} to split antimeridian-crossing polygons.
#'   Set to FALSE for orthographic/globe projections where wrapping creates gaps.
#'
#' @return If return_sf = TRUE: sf object with columns:
#'   \item{cell_id}{Cell identifier}
#'   \item{geometry}{POLYGON geometry (sfc_POLYGON)}
#'
#'   If return_sf = FALSE: data frame with columns:
#'   \item{cell_id}{Cell identifier}
#'   \item{lon}{Vertex longitude}
#'   \item{lat}{Vertex latitude}
#'   \item{order}{Vertex order (1-7, 7 closes the polygon)}
#'
#' @details
#' This function uses a native C++ implementation that is significantly faster
#' than 'dggridR' polygon generation, especially for large numbers of cells.
#'
#' For the recommended S4 interface, use \code{\link{cell_to_sf}} instead.
#'
#' @family sf conversion
#' @seealso \code{\link{cell_to_sf}} for the recommended S4 interface
#' @keywords internal
#' @export
#' @examples
#' library(hexify)
#'
#' # Generate some data with hex cells
#' df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 45))
#' result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
#'
#' # Get polygons as sf object (using HexData)
#' polys <- cell_to_sf(grid = result)
#'
#' # Or with explicit parameters
#' polys <- hexify_cell_to_sf(result@cell_id, resolution = 10, aperture = 3)
#'
#' # Plot with sf
#' library(sf)
#' plot(st_geometry(polys), col = "lightblue", border = "blue")
hexify_cell_to_sf <- function(cell_id, resolution = NULL, aperture = NULL,
                              return_sf = TRUE, grid = NULL,
                              wrap_dateline = TRUE) {

  # Extract from grid if provided
  if (!is.null(grid)) {
    g <- extract_grid(grid)
    # H3 grids: redirect to cell_to_sf() which handles H3 boundaries
    if (is_h3_grid(g)) {
      return(cell_to_sf(cell_id, g))
    }
    resolution <- g@resolution
    aperture <- as.integer(g@aperture)
  }

  # Input validation
  if (!is.numeric(cell_id)) {
    stop("cell_id must be numeric (integer cell IDs)")
  }
  if (is.null(resolution) || is.null(aperture)) {
    stop("resolution and aperture must be provided, or supply a grid object")
  }
  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  # Remove NA values and duplicates
  cell_id <- unique(cell_id[!is.na(cell_id)])
  if (length(cell_id) == 0) {
    stop("No valid cell_id values provided")
  }

  aperture <- as.integer(aperture)
  resolution <- as.integer(resolution)

  if (return_sf) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for return_sf = TRUE. ",
           "Install with: install.packages('sf')")
    }

    # Use the list-based corner function for efficient sf construction
    corners_list <- cpp_cell_to_corners(cell_id, resolution, aperture)

    # Handle antimeridian-crossing polygons by normalizing coordinates
    polygons <- lapply(corners_list, function(coords) {
      lons <- coords[, 1]
      lon_range <- max(lons, na.rm = TRUE) - min(lons, na.rm = TRUE)

      if (lon_range > 180) {
        lons[lons < 0] <- lons[lons < 0] + 360
        coords[, 1] <- lons
        mean_lon <- mean(lons)
        if (mean_lon > 180) {
          coords[, 1] <- coords[, 1] - 360
        }
      }

      sf::st_polygon(list(coords))
    })

    sfc <- sf::st_sfc(polygons, crs = 4326)
    sfc <- sf::st_make_valid(sfc)
    result_sf <- sf::st_sf(cell_id = cell_id, geometry = sfc)
    if (wrap_dateline) {
      result_sf <- sf::st_wrap_dateline(result_sf,
        options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
    }
    result_sf

  } else {
    # Return data frame format
    result <- cpp_cell_to_polygon(cell_id, resolution, aperture)
    names(result) <- c("cell_id", "lon", "lat", "order")
    result$cell_id <- as.integer(result$cell_id)
    result
  }
}


# =============================================================================
# GRID GENERATION
# =============================================================================

#' Generate a rectangular grid of hexagon polygons
#'
#' Creates hexagon polygons covering a rectangular geographic region.
#'
#' @param minlon,maxlon Longitude bounds
#' @param minlat,maxlat Latitude bounds
#' @param area Target cell area in km^2
#' @param aperture Grid aperture: 3, 4, or 7
#' @param resround Resolution rounding: "nearest", "up", or "down"
#'
#' @return sf object with hexagon polygons covering the specified region
#'
#' @family sf conversion
#' @seealso \code{\link{grid_rect}} for the recommended S4 interface,
#'   \code{\link{hexify_grid_global}} for global grids
#' @keywords internal
#' @export
#' @examples
#' library(hexify)
#' library(sf)
#'
#' grid <- hexify_grid_rect(
#'   minlon = -10, maxlon = 20,
#'   minlat = 35, maxlat = 60,
#'   area = 5000
#' )
#' plot(st_geometry(grid), border = "gray")
hexify_grid_rect <- function(minlon, maxlon, minlat, maxlat,
                             area, aperture = 3L, resround = "nearest") {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Create grid of sample points
  diagonal <- sqrt(area * 2 / sqrt(3))
  spacing_deg <- diagonal / KM_PER_DEGREE * 0.8

  lons <- seq(minlon, maxlon, by = spacing_deg)
  lats <- seq(minlat, maxlat, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  # Assign to hexes and get polygons
  result <- hexify(grid_pts, lon = "lon", lat = "lat", area_km2 = area,
                   aperture = aperture, resround = resround)

  cell_to_sf(grid = result)
}

#' Generate a global grid of hexagon polygons
#'
#' Creates hexagon polygons covering the entire Earth.
#'
#' @param area Target cell area in km^2
#' @param aperture Grid aperture: 3, 4, or 7
#' @param resround Resolution rounding: "nearest", "up", or "down"
#'
#' @return sf object with hexagon polygons covering the globe
#'
#' @family sf conversion
#' @seealso \code{\link{grid_global}} for the recommended S4 interface,
#'   \code{\link{hexify_grid_rect}} for regional grids
#' @keywords internal
#' @export
#' @examples
#' library(hexify)
#' library(sf)
#'
#' # Coarse global grid (~100,000 km^2 cells)
#' global_grid <- hexify_grid_global(area = 100000)
#' plot(st_geometry(global_grid), border = "gray")
hexify_grid_global <- function(area, aperture = 3L, resround = "nearest") {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Estimate cell count (uses EARTH_SURFACE_KM2 from constants.R)
  approx_cells <- EARTH_SURFACE_KM2 / area

  if (approx_cells > 100000) {
    warning(sprintf(
      "This will generate approximately %.0f cells. Consider larger area.",
      approx_cells
    ))
  }

  # Generate dense sample points
  diagonal <- sqrt(area * 2 / sqrt(3))
  spacing_deg <- diagonal / KM_PER_DEGREE * 0.7

  lons <- seq(-180, 180, by = spacing_deg)
  lats <- seq(-85, 85, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  result <- hexify(grid_pts, lon = "lon", lat = "lat", area_km2 = area,
                   aperture = aperture, resround = resround)

  cell_to_sf(grid = result)
}

# =============================================================================
# LOW-LEVEL SF HELPERS
# =============================================================================

#' Build an sf POLYGON from six (lon, lat) corner pairs
#'
#' Low-level helper to create a single hexagon polygon from corner coordinates.
#' Most users should use \code{\link{cell_to_sf}} instead.
#'
#' @param lon numeric vector of length 6 (longitude)
#' @param lat numeric vector of length 6 (latitude)
#' @param crs integer CRS (default 4326)
#' @return sf object with one POLYGON geometry
#'
#' @family sf conversion
#' @keywords internal
#' @export
hex_corners_to_sf <- function(lon, lat, crs = 4326) {
  stopifnot(length(lon) == 6L, length(lat) == 6L)
  pts <- cbind(lon, lat)
  ring <- rbind(pts, pts[1, , drop = FALSE])  # close polygon
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = crs))
}
