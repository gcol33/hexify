# hex_polygons.R
# Native polygon generation for hexify
#
# This file provides efficient polygon generation from hex_id values,
# with sf integration for modern spatial workflows.

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Get resolution from area (internal)
#' @noRd
.resolution_from_area <- function(area, aperture) {
  grid <- hexify_construct(area = area, aperture = as.integer(aperture),
                           resround = "nearest")
  grid$resolution
}

# =============================================================================
# CORE POLYGON GENERATION
# =============================================================================

#' Generate hexagon polygons from cell IDs
#'
#' Creates polygon geometries for hexagonal grid cells from their SEQNUM
#' identifiers. Returns an sf object by default, or a data frame for
#' lightweight workflows.
#'
#' @param hex_id Integer vector of cell identifiers (SEQNUM values)
#' @param resolution Grid resolution level
#' @param aperture Grid aperture: 3, 4, or 7
#' @param return_sf Logical. If TRUE (default), returns sf object with polygon
#'   geometries. If FALSE, returns data frame with vertex coordinates.
#'
#' @return If return_sf = TRUE: sf object with columns:
#'   \item{hex_id}{Cell identifier}
#'   \item{geometry}{POLYGON geometry (sfc_POLYGON)}
#'
#'   If return_sf = FALSE: data frame with columns:
#'   \item{hex_id}{Cell identifier}
#'   \item{lon}{Vertex longitude}
#'   \item{lat}{Vertex latitude}
#'   \item{order}{Vertex order (1-7, 7 closes the polygon)}
#'
#' @details
#' This function uses a native C++ implementation that is significantly faster
#' than dggridR's polygon generation, especially for large numbers of cells.
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#'
#' # Generate some data with hex cells
#' df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 45))
#' result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
#'
#' # Get polygons as sf object
#' polys <- hex_polygons(result$hex_id, resolution = 10, aperture = 3)
#'
#' # Plot with sf
#' library(sf)
#' plot(st_geometry(polys), col = "lightblue", border = "blue")
#' }
hex_polygons <- function(hex_id, resolution, aperture = 3L,
                         return_sf = TRUE) {

  # Input validation
  if (!is.numeric(hex_id)) {
    stop("hex_id must be numeric (integer SEQNUM values)")
  }
  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  # Remove NA values and duplicates
  hex_id <- unique(hex_id[!is.na(hex_id)])
  if (length(hex_id) == 0) {
    stop("No valid hex_id values provided")
  }

  aperture <- as.integer(aperture)
  resolution <- as.integer(resolution)

  if (return_sf) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for return_sf = TRUE. ",
           "Install with: install.packages('sf')")
    }

    # Use the list-based corner function for efficient sf construction
    corners_list <- cpp_seqnum_to_corners_dggrid(hex_id, resolution, aperture)

    # Convert each matrix to sf polygon
    polygons <- lapply(corners_list, function(coords) {
      sf::st_polygon(list(coords))
    })

    sfc <- sf::st_sfc(polygons, crs = 4326)
    sf::st_sf(hex_id = hex_id, geometry = sfc)

  } else {
    # Return data frame format
    result <- cpp_seqnum_to_polygon_dggrid(hex_id, resolution, aperture)
    names(result) <- c("hex_id", "lon", "lat", "order")
    result$hex_id <- as.integer(result$hex_id)
    result
  }
}

# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

#' Generate polygons directly from hexify result
#'
#' Convenience function that extracts resolution from a hexify result and
#' generates polygons. Resolution is auto-detected from the hex_area column.
#'
#' @param data Data frame returned by hexify() containing hex_id and hex_area
#' @param aperture Grid aperture (default 3)
#' @param return_sf Logical. If TRUE (default), returns sf object
#'
#' @return sf object or data frame with polygon geometries (see hex_polygons)
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(sf)
#'
#' # Simple workflow: hexify then get polygons
#' df <- data.frame(lon = runif(100, -10, 10), lat = runif(100, 40, 50))
#' result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
#'
#' # Get polygons - resolution auto-detected
#' polys <- hex_to_polygons(result)
#' plot(st_geometry(polys))
#' }
hex_to_polygons <- function(data, aperture = 3L, return_sf = TRUE) {

  if (!"hex_id" %in% names(data)) {
    stop("data must contain 'hex_id' column (output from hexify())")
  }
  if (!"hex_area" %in% names(data)) {
    stop("data must contain 'hex_area' column (output from hexify()). ",
         "Use hex_polygons() directly if you know the resolution.")
  }

  resolution <- .resolution_from_area(data$hex_area[1], aperture)

  hex_polygons(
    hex_id = unique(data$hex_id),
    resolution = resolution,
    aperture = as.integer(aperture),
    return_sf = return_sf
  )
}

#' Quick plot of hexify results
#'
#' Simple plotting function for hexify output using base R graphics.
#' For ggplot2/tmap, use hex_to_polygons() instead.
#'
#' @param data Data frame returned by hexify()
#' @param aperture Grid aperture (default 3)
#' @param col Fill color for hexagons
#' @param border Border color for hexagons
#' @param add If TRUE, add to existing plot
#' @param ... Additional arguments passed to polygon()
#'
#' @return NULL (invisibly). Creates a plot as side effect.
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#'
#' df <- data.frame(lon = c(16.37, 2.35, -3.70), lat = c(48.21, 48.86, 40.42))
#' result <- hexify(df, lon = "lon", lat = "lat", area = 5000)
#'
#' # Quick plot
#' hex_plot(result, col = "lightblue", border = "darkblue")
#' points(df$lon, df$lat, pch = 19, col = "red")
#' }
hex_plot <- function(data, aperture = 3L, col = "lightgray",
                     border = "black", add = FALSE, ...) {

  if (!"hex_id" %in% names(data)) {
    stop("data must contain 'hex_id' column (output from hexify())")
  }
  if (!"hex_area" %in% names(data)) {
    stop("data must contain 'hex_area' column (output from hexify())")
  }

  # Get polygon coordinates
  polys_df <- hex_to_polygons(data, aperture = aperture, return_sf = FALSE)
  unique_ids <- unique(polys_df$hex_id)

  if (!add) {
    plot(range(polys_df$lon), range(polys_df$lat), type = "n",
         xlab = "Longitude", ylab = "Latitude",
         asp = 1 / cos(mean(polys_df$lat) * pi / 180))
  }

  for (id in unique_ids) {
    idx <- polys_df$hex_id == id
    polygon(polys_df$lon[idx], polys_df$lat[idx],
            col = col, border = border, ...)
  }

  invisible(NULL)
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
#' @param area Target cell area in km²
#' @param aperture Grid aperture: 3, 4, or 7
#' @param resround Resolution rounding: "nearest", "up", or "down"
#'
#' @return sf object with hexagon polygons covering the specified region
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(sf)
#'
#' grid <- hex_grid_rect(
#'   minlon = -10, maxlon = 20,
#'   minlat = 35, maxlat = 60,
#'   area = 5000
#' )
#' plot(st_geometry(grid), border = "gray")
#' }
hex_grid_rect <- function(minlon, maxlon, minlat, maxlat,
                          area, aperture = 3L, resround = "nearest") {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Create grid of sample points
  diagonal <- sqrt(area * 2 / sqrt(3))
  spacing_deg <- diagonal / 111 * 0.8

  lons <- seq(minlon, maxlon, by = spacing_deg)
  lats <- seq(minlat, maxlat, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  # Assign to hexes and get polygons
  result <- hexify(grid_pts, lon = "lon", lat = "lat", area = area,
                   aperture = aperture, resround = resround)

  hex_to_polygons(result, aperture = as.integer(aperture))
}

#' Generate a global grid of hexagon polygons
#'
#' Creates hexagon polygons covering the entire Earth.
#'
#' @param area Target cell area in km²
#' @param aperture Grid aperture: 3, 4, or 7
#' @param resround Resolution rounding: "nearest", "up", or "down"
#'
#' @return sf object with hexagon polygons covering the globe
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(sf)
#'
#' # Coarse global grid (~100,000 km² cells)
#' global_grid <- hex_grid_global(area = 100000)
#' plot(st_geometry(global_grid), border = "gray")
#' }
hex_grid_global <- function(area, aperture = 3L, resround = "nearest") {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Estimate cell count

  EARTH_SURFACE_KM2 <- 510065623
  approx_cells <- EARTH_SURFACE_KM2 / area

  if (approx_cells > 100000) {
    warning(sprintf(
      "This will generate approximately %.0f cells. Consider larger area.",
      approx_cells
    ))
  }

  # Generate dense sample points
  diagonal <- sqrt(area * 2 / sqrt(3))
  spacing_deg <- diagonal / 111 * 0.7

  lons <- seq(-180, 180, by = spacing_deg)
  lats <- seq(-85, 85, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  result <- hexify(grid_pts, lon = "lon", lat = "lat", area = area,
                   aperture = aperture, resround = resround)

  hex_to_polygons(result, aperture = as.integer(aperture))
}
