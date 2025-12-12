# hexify_cell_to_sf.R
# Cell ID to sf polygon conversion and grid generation
#
# This file provides efficient polygon generation from cell IDs (SEQNUM values),
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
#' polys <- hexify_cell_to_sf(result$hex_id, resolution = 10, aperture = 3)
#'
#' # Plot with sf
#' library(sf)
#' plot(st_geometry(polys), col = "lightblue", border = "blue")
#' }
hexify_cell_to_sf <- function(hex_id, resolution, aperture = 3L,
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
    corners_list <- cpp_cell_to_corners(hex_id, resolution, aperture)

    # Convert each matrix to sf polygon
    polygons <- lapply(corners_list, function(coords) {
      sf::st_polygon(list(coords))
    })

    sfc <- sf::st_sfc(polygons, crs = 4326)
    sf::st_sf(hex_id = hex_id, geometry = sfc)

  } else {
    # Return data frame format
    result <- cpp_cell_to_polygon(hex_id, resolution, aperture)
    names(result) <- c("hex_id", "lon", "lat", "order")
    result$hex_id <- as.integer(result$hex_id)
    result
  }
}

#' Backwards-compatible alias for hexify_cell_to_sf
#'
#' @rdname hexify_cell_to_sf
#' @export
hexify_polygons <- function(hex_id, resolution, aperture = 3L,
                            return_sf = TRUE) {
  hexify_cell_to_sf(hex_id, resolution, aperture, return_sf)
}

# =============================================================================
# HEXIFY RESULT TO SF CONVERSION
# =============================================================================

#' Convert hexify result to sf object
#'
#' Converts the output of \code{\link{hexify}()} to an sf spatial object.
#' Can return either point geometries (cell centers) or polygon geometries
#' (cell boundaries).
#'
#' @param data Data frame returned by \code{\link{hexify}()} containing
#'   hex_id, hex_cen_lon, hex_cen_lat columns
#' @param geometry Type of geometry: "point" for cell centers (default),
#'   "polygon" for cell boundaries
#' @param aperture Grid aperture: 3, 4, or 7. Only needed for polygon geometry.
#' @param crs Coordinate reference system (default 4326 = WGS84)
#'
#' @return sf object with:
#'   - All original columns from data
#'   - geometry column (sfc_POINT or sfc_POLYGON)
#'
#' @details
#' This is the recommended way to convert hexify output to sf for spatial
#' operations or plotting with ggplot2/tmap.
#'
#' For \code{geometry = "point"}: Creates point geometries at cell centers
#' using hex_cen_lon and hex_cen_lat. Fast, suitable for large datasets.
#'
#' For \code{geometry = "polygon"}: Creates polygon geometries for cell
#' boundaries. Slower but useful for choropleth maps. Duplicate cells are
#' automatically handled (each cell boundary appears once).
#'
#' @seealso \code{\link{hexify}} for the main function,
#'   \code{\link{hexify_to_polygons}} for polygon-only conversion
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(sf)
#'
#' # Create sample data
#' df <- data.frame(
#'   site = c("Paris", "Vienna", "Madrid"),
#'   lon = c(2.35, 16.37, -3.70),
#'   lat = c(48.86, 48.21, 40.42),
#'   value = c(100, 200, 150)
#' )
#'
#' # Hexify the data
#' result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
#'
#' # Convert to sf points (cell centers)
#' sf_points <- hexify_to_sf(result)
#' plot(sf_points["value"])
#'
#' # Convert to sf polygons (cell boundaries)
#' sf_polys <- hexify_to_sf(result, geometry = "polygon")
#' plot(sf_polys["value"])
#'
#' # Use with ggplot2
#' library(ggplot2)
#' ggplot(sf_polys) +
#'   geom_sf(aes(fill = value)) +
#'   theme_minimal()
#' }
hexify_to_sf <- function(data, geometry = c("point", "polygon"),
                         aperture = 3L, crs = 4326) {

  geometry <- match.arg(geometry)

  # Validate input
  if (!"hex_id" %in% names(data)) {
    stop("data must contain 'hex_id' column (output from hexify())")
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }


  if (geometry == "point") {
    # Point geometry from cell centers
    if (!"hex_cen_lon" %in% names(data) || !"hex_cen_lat" %in% names(data)) {
      stop("data must contain 'hex_cen_lon' and 'hex_cen_lat' columns ",
           "(output from hexify())")
    }

    sf::st_as_sf(data,
                 coords = c("hex_cen_lon", "hex_cen_lat"),
                 crs = crs)

  } else {
    # Polygon geometry from cell boundaries
    if (!"hex_area" %in% names(data)) {
      stop("data must contain 'hex_area' column for polygon geometry ",
           "(output from hexify())")
    }

    resolution <- .resolution_from_area(data$hex_area[1], aperture)

    # Get unique cells for polygon generation
    unique_ids <- unique(data$hex_id)
    polys_sf <- hexify_cell_to_sf(unique_ids, resolution,
                                   as.integer(aperture), return_sf = TRUE)

    # Merge original data attributes with polygons
    # Keep all rows from data, join polygon geometry by hex_id
    data_no_geom <- data
    result <- merge(polys_sf, data_no_geom, by = "hex_id", all.y = TRUE)

    # Ensure sf class is preserved
    sf::st_as_sf(result)
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
#' @return sf object or data frame with polygon geometries (see hexify_cell_to_sf)
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
#' polys <- hexify_to_polygons(result)
#' plot(st_geometry(polys))
#' }
hexify_to_polygons <- function(data, aperture = 3L, return_sf = TRUE) {

  if (!"hex_id" %in% names(data)) {
    stop("data must contain 'hex_id' column (output from hexify())")
  }
  if (!"hex_area" %in% names(data)) {
    stop("data must contain 'hex_area' column (output from hexify()). ",
         "Use hexify_polygons() directly if you know the resolution.")
  }

  resolution <- .resolution_from_area(data$hex_area[1], aperture)

  hexify_cell_to_sf(
    hex_id = unique(data$hex_id),
    resolution = resolution,
    aperture = as.integer(aperture),
    return_sf = return_sf
  )
}

#' Quick plot of hexify results
#'
#' Simple plotting function for hexify output using base R graphics.
#' For ggplot2/tmap, use hexify_to_polygons() instead.
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
#' hexify_plot(result, col = "lightblue", border = "darkblue")
#' points(df$lon, df$lat, pch = 19, col = "red")
#' }
hexify_plot <- function(data, aperture = 3L, col = "lightgray",
                        border = "black", add = FALSE, ...) {

  if (!"hex_id" %in% names(data)) {
    stop("data must contain 'hex_id' column (output from hexify())")
  }
  if (!"hex_area" %in% names(data)) {
    stop("data must contain 'hex_area' column (output from hexify())")
  }

  # Get polygon coordinates
  polys_df <- hexify_to_polygons(data, aperture = aperture, return_sf = FALSE)
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
#' grid <- hexify_grid_rect(
#'   minlon = -10, maxlon = 20,
#'   minlat = 35, maxlat = 60,
#'   area = 5000
#' )
#' plot(st_geometry(grid), border = "gray")
#' }
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
  result <- hexify(grid_pts, lon = "lon", lat = "lat", area = area,
                   aperture = aperture, resround = resround)

  hexify_to_polygons(result, aperture = as.integer(aperture))
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
#' global_grid <- hexify_grid_global(area = 100000)
#' plot(st_geometry(global_grid), border = "gray")
#' }
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

  result <- hexify(grid_pts, lon = "lon", lat = "lat", area = area,
                   aperture = aperture, resround = resround)

  hexify_to_polygons(result, aperture = as.integer(aperture))
}

# =============================================================================
# LOW-LEVEL SF HELPERS
# =============================================================================

#' Build an sf POLYGON from six (lon, lat) corner pairs
#'
#' @param lon numeric vector of length 6 (longitude)
#' @param lat numeric vector of length 6 (latitude)
#' @param crs integer CRS (default 4326)
#' @return sf object with one POLYGON geometry
#' @export
hex_corners_to_sf <- function(lon, lat, crs = 4326) {
  stopifnot(length(lon) == 6L, length(lat) == 6L)
  pts <- cbind(lon, lat)
  ring <- rbind(pts, pts[1, , drop = FALSE])  # close polygon
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(ring)), crs = crs))
}
