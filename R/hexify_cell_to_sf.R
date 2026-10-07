# hexify_cell_to_sf.R
# Cell ID to sf polygon conversion and grid generation
#
# This file provides efficient polygon generation from cell IDs,
# with sf integration for modern spatial workflows.

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
#' @param aperture Grid aperture: 3, 4, 7, or a mixed spelling as
#'   \code{\link{hex_grid}} takes it. Can be omitted if grid is provided.
#' @param return_sf Logical. If TRUE (default), returns sf object with polygon
#'   geometries. If FALSE, returns data frame with vertex coordinates.
#' @param grid Optional HexGridInfo object. If provided, resolution and aperture
#'   are extracted from it.
#' @inheritParams cell_to_sf
#'
#' @return If return_sf = TRUE: sf object with columns:
#'   \item{cell_id}{Cell identifier}
#'   \item{geometry}{POLYGON geometry (sfc_POLYGON)}
#'
#'   If return_sf = FALSE: data frame with columns:
#'   \item{cell_id}{Cell identifier}
#'   \item{lon}{Vertex longitude}
#'   \item{lat}{Vertex latitude}
#'   \item{order}{Vertex order along the ring; the last vertex repeats the
#'     first, closing the polygon}
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

  if (!is.null(grid)) {
    g <- extract_grid(grid)
  } else {
    if (is.null(resolution) || is.null(aperture)) {
      stop("resolution and aperture must be provided, or supply a grid object")
    }
    g <- hex_grid(resolution = resolution, aperture = aperture)
  }

  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  if (is_h3_grid(g) || return_sf) {
    return(cell_to_sf(cell_id, g, wrap_dateline = wrap_dateline))
  }

  cell_id <- unique(cell_id[!is.na(cell_id)])
  if (length(cell_id) == 0) {
    stop("No valid cell_id values provided")
  }

  rings <- isea_cell_rings(cell_id, g@resolution, g@aperture, icosa_arg(g))
  n_vertices <- vapply(rings, nrow, integer(1))
  vertices <- do.call(rbind, rings)
  data.frame(
    cell_id = rep(cell_id, n_vertices),
    lon = vertices[, 1],
    lat = vertices[, 2],
    order = sequence(n_vertices)
  )
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
#' @param radius_km Radius of the body, in kilometers, or a body name such as
#'   "mars" (default Earth). See \code{\link{hex_grid}}.
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
                             area, aperture = 3L, resround = "nearest",
                             radius_km = EARTH_RADIUS_KM) {
  grid <- hex_grid(area_km2 = area, aperture = aperture, resround = resround,
                   radius_km = radius_km)
  grid_rect(c(minlon, minlat, maxlon, maxlat), grid)
}

#' Generate a global grid of hexagon polygons
#'
#' Creates hexagon polygons covering a whole body.
#'
#' @param area Target cell area in km^2
#' @param aperture Grid aperture: 3, 4, or 7
#' @param resround Resolution rounding: "nearest", "up", or "down"
#' @param radius_km Radius of the body, in kilometers, or a body name such as
#'   "mars" (default Earth). See \code{\link{hex_grid}}.
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
hexify_grid_global <- function(area, aperture = 3L, resround = "nearest",
                               radius_km = EARTH_RADIUS_KM) {
  grid <- hex_grid(area_km2 = area, aperture = aperture, resround = resround,
                   radius_km = radius_km)
  grid_global(grid)
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
