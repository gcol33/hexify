# R/hex_extract.R
# Raster extraction at hex cell centers

#' Extract Raster Values at Hex Cell Centers
#'
#' Samples raster values at hexagonal cell centers. Faster than
#' [hex_zonal()] because it only queries cell center points, not full
#' polygons.
#'
#' @param raster A `terra::SpatRaster` object.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#' @param cells Optional cell IDs to extract. If `NULL` (default), extracts
#'   at all cell centers from `grid` (only works for HexData objects or
#'   if a `boundary` is provided).
#' @param boundary Optional sf polygon to limit extraction extent.
#'
#' @return A data.frame with columns `cell_id`, plus one column per raster
#'   layer.
#'
#' @details
#' Requires the `terra` package (in Suggests). The function:
#' 1. Generates cell centers for the raster extent
#' 2. Calls `terra::extract(raster, cell_center_matrix)`
#' 3. Attaches cell IDs
#'
#' For full zonal statistics (aggregating all pixels within each hex polygon),
#' use [hex_zonal()] instead.
#'
#' @seealso [hex_zonal()] for polygon-based zonal statistics,
#'   [hexify()] for creating HexData objects
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   # Create a small synthetic raster
#'   r <- terra::rast(nrows = 10, ncols = 10,
#'                    xmin = -10, xmax = 10, ymin = 40, ymax = 55)
#'   terra::values(r) <- runif(100)
#'   names(r) <- "temperature"
#'
#'   # Extract at hex cell centers
#'   g <- hex_grid(area_km2 = 500)
#'   df <- data.frame(lon = c(0, 5), lat = c(45, 50))
#'   hd <- hexify(df, lon = "lon", lat = "lat", grid = g)
#'   hex_extract(r, hd)
#' }
#' }
hex_extract <- function(raster, grid, cells = NULL, boundary = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for hex_extract()")
  }

  g <- extract_grid(grid)

  u_cell_ids <- raster_target_cells(grid, g, cells, boundary)
  ll <- cell_to_lonlat(u_cell_ids, g)

  # Extract raster values at cell centers, read in the raster's own CRS
  pts <- sf::st_as_sf(data.frame(lon = ll$lon_deg, lat = ll$lat_deg),
                      coords = c("lon", "lat"), crs = grid_crs(g))
  pts <- geometry_in_crs(pts, raster_crs(raster), "hex_extract",
                         "the cell-center geometry", "raster")
  extracted <- terra::extract(raster, terra::vect(pts))

  # Build result
  result <- data.frame(cell_id = u_cell_ids, stringsAsFactors = FALSE)
  # Drop the ID column from terra::extract output
  if ("ID" %in% names(extracted)) {
    extracted$ID <- NULL
  }
  result <- cbind(result, extracted)
  result
}

#' Cells a raster summary runs over
#'
#' Explicit `cells` come first, then the cells meeting `boundary`, then a
#' HexData object's own cells.
#'
#' @param grid The HexGridInfo or HexData object the caller was given
#' @param g The HexGridInfo extracted from `grid`
#' @param cells Optional cell IDs
#' @param boundary Optional sf polygon
#' @return Unique, non-missing cell IDs
#' @noRd
raster_target_cells <- function(grid, g, cells = NULL, boundary = NULL) {
  ids <- if (!is.null(cells)) {
    cells
  } else if (!is.null(boundary)) {
    grid_clip(boundary, g)$cell_id
  } else if (is_hex_data(grid)) {
    grid@cell_id
  } else {
    stop("Provide a HexData object, cell IDs via 'cells', or a 'boundary' polygon")
  }
  unique(ids[!is.na(ids)])
}
