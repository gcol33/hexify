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

  # Get cell centers. `cells`/`boundary`, when supplied, take precedence over
  # a HexData grid's own already-assigned cells (matching hex_zonal()) --
  # otherwise they'd be silently ignored whenever `grid` is a HexData object.
  if (!is.null(cells)) {
    cell_ids <- cells
    ll <- cell_to_lonlat(cells, g)
    centers <- cbind(lon = ll$lon_deg, lat = ll$lat_deg)
  } else if (!is.null(boundary)) {
    # Generate cells within boundary. grid_clip() returns polygon geometries
    # (a cell_id column, no lon/lat), so cell centers are recomputed the same
    # way as the 'cells' branch above.
    cell_data <- grid_clip(boundary, g)
    cell_ids <- cell_data$cell_id
    ll <- cell_to_lonlat(cell_ids, g)
    centers <- cbind(lon = ll$lon_deg, lat = ll$lat_deg)
  } else if (is_hex_data(grid)) {
    cell_ids <- grid@cell_id
    centers <- grid@cell_center
  } else {
    stop("Provide a HexData object, cell IDs via 'cells', or a 'boundary' polygon")
  }

  # Keep unique cells only
  uid <- !duplicated(cell_ids)
  u_cell_ids <- cell_ids[uid]
  u_centers <- centers[uid, , drop = FALSE]

  # Extract raster values at cell centers
  pts <- terra::vect(u_centers, crs = grid_crs_wkt(g))
  extracted <- terra::extract(raster, pts)

  # Build result
  result <- data.frame(cell_id = u_cell_ids, stringsAsFactors = FALSE)
  # Drop the ID column from terra::extract output
  if ("ID" %in% names(extracted)) {
    extracted$ID <- NULL
  }
  result <- cbind(result, extracted)
  result
}
