# R/hex_zonal.R
# Full zonal statistics: aggregate raster pixels within hex polygons

#' Zonal Statistics for Hex Cells
#'
#' Computes zonal statistics by aggregating all raster pixels falling
#' within each hexagonal cell polygon. More accurate than [hex_extract()]
#' but slower because it requires polygon geometries.
#'
#' @param raster A `terra::SpatRaster` object.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#' @param fun Summary function name: `"mean"` (default), `"sum"`, `"min"`,
#'   `"max"`, `"sd"`, or `"count"`.
#' @param boundary Optional sf polygon to limit the analysis extent.
#' @param cells Optional cell IDs. If provided, only these cells are included.
#'
#' @return A data.frame with columns `cell_id`, plus one column per raster
#'   layer containing the aggregated values.
#'
#' @details
#' Requires the `terra` package (in Suggests). The function:
#' 1. Generates hex polygon geometries via [cell_to_sf()]
#' 2. Calls `terra::extract(raster, polygons, fun = fun)`
#' 3. Joins results back to cell IDs
#'
#' For point-based extraction (faster, at cell centers only), use
#' [hex_extract()] instead.
#'
#' @seealso [hex_extract()] for cell-center extraction,
#'   [cell_to_sf()] for hex polygon generation
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   r <- terra::rast(nrows = 100, ncols = 100,
#'                    xmin = -10, xmax = 10, ymin = 40, ymax = 55)
#'   terra::values(r) <- runif(10000)
#'   names(r) <- "temperature"
#'
#'   g <- hex_grid(area_km2 = 500)
#'   df <- data.frame(lon = c(0, 5), lat = c(45, 50))
#'   hd <- hexify(df, lon = "lon", lat = "lat", grid = g)
#'   hex_zonal(r, hd, fun = "mean")
#' }
#' }
hex_zonal <- function(raster, grid, fun = "mean", boundary = NULL,
                       cells = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for hex_zonal()")
  }

  g <- extract_grid(grid)

  # Get cell IDs
  if (!is.null(cells)) {
    cell_ids <- unique(cells[!is.na(cells)])
  } else if (is_hex_data(grid)) {
    cell_ids <- unique(grid@cell_id)
  } else if (!is.null(boundary)) {
    cell_data <- grid_clip(boundary, g)
    cell_ids <- unique(cell_data$cell_id)
  } else {
    stop("Provide a HexData object, cell IDs via 'cells', or a 'boundary' polygon")
  }

  # Generate hex polygons, read in the raster's own CRS
  hex_sf <- cell_to_sf(cell_ids, g)
  hex_vect <- terra::vect(geometry_in_crs(hex_sf, raster_crs(raster),
                                          "hex_zonal", "the grid's cell geometry",
                                          "raster"))

  # Extract zonal statistics. terra::extract() has no built-in "count"
  # function, so translate it to a count of non-NA pixels per cell.
  fun_actual <- if (identical(fun, "count")) {
    function(x, ...) sum(!is.na(x))
  } else {
    fun
  }
  extracted <- terra::extract(raster, hex_vect, fun = fun_actual, ID = TRUE)

  # Build result. Key off hex_sf$cell_id (not the original cell_ids) since
  # cell_to_sf() deduplicates internally and hex_vect/extracted follow its
  # row order -- using the pre-dedup cell_ids here would silently misalign
  # rows whenever the input contained duplicates.
  result <- data.frame(cell_id = hex_sf$cell_id, stringsAsFactors = FALSE)
  if ("ID" %in% names(extracted)) {
    extracted$ID <- NULL
  }
  result <- cbind(result, extracted)
  result
}
