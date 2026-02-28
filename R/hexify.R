# hexify.R
# Main user-facing convenience function
#
# This is the primary entry point for the hexify package.

#' Assign hexagonal DGGS cell IDs to geographic points
#'
#' Takes a data.frame or sf object with geographic coordinates and returns
#' a HexData object that stores the original data plus cell assignments.
#' The original data is preserved unchanged; cell IDs and centers are stored
#' in separate slots.
#'
#' @param data A data.frame or sf object containing coordinates
#' @param grid A HexGridInfo object from \code{hex_grid()}. If provided, overrides
#'   area_km2, resolution, and aperture parameters.
#' @param lon Column name for longitude (ignored if data is sf)
#' @param lat Column name for latitude (ignored if data is sf)
#' @param area_km2 Target cell area in km^2 (mutually exclusive with diagonal).
#' @param diagonal Target cell diagonal (long diagonal) in km
#' @param resolution Grid resolution (0-30). Alternative to area_km2.
#' @param aperture Grid aperture: 3, 4, 7, or "4/3" for mixed (default 3)
#' @param resround How to round resolution: "nearest", "up", or "down"
#'
#' @return A HexData object containing:
#'   \itemize{
#'     \item \code{data}: The original input data (unchanged)
#'     \item \code{grid}: The HexGridInfo specification
#'     \item \code{cell_id}: Numeric vector of cell IDs for each row
#'     \item \code{cell_center}: Matrix of cell center coordinates (lon, lat)
#'   }
#'
#'   Use \code{as.data.frame(result)} to extract the original data.
#'   Use \code{cells(result)} to get unique cell IDs.
#'   Use \code{result@@cell_id} to get all cell IDs.
#'   Use \code{result@@cell_center} to get cell center coordinates.
#'
#' @details
#' For sf objects, coordinates are automatically extracted and transformed to
#' 'WGS84' (EPSG:4326) if needed. The geometry column is preserved.
#'
#' Either \code{area_km2} (or \code{area}), \code{diagonal}, or \code{resolution}
#' must be provided unless a \code{grid} object is supplied.
#'
#' The HexData return type (default) stores the grid specification so downstream
#' functions like \code{plot()}, \code{hexify_cell_to_sf()}, etc. don't need
#' grid parameters repeated.
#'
#' @section Grid Specification:
#' You can create a grid specification once and reuse it:
#' \preformatted{
#' grid <- hex_grid(area_km2 = 1000)
#' result1 <- hexify(df1, grid = grid)
#' result2 <- hexify(df2, grid = grid)
#' }
#'
#' @family hexify main
#' @seealso \code{\link{hex_grid}} for grid specification,
#'   \code{\link{HexData-class}} for return object details,
#'   \code{\link{as_sf}} for converting to sf
#' @export
#' @examples
#' # Simple data.frame
#' df <- data.frame(
#'   site = c("Vienna", "Paris", "Madrid"),
#'   lon = c(16.37, 2.35, -3.70),
#'   lat = c(48.21, 48.86, 40.42)
#' )
#'
#' # New recommended workflow: use grid object
#' grid <- hex_grid(area_km2 = 1000)
#' result <- hexify(df, grid = grid, lon = "lon", lat = "lat")
#' print(result)  # Shows grid info
#' plot(result)   # Plot with default styling
#'
#' # Direct area specification (grid created internally)
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Extract plain data.frame
#' df_result <- as.data.frame(result)
#'
#' # With sf object (any CRS)
#' library(sf)
#' pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
#' result_sf <- hexify(pts, area_km2 = 1000)
#'
#' # Different apertures
#' result_ap4 <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 4)
#'
#' # Mixed aperture (ISEA43H)
#' result_mixed <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = "4/3")
hexify <- function(data,
                   grid = NULL,
                   lon = "lon",
                   lat = "lat",
                   area_km2 = NULL,
                   diagonal = NULL,
                   resolution = NULL,
                   aperture = 3,
                   resround = "nearest") {

  # -------------------------------------------------------------------------
  # Extract or build grid specification
  # -------------------------------------------------------------------------
  if (!is.null(grid)) {
    # Grid object provided - extract parameters
    if (is_hex_grid(grid)) {
      hex_grid_obj <- grid
    } else if (inherits(grid, "hexify_grid")) {
      # Legacy S3 grid object
      hex_grid_obj <- hexify_grid_to_HexGridInfo(grid)
    } else {
      stop("grid must be a HexGridInfo object from hex_grid() or legacy hexify_grid")
    }
  } else {
    # Build grid from parameters
    if (is.null(area_km2) && is.null(diagonal) && is.null(resolution)) {
      stop("Either 'grid', 'area_km2', 'diagonal', or 'resolution' must be provided")
    }
    if (!is.null(area_km2) && !is.null(diagonal)) {
      stop("Provide either 'area_km2' or 'diagonal', not both")
    }

    # Convert diagonal to area if provided
    if (!is.null(diagonal)) {
      area_km2 <- diagonal^2 * sqrt(3) / 2
    }

    # Create HexGridInfo object (hex_grid handles aperture parsing)
    hex_grid_obj <- hex_grid(
      area_km2 = area_km2,
      resolution = resolution,
      aperture = aperture,
      resround = resround
    )
  }

  # -------------------------------------------------------------------------
  # Extract coordinates from data
  # -------------------------------------------------------------------------
  is_sf <- inherits(data, "sf")
  mapping <- list()

  if (is_sf) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to process sf objects")
    }

    # Get coordinates, transforming to WGS84 if needed
    if (sf::st_crs(data)$epsg != 4326 && !is.na(sf::st_crs(data)$epsg)) {
      coords_sf <- sf::st_transform(data, 4326)
    } else {
      coords_sf <- data
    }

    coords <- sf::st_coordinates(coords_sf)
    lon_vec <- coords[, 1]
    lat_vec <- coords[, 2]
    mapping$geometry <- attr(data, "sf_column")
  } else {
    # Regular data.frame
    if (!lon %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", lon))
    }
    if (!lat %in% names(data)) {
      stop(sprintf("Column '%s' not found in data", lat))
    }

    lon_vec <- data[[lon]]
    lat_vec <- data[[lat]]
    mapping$lon <- lon
    mapping$lat <- lat
  }

  # Validate coordinates
  if (!is.numeric(lon_vec) || !is.numeric(lat_vec)) {
    stop("Coordinates must be numeric")
  }

  na_mask <- is.na(lon_vec) | is.na(lat_vec)
  if (all(na_mask)) {
    stop("All coordinates are NA")
  }
  if (any(na_mask)) {
    warning(sprintf("%d coordinate pairs contain NA values and will be skipped",
                    sum(na_mask)))
  }

  # -------------------------------------------------------------------------
  # Perform hexification
  # -------------------------------------------------------------------------
  aperture_str <- hex_grid_obj@aperture
  res <- hex_grid_obj@resolution

  if (is_h3_grid(hex_grid_obj)) {
    # H3 path: use native C backend
    cell_ids <- cpp_h3_latLngToCell(lon_vec, lat_vec, res)
    center_df <- cpp_h3_cellToLatLng(cell_ids)
    centers <- list(lon_deg = center_df$lon, lat_deg = center_df$lat)
  } else if (aperture_str == "4/3") {
    level <- as.integer(res / 2)
    cell_ids <- cpp_lonlat_to_cell_ap43(lon_vec, lat_vec, res, level)
    centers <- cpp_cell_to_lonlat_ap43(cell_ids, res, level)
  } else {
    aperture_num <- as.integer(aperture_str)
    cell_ids <- cpp_lonlat_to_cell(lon_vec, lat_vec, res, aperture_num)
    centers <- cpp_cell_to_lonlat(cell_ids, res, aperture_num)
  }

  # Build cell center matrix
  cell_center <- matrix(
    c(centers$lon_deg, centers$lat_deg),
    ncol = 2,
    dimnames = list(NULL, c("lon", "lat"))
  )

  # -------------------------------------------------------------------------
  # Return HexData object (original data unchanged)
  # -------------------------------------------------------------------------
  new_hex_data(
    data = data,
    grid = hex_grid_obj,
    cell_id = cell_ids,
    cell_center = cell_center
  )
}

