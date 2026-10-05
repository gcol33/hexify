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
#'   area_km2, resolution, aperture and radius_km parameters.
#' @param lon Column name for longitude (ignored if data is sf)
#' @param lat Column name for latitude (ignored if data is sf)
#' @param area_km2 Target cell area in km^2 (mutually exclusive with diagonal).
#' @param diagonal Target cell spacing in km: the short (flat-to-flat)
#'   diagonal, equal to the distance between neighbouring cell centres
#' @param resolution Grid resolution (0-30). Alternative to area_km2.
#' @param aperture Grid aperture: 3, 4, 7, a mixed family such as "4/3" or
#'   "4/7", or one aperture per resolution level, e.g. \code{c(4, 4, 7, 3)}
#'   (default 3)
#' @param resround How to round resolution: "nearest", "up", or "down"
#' @param radius_km Radius of the body the grid covers, in kilometers, or a body
#'   name such as "mars" (default Earth). See \code{\link{hex_grid}}.
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
#' Either \code{area_km2}, \code{diagonal}, or \code{resolution}
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
#'   \code{\link[=st_as_sf.HexData]{st_as_sf}} for converting to sf
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
                   resround = "nearest",
                   radius_km = EARTH_RADIUS_KM) {

  # -------------------------------------------------------------------------
  # Extract or build grid specification
  # -------------------------------------------------------------------------
  if (!is.null(grid)) {
    hex_grid_obj <- extract_grid(grid)
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
      area_km2 <- hex_area_from_spacing(diagonal)
    }

    # Create HexGridInfo object (hex_grid handles aperture parsing)
    hex_grid_obj <- hex_grid(
      area_km2 = area_km2,
      resolution = resolution,
      aperture = aperture,
      resround = resround,
      radius_km = radius_km
    )
  }

  # -------------------------------------------------------------------------
  # Extract coordinates from data
  # -------------------------------------------------------------------------
  is_sf <- inherits(data, "sf")

  if (is_sf) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to process sf objects")
    }

    # The quantizer reads longitude and latitude on the grid's own body, so
    # that is the frame the coordinates are put in: WGS84 for an Earth grid,
    # and the body's own sphere elsewhere. Data carrying no CRS is taken as
    # already being in it.
    coords_sf <- geometry_in_crs(data, grid_crs(hex_grid_obj), "hexify",
                                 "the data", "grid")

    coords <- sf::st_coordinates(coords_sf)
    lon_vec <- coords[, 1]
    lat_vec <- coords[, 2]
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
  }

  # Validate coordinates
  if (!is.numeric(lon_vec) || !is.numeric(lat_vec)) {
    stop("Coordinates must be numeric")
  }
  validate_lon(lon_vec)
  validate_lat(lat_vec)

  na_mask <- is.na(lon_vec) | is.na(lat_vec)
  if (all(na_mask)) {
    stop("All coordinates are NA")
  }
  if (any(na_mask)) {
    warning(sprintf("%d coordinate pairs contain NA values and will be skipped",
                    sum(na_mask)))
    data <- data[!na_mask, , drop = FALSE]
    lon_vec <- lon_vec[!na_mask]
    lat_vec <- lat_vec[!na_mask]
  }

  # -------------------------------------------------------------------------
  # Perform hexification
  # -------------------------------------------------------------------------
  cell_ids <- lonlat_to_cell(lon_vec, lat_vec, hex_grid_obj)
  centers <- cell_to_lonlat(cell_ids, hex_grid_obj)

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

