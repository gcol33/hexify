# hexify.R
# Main user-facing convenience function
#
# This is the primary entry point for the hexify package.

#' Assign hexagonal DGGS cell IDs to geographic points
#'
#' Takes a data.frame or sf object with geographic coordinates and returns
#' the same data with additional columns for hex cell ID and center coordinates.
#'
#' @param data A data.frame or sf object containing coordinates
#' @param lon Column name for longitude (ignored if data is sf)
#' @param lat Column name for latitude (ignored if data is sf)
#' @param area Target cell area in km² (mutually exclusive with diagonal)
#' @param diagonal Target cell diagonal (long diagonal) in km (mutually exclusive with area)
#' @param aperture Grid aperture: 3, 4, 7, or "4/3" for mixed (default 3)
#' @param mixed_aperture_level For mixed aperture "4/3": number of aperture-4 levels
#'   before switching to aperture-3 (default NULL, auto-calculated as resolution/2)
#' @param resround How to round resolution: "nearest", "up", or "down"
#'
#' @return The input data with additional columns:
#'   \item{cell_id}{Stable DGGS cell identifier (integer)}
#'   \item{cell_cen_lon}{Longitude of cell center in degrees}
#'   \item{cell_cen_lat}{Latitude of cell center in degrees}
#'   \item{cell_area}{Actual cell area in km² (based on matched resolution)}
#'   \item{cell_diag}{Actual cell long diagonal in km}
#'
#' @details
#' For sf objects, coordinates are automatically extracted and transformed to
#' WGS84 (EPSG:4326) if needed. The geometry column is preserved.
#'
#' Either `area` or `diagonal` must be provided, but not both:
#' - `area`: Target cell area in square kilometers
#' - `diagonal`: Long diagonal of hexagon in kilometers
#'
#' The diagonal relates to area approximately as:
#' area ≈ (3 * sqrt(3) / 2) * (diagonal / 2)² ≈ 0.6495 * diagonal²
#'
#' Supported apertures:
#' - `aperture = 3`: ISEA3H (default, compatible with dggridR).
#'   Cell count = 10 * 3^res + 2.
#' - `aperture = 4`: ISEA4H. Cell count = 10 * 4^res + 2.
#' - `aperture = 7`: ISEA7H. Cell count = 10 * 7^res + 2.
#' - `aperture = "4/3"`: ISEA43H mixed aperture.
#'   Uses aperture-4 for first `mixed_aperture_level` levels, then aperture-3.
#'   Cell count = 10 * 4^mixed_level * 3^(res - mixed_level) + 2.
#'
#' @export
#' @examples
#' \dontrun{
#' # Simple data.frame
#' df <- data.frame(
#'   site = c("Vienna", "Paris", "Madrid"),
#'   lon = c(16.37, 2.35, -3.70),
#'   lat = c(48.21, 48.86, 40.42)
#' )
#' result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
#'
#' # With sf object (any CRS)
#' library(sf)
#' pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
#' result_sf <- hexify(pts, area = 1000)
#'
#' # Using diagonal instead of area
#' result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)
#'
#' # Different apertures
#' result_ap4 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)
#' result_ap7 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)
#'
#' # Mixed aperture (ISEA43H)
#' result_mixed <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "4/3")
#' }
hexify <- function(data,
                   lon = "lon",
                   lat = "lat",
                   area = NULL,
                   diagonal = NULL,
                   aperture = 3L,
                   mixed_aperture_level = NULL,
                   resround = "nearest") {

  # Validate area/diagonal

  if (is.null(area) && is.null(diagonal)) {
    stop("Either 'area' (km²) or 'diagonal' (long diagonal in km) must be provided")
  }
  if (!is.null(area) && !is.null(diagonal)) {
    stop("Provide either 'area' or 'diagonal', not both")
  }

  # Parse aperture - can be integer (3, 4, 7) or string ("4/3")
  is_mixed_aperture <- FALSE
  if (is.character(aperture)) {
    # Parse mixed aperture string
    if (aperture == "4/3") {
      is_mixed_aperture <- TRUE
      # "4/3" means "aperture 4 first, then 3"
      aperture_num <- 3L  # Base aperture for resolution calculation
    } else {
      stop("Aperture must be 3, 4, 7, or '4/3' (mixed aperture)")
    }
  } else {
    aperture_num <- as.integer(aperture)
    if (!aperture_num %in% c(3L, 4L, 7L)) {
      stop("Aperture must be 3, 4, or 7")
    }
  }

  # Convert diagonal to area if provided
  # Long diagonal = 2 * (sqrt(3)/2) * side = sqrt(3) * side
  # Area = (3 * sqrt(3) / 2) * side²
  # Therefore: area = (3 * sqrt(3) / 2) * (diagonal / sqrt(3))² = diagonal² / 2 * sqrt(3)
  if (!is.null(diagonal)) {
    area <- diagonal^2 * sqrt(3) / 2
  }

  # Handle sf objects
  is_sf <- inherits(data, "sf")

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

    # Extract coordinates
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

  # Check for NA values
  na_mask <- is.na(lon_vec) | is.na(lat_vec)
  if (all(na_mask)) {
    stop("All coordinates are NA")
  }
  if (any(na_mask)) {
    warning(sprintf("%d coordinate pairs contain NA values and will be skipped",
                    sum(na_mask)))
  }

  # Build grid
  grid <- hexify_grid(
    area = area,
    aperture = aperture_num,
    resround = resround
  )

  # Handle mixed aperture
  if (is_mixed_aperture) {
    # Default mixed_aperture_level to half the resolution if not specified
    if (is.null(mixed_aperture_level)) {
      mixed_aperture_level <- as.integer(grid$resolution / 2)
    }
    mixed_aperture_level <- as.integer(mixed_aperture_level)

    if (mixed_aperture_level < 0 || mixed_aperture_level > grid$resolution) {
      stop("mixed_aperture_level must be between 0 and resolution")
    }

    # Get cell assignments using mixed aperture function
    cell_ids <- cpp_lonlat_to_cell_ap43(lon_vec, lat_vec,
                                         grid$resolution, mixed_aperture_level)

    # Get cell centers
    centers <- cpp_cell_to_lonlat_ap43(cell_ids, grid$resolution, mixed_aperture_level)

    # Calculate actual cell count for mixed aperture
    # N = 10 * 4^mixed_level * 3^(res - mixed_level) + 2
    n_cells <- 10 * (4^mixed_aperture_level) * (3^(grid$resolution - mixed_aperture_level)) + 2
  } else {
    # Get cell assignments as integer cell ID
    cell_ids <- cpp_lonlat_to_cell(lon_vec, lat_vec, grid$resolution, grid$aperture)

    # Get cell centers
    centers <- cpp_cell_to_lonlat(cell_ids, grid$resolution, grid$aperture)

    # Calculate cell count: N = 10 * aperture^res + 2
    n_cells <- 10 * (grid$aperture^grid$resolution) + 2
  }

  # Calculate actual cell area and spacing from the grid
  # Uses EARTH_SURFACE_KM2 from constants.R (510072000 km²)
  actual_area <- EARTH_SURFACE_KM2 / n_cells
  # Long diagonal = sqrt(2 * area / sqrt(3)) * sqrt(2) = sqrt(4 * area / sqrt(3))
  actual_spacing <- sqrt(actual_area * 2 / sqrt(3))

  # Add columns to data (cell ID as numeric to handle high-resolution grids)
  # At resolution 20 with aperture 3: N = 10 * 3^20 + 2 ≈ 34 billion (exceeds R integer max)
  data$cell_id <- cell_ids
  data$cell_cen_lon <- centers$lon_deg
  data$cell_cen_lat <- centers$lat_deg
  data$cell_area <- actual_area
  data$cell_diag <- actual_spacing

  return(data)
}
