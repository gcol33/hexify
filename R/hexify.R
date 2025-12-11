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
#' @param area Target cell area in km² (mutually exclusive with spacing)
#' @param spacing Target cell spacing (long diagonal) in km (mutually exclusive with area)
#' @param aperture Grid aperture: 3, 4, or 7 (default 3)
#' @param resround How to round resolution: "nearest", "up", or "down"
#'
#' @return The input data with additional columns:
#'   \item{hex_id}{Stable DGGS cell identifier (character)}
#'   \item{hex_cen_lon}{Longitude of cell center in degrees}
#'   \item{hex_cen_lat}{Latitude of cell center in degrees}
#'
#' @details
#' For sf objects, coordinates are automatically extracted and transformed to
#' WGS84 (EPSG:4326) if needed. The geometry column is preserved.
#'
#' Either `area` or `spacing` must be provided, but not both:
#' - `area`: Target cell area in square kilometers
#' - `spacing`: Long diagonal of hexagon in kilometers
#'
#' The spacing (long diagonal) relates to area approximately as:
#' area ≈ (3 * sqrt(3) / 2) * (spacing / 2)² ≈ 0.6495 * spacing²
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
#' # Using spacing instead of area
#' result <- hexify(df, lon = "lon", lat = "lat", spacing = 50)
#' }
hexify <- function(data,
                   lon = "lon",
                   lat = "lat",
                   area = NULL,
                   spacing = NULL,
                   aperture = 3L,
                   resround = "nearest") {

  # Validate area/spacing

  if (is.null(area) && is.null(spacing)) {
    stop("Either 'area' (km²) or 'spacing' (long diagonal in km) must be provided")
  }
  if (!is.null(area) && !is.null(spacing)) {
    stop("Provide either 'area' or 'spacing', not both")
  }

  # Convert spacing to area if provided
  # Long diagonal = 2 * (sqrt(3)/2) * side = sqrt(3) * side
  # Area = (3 * sqrt(3) / 2) * side²
  # Therefore: area = (3 * sqrt(3) / 2) * (spacing / sqrt(3))² = spacing² / 2 * sqrt(3)
  if (!is.null(spacing)) {
    area <- spacing^2 * sqrt(3) / 2
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

  # Build grid
  grid <- hexify_construct(
    area = area,
    aperture = as.integer(aperture),
    resround = resround
  )

  # Get cell assignments as integer SEQNUM (dggridR-compatible)
  seqnum <- cpp_lonlat_to_seqnum_dggrid(lon_vec, lat_vec, grid$resolution, grid$aperture)

  # Get cell centers
  centers <- cpp_seqnum_to_lonlat_dggrid(seqnum, grid$resolution, grid$aperture)

  # Add columns to data (using integer SEQNUM like dggridR)
  data$hex_id <- as.integer(seqnum)
  data$hex_cen_lon <- centers$lon_deg
  data$hex_cen_lat <- centers$lat_deg

  return(data)
}
