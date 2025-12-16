# hexify.R
# Main user-facing convenience function
#
# This is the primary entry point for the hexify package.

#' Assign hexagonal DGGS cell IDs to geographic points
#'
#' Takes a data.frame or sf object with geographic coordinates and returns
#' the data with additional columns for hex cell ID and center coordinates.
#' By default returns a HexData object that stores the grid specification
#' for use in downstream operations.
#'
#' @param data A data.frame or sf object containing coordinates
#' @param grid A HexGrid object from \code{hex_grid()}. If provided, overrides
#'   area_km2, resolution, and aperture parameters.
#' @param lon Column name for longitude (ignored if data is sf)
#' @param lat Column name for latitude (ignored if data is sf)
#' @param area_km2 Target cell area in km² (mutually exclusive with diagonal).
#'   Alias: \code{area} for backwards compatibility.
#' @param diagonal Target cell diagonal (long diagonal) in km
#' @param resolution Grid resolution (0-30). Alternative to area_km2.
#' @param aperture Grid aperture: 3, 4, 7, or "4/3" for mixed (default 3)
#' @param mixed_aperture_level For mixed aperture "4/3": number of aperture-4 levels
#'   before switching to aperture-3 (default NULL, auto-calculated as resolution/2)
#' @param resround How to round resolution: "nearest", "up", or "down"
#'
#' @return A HexData object containing:
#'   \itemize{
#'     \item The input data with cell assignment columns
#'     \item The grid specification for downstream operations
#'   }
#'
#'   Use \code{as.data.frame(result)} to extract a plain data.frame.
#'   Use \code{as_sf(result)} to convert to sf object.
#'
#'   Added columns:
#'   \item{cell_id}{Stable DGGS cell identifier}
#'   \item{cell_cen_lon}{Longitude of cell center in degrees}
#'   \item{cell_cen_lat}{Latitude of cell center in degrees}
#'   \item{cell_area_km2}{Actual cell area in km²}
#'   \item{cell_diag_km}{Actual cell long diagonal in km}
#'
#' @details
#' For sf objects, coordinates are automatically extracted and transformed to
#' WGS84 (EPSG:4326) if needed. The geometry column is preserved.
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
#' \dontrun{
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
#' }
hexify <- function(data,
                   grid = NULL,
                   lon = "lon",
                   lat = "lat",
                   area_km2 = NULL,
                   diagonal = NULL,
                   resolution = NULL,
                   aperture = 3L,
                   mixed_aperture_level = NULL,
                   resround = "nearest") {

  # -------------------------------------------------------------------------
  # Extract or build grid specification
  # -------------------------------------------------------------------------
  if (!is.null(grid)) {
    # Grid object provided - extract parameters
    if (is_hex_grid(grid)) {
      hex_grid_obj <- grid
      aperture_num <- grid@aperture
      is_mixed_aperture <- grid@mixed_aperture
      mixed_aperture_level <- grid@mixed_aperture_level
      res <- grid@resolution
    } else if (inherits(grid, "hexify_grid")) {
      # Legacy S3 grid object
      hex_grid_obj <- hexify_grid_to_HexGrid(grid)
      aperture_num <- grid$aperture
      is_mixed_aperture <- FALSE
      res <- grid$resolution
    } else {
      stop("grid must be a HexGrid object from hex_grid() or legacy hexify_grid")
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

    # Parse aperture
    is_mixed_aperture <- FALSE
    if (is.character(aperture)) {
      if (aperture == "4/3") {
        is_mixed_aperture <- TRUE
        aperture_num <- 3L
      } else {
        stop("Aperture must be 3, 4, 7, or '4/3' (mixed aperture)")
      }
    } else {
      aperture_num <- as.integer(aperture)
      if (!aperture_num %in% c(3L, 4L, 7L)) {
        stop("Aperture must be 3, 4, or 7")
      }
    }

    # Create HexGrid object
    hex_grid_obj <- hex_grid(
      area_km2 = area_km2,
      resolution = resolution,
      aperture = if (is_mixed_aperture) "4/3" else aperture_num,
      resround = resround,
      mixed_aperture_level = mixed_aperture_level
    )
    res <- hex_grid_obj@resolution
    if (is_mixed_aperture) {
      mixed_aperture_level <- hex_grid_obj@mixed_aperture_level
    }
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
  if (is_mixed_aperture) {
    cell_ids <- cpp_lonlat_to_cell_ap43(lon_vec, lat_vec, res, mixed_aperture_level)
    centers <- cpp_cell_to_lonlat_ap43(cell_ids, res, mixed_aperture_level)
    n_cells <- 10 * (4^mixed_aperture_level) * (3^(res - mixed_aperture_level)) + 2
  } else {
    cell_ids <- cpp_lonlat_to_cell(lon_vec, lat_vec, res, aperture_num)
    centers <- cpp_cell_to_lonlat(cell_ids, res, aperture_num)
    n_cells <- 10 * (aperture_num^res) + 2
  }

  actual_area <- EARTH_SURFACE_KM2 / n_cells
  actual_spacing <- sqrt(actual_area * 2 / sqrt(3))

  # -------------------------------------------------------------------------
  # Add columns to data
  # -------------------------------------------------------------------------
  data$cell_id <- cell_ids
  data$cell_cen_lon <- centers$lon_deg
  data$cell_cen_lat <- centers$lat_deg
  data$cell_area_km2 <- actual_area
  data$cell_diag_km <- actual_spacing

  # -------------------------------------------------------------------------
  # Return HexData object
  # -------------------------------------------------------------------------
  new_hex_data(
    data = data,
    grid = hex_grid_obj,
    mapping = mapping,
    kind = "points",
    meta = list()
  )
}

#' Extract plain data frame from hexify result
#'
#' Convenience function to get a plain data.frame from hexify() result.
#' Alias for \code{as.data.frame()}.
#'
#' @param x Result from \code{hexify()}
#' @return A data.frame
#'
#' @export
#' @examples
#' \dontrun{
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#' df_plain <- hexify_df(result)
#' }
hexify_df <- function(x) {
  if (is_hex_data(x)) {
    as.data.frame(x)
  } else if (inherits(x, "sf")) {
    as.data.frame(sf::st_drop_geometry(x))
  } else {
    as.data.frame(x)
  }
}
