# constructors.R
# Constructor functions for HexGridInfo and HexData S4 classes
#
# These functions provide user-friendly interfaces for creating
# grid specifications and hexified data objects.

# =============================================================================
# HexGridInfo CONSTRUCTOR
# =============================================================================

#' Create a Hexagonal Grid Specification
#'
#' Creates a HexGridInfo object that stores all parameters needed for hexagonal
#' grid operations. Use this to define the grid once and pass it to all
#' downstream functions.
#'
#' @param area_km2 Target cell area in square kilometers. Mutually exclusive
#'   with \code{resolution}.
#' @param resolution Grid resolution level (0-30 for ISEA, 0-15 for H3).
#'   Mutually exclusive with \code{area_km2}. For H3, typical use cases by
#'   resolution:
#'   \itemize{
#'     \item 0-3: continental/country scale
#'     \item 4-7: regional/city scale
#'     \item 8-10: neighborhood/block scale (FCC uses 8-9)
#'     \item 11-15: building/sub-meter scale
#'   }
#' @param aperture Grid aperture: 3 (default), 4, 7, or "4/3" for mixed.
#'   Ignored for H3 grids (fixed at 7).
#' @param type Grid type: "isea" (default) or "h3".
#' @param resround Resolution rounding when using \code{area_km2}:
#'   "nearest" (default), "up", or "down".
#' @param crs Coordinate reference system EPSG code (default 4326 = 'WGS84').
#'
#' @return A HexGridInfo object containing the grid specification.
#'
#' @details
#' Exactly one of \code{area_km2} or \code{resolution} must be provided.
#'
#' When \code{area_km2} is provided, the resolution is calculated automatically
#' using the cell count formula: N = 10 * aperture^res + 2 (ISEA) or by
#' matching the closest H3 resolution.
#'
#' H3 grids use the Uber H3 hierarchical hexagonal system. Unlike ISEA grids,
#' H3 cells are NOT exactly equal-area (area varies by ~3-5\% depending on
#' location).
#'
#' @seealso \code{\link{hexify}} for assigning points to cells,
#'   \code{\link{HexGridInfo-class}} for class documentation
#'
#' @section One Grid, Many Datasets:
#'
#' A HexGridInfo acts as a shared spatial reference system - like a CRS, but
#' discrete and equal-area. Define the grid once, then attach multiple
#' datasets without repeating parameters:
#'
#' \preformatted{
#' # Step 1: Define the grid once
#' grid <- hex_grid(area_km2 = 1000)
#'
#' # Step 2: Attach multiple datasets to the same grid
#' birds <- hexify(bird_obs, lon = "longitude", lat = "latitude", grid = grid)
#' mammals <- hexify(mammal_obs, lon = "lon", lat = "lat", grid = grid)
#' climate <- hexify(weather_stations, lon = "x", lat = "y", grid = grid)
#'
#' # No aperture, resolution, or area needed after step 1 - the grid
#' # travels with the data.
#'
#' # Step 3: Work at the cell level
#' # Once hexified, lon/lat no longer matter - cell_id is the shared key
#' bird_counts <- aggregate(species ~ cell_id, data = as.data.frame(birds), length)
#' mammal_richness <- aggregate(species ~ cell_id, data = as.data.frame(mammals),
#'                              function(x) length(unique(x)))
#'
#' # Join datasets by cell_id - guaranteed to align because same grid
#' combined <- merge(bird_counts, mammal_richness, by = "cell_id")
#'
#' # Step 4: Visual confirmation
#' # All datasets produce identical grid overlays
#' plot(birds)   # See the grid
#' plot(mammals) # Same grid, different data
#' }
#'
#' @export
#' @examples
#' # Create grid by target area
#' grid <- hex_grid(area_km2 = 1000)
#' print(grid)
#'
#' # Create grid by resolution
#' grid <- hex_grid(resolution = 8, aperture = 3)
#'
#' # Create grid with different aperture
#' grid4 <- hex_grid(area_km2 = 500, aperture = 4)
#'
#' # Create mixed aperture grid
#' grid43 <- hex_grid(area_km2 = 1000, aperture = "4/3")
#'
#' # Use grid in hexify
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
hex_grid <- function(area_km2 = NULL,
                     resolution = NULL,
                     aperture = 3,
                     type = c("isea", "h3"),
                     resround = "nearest",
                     crs = 4326L) {

  type <- match.arg(type)

  # =========================================================================
  # H3 grid path
  # =========================================================================
  if (type == "h3") {
    if (!missing(aperture) && aperture != 3) {
      warning("aperture is ignored for H3 grids (H3 uses fixed aperture 7)")
    }

    # Validate: exactly one of area_km2 or resolution
    if (is.null(area_km2) && is.null(resolution)) {
      stop("Exactly one of 'area_km2' or 'resolution' must be provided")
    }
    if (!is.null(area_km2) && !is.null(resolution)) {
      stop("Provide either 'area_km2' or 'resolution', not both")
    }

    if (!is.null(area_km2)) {
      if (!is.numeric(area_km2) || area_km2 <= 0) {
        stop("area_km2 must be a positive number")
      }
      # Find closest H3 resolution by area
      resolution <- closest_h3_resolution(area_km2)
      warning(sprintf(
        "H3 cells are not exactly equal-area. Closest resolution %d has average area ~%.3f km^2 (requested %.3f km^2)",
        resolution, H3_AVG_AREA_KM2[resolution + 1L], area_km2
      ))
    } else {
      resolution <- as.integer(resolution)
      if (resolution < H3_MIN_RESOLUTION || resolution > H3_MAX_RESOLUTION) {
        stop(sprintf("H3 resolution must be between %d and %d",
                     H3_MIN_RESOLUTION, H3_MAX_RESOLUTION))
      }
      rlang::inform(
        "H3 cells are not exactly equal-area; area varies ~3-5% by latitude.",
        .frequency = "once",
        .frequency_id = "hexify_h3_not_equal_area"
      )
    }

    actual_area <- H3_AVG_AREA_KM2[resolution + 1L]
    actual_diagonal <- sqrt(actual_area * 2 / sqrt(3))

    grid <- new("HexGridInfo",
                aperture = "7",
                resolution = as.integer(resolution),
                area_km2 = as.numeric(actual_area),
                diagonal_km = as.numeric(actual_diagonal),
                crs = as.integer(crs),
                grid_type = "h3")
    return(grid)
  }

  # =========================================================================
  # ISEA grid path (original logic)
  # =========================================================================

  # -------------------------------------------------------------------------
  # Parse aperture (handle "4/3" mixed aperture)
  # -------------------------------------------------------------------------
  aperture_str <- as.character(aperture)

  if (aperture_str == "4/3") {
    aperture_num <- 3L  # Base aperture for resolution calculation
  } else if (aperture_str %in% c("3", "4", "7")) {
    aperture_num <- as.integer(aperture_str)
  } else {
    stop("Aperture must be 3, 4, 7, or '4/3' for mixed aperture")
  }

  # -------------------------------------------------------------------------
  # Validate area_km2 / resolution (exactly one required)
  # -------------------------------------------------------------------------
  if (is.null(area_km2) && is.null(resolution)) {
    stop("Exactly one of 'area_km2' or 'resolution' must be provided")
  }
  if (!is.null(area_km2) && !is.null(resolution)) {
    stop("Provide either 'area_km2' or 'resolution', not both")
  }

  # -------------------------------------------------------------------------
  # Calculate resolution from area if needed
  # -------------------------------------------------------------------------
  if (!is.null(area_km2)) {
    if (!is.numeric(area_km2) || area_km2 <= 0) {
      stop("area_km2 must be a positive number")
    }

    # Cell count formula: N = 10 * aperture^res + 2
    # Solving for res: res = log((N - 2) / 10) / log(aperture)
    # where N = EARTH_SURFACE_KM2 / area_km2
    n_cells <- EARTH_SURFACE_KM2 / area_km2
    res_exact <- log((n_cells - 2) / 10) / log(aperture_num)

    # Apply rounding
    resolution <- switch(resround,
      "nearest" = round(res_exact),
      "up" = ceiling(res_exact),
      "down" = floor(res_exact),
      stop("resround must be 'nearest', 'up', or 'down'")
    )

    # Clamp to valid range
    resolution <- max(MIN_RESOLUTION, min(MAX_RESOLUTION, resolution))
  } else {
    resolution <- as.integer(resolution)
    if (resolution < MIN_RESOLUTION || resolution > MAX_RESOLUTION) {
      stop(sprintf("Resolution must be between %d and %d",
                   MIN_RESOLUTION, MAX_RESOLUTION))
    }
  }

  # -------------------------------------------------------------------------
  # Calculate actual area and diagonal for this resolution
  # -------------------------------------------------------------------------
  if (aperture_str == "4/3") {
    level <- as.integer(resolution / 2)
    n_cells <- 10 * (4^level) * (3^(resolution - level)) + 2
  } else {
    n_cells <- 10 * (aperture_num^resolution) + 2
  }
  actual_area <- EARTH_SURFACE_KM2 / n_cells
  actual_diagonal <- sqrt(actual_area * 2 / sqrt(3))

  # -------------------------------------------------------------------------
  # Initialize icosahedron (required for C++ functions)
  # -------------------------------------------------------------------------
  cpp_build_icosa()

  # -------------------------------------------------------------------------
  # Create and validate HexGridInfo object
  # -------------------------------------------------------------------------
  grid <- new("HexGridInfo",
              aperture = aperture_str,
              resolution = as.integer(resolution),
              area_km2 = as.numeric(actual_area),
              diagonal_km = as.numeric(actual_diagonal),
              crs = as.integer(crs),
              grid_type = "isea")

  # Validation happens automatically via setValidity
  grid
}

# =============================================================================
# HexData CONSTRUCTOR (Internal)
# =============================================================================

#' Create a HexData Object (Internal)
#'
#' Internal constructor for HexData objects. Users should use \code{hexify()}
#' instead.
#'
#' @param data Data frame or sf object (original user data, untouched)
#' @param grid HexGridInfo object
#' @param cell_id Numeric vector of cell IDs for each row
#' @param cell_center Matrix with columns lon, lat for cell centers
#'
#' @return A HexData object
#' @keywords internal
new_hex_data <- function(data,
                         grid,
                         cell_id,
                         cell_center) {

  # Validate inputs
  if (!inherits(data, "data.frame") && !inherits(data, "sf")) {
    stop("data must be a data.frame or sf object")
  }

  if (!is_hex_grid(grid)) {
    stop("grid must be a HexGridInfo object")
  }

  # Ensure cell_center is a matrix with correct column names
  if (!is.matrix(cell_center)) {
    cell_center <- as.matrix(cell_center)
  }
  if (is.null(colnames(cell_center))) {
    colnames(cell_center) <- c("lon", "lat")
  }

  # For ISEA grids, coerce to numeric; for H3, keep as character
  if (is_h3_grid(grid)) {
    cell_id <- as.character(cell_id)
  } else {
    cell_id <- as.numeric(cell_id)
  }

  new("HexData",
      data = data,
      grid = grid,
      cell_id = cell_id,
      cell_center = cell_center)
}

# =============================================================================
# TIBBLE COERCION (if tibble is available)
# =============================================================================

#' Convert HexData to tibble
#'
#' @param x A HexData object
#' @param ... Additional arguments (ignored)
#' @return A tibble
#'
#' @export
as_tibble.HexData <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for as_tibble(). ",
         "Install with: install.packages('tibble')")
  }

  df <- x@data
  if (inherits(df, "sf")) {
    df <- as.data.frame(sf::st_drop_geometry(df))
  }

  tibble::as_tibble(df)
}

# =============================================================================
# SF COERCION
# =============================================================================

#' Convert HexData to sf Object
#'
#' Converts a HexData object to an sf spatial features object. Can create
#' either point geometries (cell centers) or polygon geometries (cell boundaries).
#'
#' @param x A HexData object
#' @param geometry Type of geometry: "point" (default) or "polygon"
#' @param ... Additional arguments (ignored)
#'
#' @return An sf object
#'
#' @details
#' For point geometry, cell centers (cell_cen_lon, cell_cen_lat) are used.
#' For polygon geometry, cell boundaries are computed using the grid specification.
#'
#' @export
#' @examples
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Get sf points
#' sf_pts <- as_sf(result)
#'
#' # Get sf polygons
#' sf_poly <- as_sf(result, geometry = "polygon")
as_sf <- function(x, geometry = c("point", "polygon"), ...) {
  UseMethod("as_sf")
}

#' @export
as_sf.HexData <- function(x, geometry = c("point", "polygon"), ...) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  geometry <- match.arg(geometry)
  data <- x@data
  grid <- x@grid

  # If already sf, check if we need to change geometry type
  if (inherits(data, "sf")) {
    if (geometry == "point") {
      return(data)  # Already has geometry
    }
    # For polygon, we need to replace geometry
    data <- sf::st_drop_geometry(data)
  }

  if (geometry == "point") {
    # Point geometry from cell centers (stored in cell_center slot)
    df_with_coords <- cbind(
      data,
      cell_id = x@cell_id,
      cell_cen_lon = x@cell_center[, "lon"],
      cell_cen_lat = x@cell_center[, "lat"]
    )
    sf::st_as_sf(
      df_with_coords,
      coords = c("cell_cen_lon", "cell_cen_lat"),
      crs = grid@crs
    )

  } else {
    # Polygon geometry from cell boundaries
    unique_ids <- unique(x@cell_id)

    polys_sf <- cell_to_sf(unique_ids, grid)

    # Add cell_id to data for merge
    data_with_id <- cbind(data, cell_id = x@cell_id)

    # Merge with original data
    result <- merge(polys_sf, data_with_id, by = "cell_id", all.y = TRUE)
    sf::st_as_sf(result)
  }
}

#' @export
as_sf.default <- function(x, ...) {
  stop("as_sf() is not defined for objects of class ", class(x)[1])
}
