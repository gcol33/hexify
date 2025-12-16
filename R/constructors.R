# constructors.R
# Constructor functions for HexGrid and HexData S4 classes
#
# These functions provide user-friendly interfaces for creating
# grid specifications and hexified data objects.

# =============================================================================
# HexGrid CONSTRUCTOR
# =============================================================================

#' Create a Hexagonal Grid Specification
#'
#' Creates a HexGrid object that stores all parameters needed for hexagonal
#' grid operations. Use this to define the grid once and pass it to all
#' downstream functions.
#'
#' @param area_km2 Target cell area in square kilometers. Mutually exclusive
#'   with \code{resolution}.
#' @param resolution Grid resolution level (0-30). Mutually exclusive with
#'   \code{area_km2}.
#' @param aperture Grid aperture: 3 (default), 4, 7, or "4/3" for mixed.
#' @param topology Grid topology: "H" (default) or "HEXAGON".
#' @param grid_system Grid system: "ISEA" (default, only supported option).
#' @param index_type Index format: "integer" (default) or "character".
#'   Automatically set to "character" for high resolutions to avoid overflow.
#' @param resround Resolution rounding when using \code{area_km2}:
#'   "nearest" (default), "up", or "down".
#' @param mixed_aperture_level For aperture "4/3": number of aperture-4 levels
#'   before switching to aperture-3. If NULL, defaults to resolution/2.
#' @param crs_input Input coordinate reference system (default 4326 = WGS84).
#' @param ... Additional metadata to store in the grid specification.
#'
#' @return A HexGrid object containing the grid specification.
#'
#' @details
#' Exactly one of \code{area_km2} or \code{resolution} must be provided.
#'
#' When \code{area_km2} is provided, the resolution is calculated automatically
#' using the cell count formula: N = 10 * aperture^res + 2.
#'
#' The \code{index_type} parameter controls how cell IDs are stored:
#' \itemize{
#'   \item "integer": Numeric cell IDs (faster, but overflow risk at high res)
#'   \item "character": String indices (safer for high resolution grids)
#' }
#'
#' For resolutions above 15 (aperture 3) or 12 (aperture 4), the function
#' automatically switches to "character" index type to prevent integer overflow.
#'
#' @seealso \code{\link{hexify}} for assigning points to cells,
#'   \code{\link{HexGrid-class}} for class documentation
#'
#' @section One Grid, Many Datasets:
#'
#' A HexGrid acts as a shared spatial reference system - like a CRS, but
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
#' \dontrun{
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
#' result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
#' }
hex_grid <- function(area_km2 = NULL,
                     resolution = NULL,
                     aperture = 3L,
                     topology = "H",
                     grid_system = "ISEA",
                     index_type = NULL,
                     resround = "nearest",
                     mixed_aperture_level = NULL,
                     crs_input = 4326L,
                     ...) {

  # -------------------------------------------------------------------------
  # Parse aperture (handle "4/3" mixed aperture)
  # -------------------------------------------------------------------------
  mixed_aperture <- FALSE

  if (is.character(aperture)) {
    if (aperture == "4/3") {
      mixed_aperture <- TRUE
      aperture_num <- 3L  # Base aperture for resolution calculation
    } else {
      stop("Aperture must be 3, 4, 7, or '4/3' for mixed aperture")
    }
  } else {
    aperture_num <- as.integer(aperture)
    if (!aperture_num %in% c(3L, 4L, 7L)) {
      stop("Aperture must be 3, 4, or 7")
    }
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

    # Calculate actual area for this resolution
    if (mixed_aperture && !is.null(mixed_aperture_level)) {
      n_cells <- 10 * (4^mixed_aperture_level) *
        (3^(resolution - mixed_aperture_level)) + 2
    } else {
      n_cells <- 10 * (aperture_num^resolution) + 2
    }
    area_km2 <- EARTH_SURFACE_KM2 / n_cells
  }

  # -------------------------------------------------------------------------
  # Handle mixed aperture level
  # -------------------------------------------------------------------------
  if (mixed_aperture) {
    if (is.null(mixed_aperture_level)) {
      mixed_aperture_level <- as.integer(resolution / 2)
    } else {
      mixed_aperture_level <- as.integer(mixed_aperture_level)
    }

    if (mixed_aperture_level < 0 || mixed_aperture_level > resolution) {
      stop("mixed_aperture_level must be between 0 and resolution")
    }
  } else {
    mixed_aperture_level <- NA_integer_
  }

  # -------------------------------------------------------------------------
  # Determine index_type (auto-select for high resolution)
  # -------------------------------------------------------------------------
  if (is.null(index_type)) {
    # Check for potential integer overflow
    # R integers max out at 2^31 - 1 ≈ 2.1 billion
    # Aperture 3: 10 * 3^15 ≈ 143 million (safe), 3^16 ≈ 430 million (safe),
    #             3^19 ≈ 11.6 billion (overflow)
    # Aperture 4: 10 * 4^15 ≈ 10.7 billion (overflow at 15)

    max_safe_res <- if (aperture_num == 3L) 18L
                    else if (aperture_num == 4L) 14L
                    else 11L  # aperture 7

    index_type <- if (resolution > max_safe_res) "character" else "integer"
  }

  if (!index_type %in% c("integer", "character")) {
    stop("index_type must be 'integer' or 'character'")
  }

  # -------------------------------------------------------------------------
  # Validate other parameters
  # -------------------------------------------------------------------------
  if (!topology %in% c("H", "HEXAGON")) {
    stop("topology must be 'H' or 'HEXAGON'")
  }

  if (grid_system != "ISEA") {
    stop("Only 'ISEA' grid_system is currently supported")
  }

  # -------------------------------------------------------------------------
  # Collect additional metadata
  # -------------------------------------------------------------------------
  extra_meta <- list(...)

  # -------------------------------------------------------------------------
  # Initialize icosahedron (required for C++ functions)
  # -------------------------------------------------------------------------
  cpp_build_icosa()

  # -------------------------------------------------------------------------
  # Create and validate HexGrid object
  # -------------------------------------------------------------------------
  grid <- new("HexGrid",
              aperture = aperture_num,
              resolution = as.integer(resolution),
              area_km2 = as.numeric(area_km2),
              grid_system = grid_system,
              topology = topology,
              index_type = index_type,
              crs_input = as.integer(crs_input),
              crs_work = as.integer(crs_input),
              mixed_aperture = mixed_aperture,
              mixed_aperture_level = mixed_aperture_level,
              meta = extra_meta)

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
#' @param data Data frame or sf object with cell assignments
#' @param grid HexGrid object
#' @param mapping List of column name mappings
#' @param kind Type of data: "points", "cells", or "unknown"
#' @param meta Additional metadata
#'
#' @return A HexData object
#' @keywords internal
new_hex_data <- function(data,
                         grid,
                         mapping = list(),
                         kind = "unknown",
                         meta = list()) {

  # Validate inputs
  if (!inherits(data, "data.frame") && !inherits(data, "sf")) {
    stop("data must be a data.frame or sf object")
  }

  if (!is_hex_grid(grid)) {
    stop("grid must be a HexGrid object")
  }

  if (!kind %in% c("points", "cells", "unknown")) {
    stop("kind must be 'points', 'cells', or 'unknown'")
  }

  new("HexData",
      data = data,
      grid = grid,
      mapping = mapping,
      kind = kind,
      meta = meta)
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
#' \dontrun{
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Get sf points
#' sf_pts <- as_sf(result)
#'
#' # Get sf polygons
#' sf_poly <- as_sf(result, geometry = "polygon")
#' }
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
    # Point geometry from cell centers
    if (!"cell_cen_lon" %in% names(data) || !"cell_cen_lat" %in% names(data)) {
      stop("Data must contain 'cell_cen_lon' and 'cell_cen_lat' columns")
    }

    sf::st_as_sf(data,
                 coords = c("cell_cen_lon", "cell_cen_lat"),
                 crs = grid@crs_input)

  } else {
    # Polygon geometry from cell boundaries
    unique_ids <- unique(data$cell_id)

    # Generate polygons
    polys_sf <- hexify_cell_to_sf(
      cell_id = unique_ids,
      resolution = grid@resolution,
      aperture = grid@aperture,
      return_sf = TRUE
    )

    # Merge with original data
    result <- merge(polys_sf, data, by = "cell_id", all.y = TRUE)
    sf::st_as_sf(result)
  }
}

#' @export
as_sf.default <- function(x, ...) {
  stop("as_sf() is not defined for objects of class ", class(x)[1])
}
