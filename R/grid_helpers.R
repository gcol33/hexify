# grid_helpers.R
# Helper functions that accept HexGridInfo or HexData objects
#
# These functions wrap the low-level coordinate conversion functions
# to accept grid specifications, eliminating the need to repeat
# aperture/resolution parameters.

# =============================================================================
# COORDINATE CONVERSION HELPERS
# =============================================================================

#' Convert longitude/latitude to cell ID
#'
#' Converts geographic coordinates to DGGS cell IDs using a grid specification.
#'
#' @param lon Numeric vector of longitudes in degrees
#' @param lat Numeric vector of latitudes in degrees
#' @param grid A HexGridInfo or HexData object, or legacy hexify_grid
#'
#' @return Numeric vector of cell IDs
#'
#' @details
#' This function accepts either a HexGridInfo object from \code{hex_grid()} or
#' a HexData object from \code{hexify()}. If a HexData object is provided,
#' its grid specification is extracted automatically.
#'
#' @seealso \code{\link{cell_to_lonlat}} for the inverse operation,
#'   \code{\link{hex_grid}} for creating grid specifications
#'
#' @export
#' @examples
#' grid <- hex_grid(area_km2 = 1000)
#' cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
#'
#' # Or use HexData object
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#' cells <- lonlat_to_cell(lon = 5, lat = 48, grid = result)
lonlat_to_cell <- function(lon, lat, grid) {
  g <- extract_grid(grid)

  if (is_h3_grid(g)) {
    return(cpp_h3_latLngToCell(as.numeric(lon), as.numeric(lat), g@resolution))
  }

  if (g@aperture == "4/3") {
    level <- as.integer(g@resolution / 2)
    cpp_lonlat_to_cell_ap43(
      as.numeric(lon),
      as.numeric(lat),
      g@resolution,
      level
    )
  } else {
    cpp_lonlat_to_cell(
      as.numeric(lon),
      as.numeric(lat),
      g@resolution,
      as.integer(g@aperture)
    )
  }
}

#' Convert cell ID to longitude/latitude
#'
#' Converts DGGS cell IDs back to geographic coordinates (cell centers).
#'
#' @param cell_id Numeric vector of cell IDs
#' @param grid A HexGridInfo or HexData object
#'
#' @return Data frame with lon_deg and lat_deg columns
#'
#' @seealso \code{\link{lonlat_to_cell}} for the forward operation
#'
#' @export
#' @examples
#' grid <- hex_grid(area_km2 = 1000)
#' cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
#' coords <- cell_to_lonlat(cells, grid)
cell_to_lonlat <- function(cell_id, grid) {
  g <- extract_grid(grid)

  if (is_h3_grid(g)) {
    result <- cpp_h3_cellToLatLng(as.character(cell_id))
    return(data.frame(lon_deg = result$lon, lat_deg = result$lat))
  }

  if (g@aperture == "4/3") {
    level <- as.integer(g@resolution / 2)
    cpp_cell_to_lonlat_ap43(
      as.numeric(cell_id),
      g@resolution,
      level
    )
  } else {
    cpp_cell_to_lonlat(
      as.numeric(cell_id),
      g@resolution,
      as.integer(g@aperture)
    )
  }
}

#' Convert cell IDs to sf polygons
#'
#' Creates sf polygon geometries for hexagonal grid cells.
#'
#' @param cell_id Numeric vector of cell IDs. If NULL and x is HexData,
#'   uses cells from x.
#' @param grid A HexGridInfo or HexData object. If HexData and cell_id is NULL,
#'   polygons are generated for all cells in the data.
#' @param wrap_dateline Logical. If TRUE (default), calls
#'   \code{sf::st_wrap_dateline()} to split antimeridian-crossing polygons.
#'   Set to FALSE for orthographic/globe projections where wrapping creates gaps.
#'
#' @return sf object with cell_id and geometry columns
#'
#' @details
#' When called with a HexData object and no cell_id argument, this function
#' generates polygons for all unique cells in the data, which is useful for
#' plotting.
#'
#' @seealso \code{\link{hex_grid}} for grid specifications,
#'   \code{\link{as_sf}} for converting HexData to sf
#'
#' @export
#' @examples
#' # From grid specification
#' grid <- hex_grid(area_km2 = 1000)
#' cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
#' polys <- cell_to_sf(cells, grid)
#'
#' # From HexData (all cells)
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#' polys <- cell_to_sf(grid = result)
cell_to_sf <- function(cell_id = NULL, grid, wrap_dateline = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Handle HexData input
  if (is_hex_data(grid)) {
    if (is.null(cell_id)) {
      cell_id <- unique(grid@cell_id)
    }
    g <- grid@grid
  } else {
    g <- extract_grid(grid)
    if (is.null(cell_id)) {
      stop("cell_id required when grid is not HexData")
    }
  }

  # Remove NA and duplicates
  cell_id <- unique(cell_id[!is.na(cell_id)])
  if (length(cell_id) == 0) {
    stop("No valid cell_id values")
  }

  # H3 path: use native C backend for boundaries
  if (is_h3_grid(g)) {
    boundaries <- cpp_h3_cellToBoundary(as.character(cell_id))
    polygons <- lapply(boundaries, function(coords) {
      if (nrow(coords) == 0) return(sf::st_polygon())
      lons <- coords[, 1]
      lon_range <- max(lons, na.rm = TRUE) - min(lons, na.rm = TRUE)
      if (lon_range > 180) {
        lons[lons < 0] <- lons[lons < 0] + 360
        coords[, 1] <- lons
        mean_lon <- mean(lons)
        if (mean_lon > 180) {
          coords[, 1] <- coords[, 1] - 360
        }
      }
      sf::st_polygon(list(coords))
    })
    sfc <- sf::st_sfc(polygons, crs = g@crs)
    sfc <- suppressWarnings(sf::st_make_valid(sfc))
    result_sf <- sf::st_sf(cell_id = as.character(cell_id), geometry = sfc)
    if (wrap_dateline) {
      result_sf <- sf::st_wrap_dateline(result_sf,
        options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
    }
    return(result_sf)
  }

  # ISEA path: generate polygons using C++ function
  # Convert aperture to integer for C++ (mixed aperture "4/3" uses 3)
  aperture_int <- if (g@aperture == "4/3") 3L else as.integer(g@aperture)

  corners_list <- cpp_cell_to_corners(
    as.numeric(cell_id),
    g@resolution,
    aperture_int
  )


  # Handle antimeridian-crossing polygons: normalize coordinates so each
  # polygon is contiguous. When wrap_dateline = TRUE, st_wrap_dateline
  # then splits at ±180° for correct flat-map rendering. For globe/
  # orthographic projections, pass wrap_dateline = FALSE to keep cells intact.

  polygons <- lapply(corners_list, function(coords) {
    lons <- coords[, 1]
    lon_range <- max(lons, na.rm = TRUE) - min(lons, na.rm = TRUE)

    if (lon_range > 180) {
      # Polygon crosses antimeridian - normalize to be contiguous
      # Shift negative lons to 0-360 range
      lons[lons < 0] <- lons[lons < 0] + 360
      coords[, 1] <- lons

      # Now shift back to standard range, but keeping contiguity
      # If centroid is > 180, shift everything by -360
      mean_lon <- mean(lons)
      if (mean_lon > 180) {
        coords[, 1] <- coords[, 1] - 360
      }
    }

    sf::st_polygon(list(coords))
  })

  sfc <- sf::st_sfc(polygons, crs = g@crs)

  # Fix any invalid geometries (self-intersecting polygons, etc.)
  # suppressWarnings: antimeridian normalization may temporarily produce
  # out-of-range longitudes that st_wrap_dateline corrects below
  sfc <- suppressWarnings(sf::st_make_valid(sfc))

  result_sf <- sf::st_sf(cell_id = cell_id, geometry = sfc)
  if (wrap_dateline) {
    result_sf <- sf::st_wrap_dateline(result_sf,
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
  }
  result_sf
}

# =============================================================================
# GRID GENERATION HELPERS
# =============================================================================

#' Generate a rectangular grid of hexagons
#'
#' Creates hexagon polygons covering a rectangular geographic region.
#' For H3 grids, all cells that overlap the bounding box are included
#' (not just cells whose center falls inside), ensuring full spatial coverage.
#'
#' @param bbox Bounding box as c(xmin, ymin, xmax, ymax), or an sf/sfc object
#' @param grid A HexGridInfo object specifying the grid parameters
#'
#' @return sf object with hexagon polygons
#'
#' @seealso \code{\link{grid_global}} for global grids
#'
#' @export
#' @examples
#' grid <- hex_grid(area_km2 = 5000)
#' europe <- grid_rect(c(-10, 35, 30, 60), grid)
#' plot(europe)
grid_rect <- function(bbox, grid) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required")
  }

  g <- extract_grid(grid)

  # Handle sf/sfc bbox input
  if (inherits(bbox, c("sf", "sfc", "bbox"))) {
    bbox <- as.numeric(sf::st_bbox(bbox))
  }

  # H3 path: fill bbox with cells using native C backend
  if (is_h3_grid(g)) {
    bbox_coords <- matrix(c(
      bbox[1], bbox[2],
      bbox[3], bbox[2],
      bbox[3], bbox[4],
      bbox[1], bbox[4],
      bbox[1], bbox[2]
    ), ncol = 2, byrow = TRUE)
    cell_ids <- cpp_h3_polygonToCells(bbox_coords, g@resolution)
    if (length(cell_ids) == 0) {
      stop("No H3 cells found in the specified bounding box at resolution ", g@resolution)
    }
    return(cell_to_sf(cell_ids, g))
  }

  minlon <- bbox[1]
  minlat <- bbox[2]
  maxlon <- bbox[3]
  maxlat <- bbox[4]

  # Create sampling grid - use diagonal_km from grid if available
  diagonal <- if (!is.na(g@diagonal_km)) g@diagonal_km else sqrt(g@area_km2 * 2 / sqrt(3))
  spacing_deg <- diagonal / KM_PER_DEGREE * 0.8

  lons <- seq(minlon, maxlon, by = spacing_deg)
  lats <- seq(minlat, maxlat, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  # Get unique cells covering the region
  cell_ids <- lonlat_to_cell(grid_pts$lon, grid_pts$lat, g)
  unique_cells <- unique(cell_ids)

  cell_to_sf(unique_cells, g)
}

#' Generate a global hexagon grid
#'
#' Creates hexagon polygons covering the entire Earth.
#'
#' @param grid A HexGridInfo object specifying the grid parameters
#' @param wrap_dateline Logical. If TRUE (default), antimeridian-crossing
#'   polygons are split at +/-180 degrees. Set to FALSE for orthographic/globe
#'   projections where wrapping creates gaps.
#'
#' @return sf object with hexagon polygons
#'
#' @details
#' This function generates a complete global grid by sampling points
#' densely across the globe. For large grids (many small cells),
#' consider using \code{grid_rect()} to generate regional subsets.
#'
#' @seealso \code{\link{grid_rect}} for regional grids
#'
#' @export
#' @examples
#' # Coarse global grid
#' grid <- hex_grid(area_km2 = 100000)
#' global <- grid_global(grid)
#' plot(global)
grid_global <- function(grid, wrap_dateline = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required")
  }

  g <- extract_grid(grid)

  # H3 path: fill globe using native C backend
  if (is_h3_grid(g)) {
    h3_n_cells <- 2 + 120 * 7^g@resolution
    if (h3_n_cells > 2e6) {
      warning(sprintf(
        "H3 global grid at res %d has ~%.0f cells. This may take a while.",
        g@resolution, h3_n_cells
      ))
    }
    # Split globe into quadrants for polygonToCells
    quads <- list(
      matrix(c(-180, 0, 0, 0, 0, 90, -180, 90, -180, 0), ncol = 2, byrow = TRUE),
      matrix(c(0, 0, 180, 0, 180, 90, 0, 90, 0, 0), ncol = 2, byrow = TRUE),
      matrix(c(-180, -90, 0, -90, 0, 0, -180, 0, -180, -90), ncol = 2, byrow = TRUE),
      matrix(c(0, -90, 180, -90, 180, 0, 0, 0, 0, -90), ncol = 2, byrow = TRUE)
    )
    all_cells <- character(0)
    for (q in quads) {
      quad_cells <- cpp_h3_polygonToCells(q, g@resolution)
      all_cells <- c(all_cells, quad_cells)
    }
    cell_ids <- unique(all_cells)
    return(cell_to_sf(cell_ids, g, wrap_dateline = wrap_dateline))
  }

  # Estimate cell count for warning (ISEA)
  if (g@aperture == "4/3") {
    level <- as.integer(g@resolution / 2)
    n_cells <- 10 * (4^level) * (3^(g@resolution - level)) + 2
  } else {
    ap <- as.integer(g@aperture)
    n_cells <- 10 * (ap^g@resolution) + 2
  }
  if (n_cells > 100000) {
    warning(sprintf(
      "This will generate approximately %.0f cells. Consider larger area_km2.",
      n_cells
    ))
  }

  # Dense sampling - use diagonal_km from grid if available
  diagonal <- if (!is.na(g@diagonal_km)) g@diagonal_km else sqrt(g@area_km2 * 2 / sqrt(3))
  spacing_deg <- diagonal / KM_PER_DEGREE * 0.7

  lons <- seq(-180, 180, by = spacing_deg)
  lats <- seq(-85, 85, by = spacing_deg)
  grid_pts <- expand.grid(lon = lons, lat = lats)

  # Add polar cap sampling (±85 to ±90 degrees)
  # The regular grid misses polar cells because lat stops at ±85
  # Near poles, longitude spacing must be DENSER not coarser - at 89°N,
  # the entire circumference is only ~6.3° of longitude-equivalent distance.
  # Use the same spacing_deg (or denser) to ensure we catch all cells.
  polar_lon_spacing <- min(spacing_deg, 15)  # At most 15°, or cell-based spacing
  polar_lons <- seq(-180, 180, by = polar_lon_spacing)
  polar_lats <- c(seq(85.5, 89.99, by = 0.5), seq(-89.99, -85.5, by = 0.5))
  polar_pts <- expand.grid(lon = polar_lons, lat = polar_lats)

  # Combine main grid with polar samples
  grid_pts <- rbind(grid_pts, polar_pts)

  cell_ids <- lonlat_to_cell(grid_pts$lon, grid_pts$lat, g)
  unique_cells <- unique(cell_ids)

  cell_to_sf(unique_cells, g, wrap_dateline = wrap_dateline)
}

#' Clip hexagon grid to polygon boundary
#'
#' Creates hexagon polygons clipped to a given polygon boundary. This is useful
#' for generating grids that conform to country borders, study areas, or other
#' irregular boundaries.
#'
#' @param boundary An sf/sfc polygon to clip to. Can be a country boundary,
#'   study area, or any polygon geometry.
#' @param grid A HexGridInfo object specifying the grid parameters
#' @param crop If TRUE (default), cells are cropped to the boundary. If FALSE,
#'   only cells whose centroids fall within the boundary are kept (no cropping).
#'
#' @return sf object with hexagon polygons clipped to the boundary
#'
#' @details
#' The function first generates cells covering the boundary polygon, then
#' clips or filters them. For H3 grids, all cells that overlap the boundary
#' are included (not just cells whose center falls inside), ensuring full
#' spatial coverage with no gaps along the boundary edge.
#'
#' When \code{crop = TRUE}, hexagons are geometrically intersected with the
#' boundary, which may produce partial hexagons at the edges. When
#' \code{crop = FALSE}, only complete hexagons whose centroids fall within
#' the boundary are returned.
#'
#' @seealso \code{\link{grid_rect}} for rectangular grids,
#'   \code{\link{grid_global}} for global grids
#'
#' @export
#' @examples
#' \donttest{
#' # Get France boundary from built-in world map
#' france <- hexify_world[hexify_world$name == "France", ]
#'
#' # Create grid clipped to France
#' grid <- hex_grid(area_km2 = 2000)
#' france_grid <- grid_clip(france, grid)
#'
#' # Plot result
#' library(ggplot2)
#' ggplot() +
#'   geom_sf(data = france, fill = "gray95") +
#'   geom_sf(data = france_grid, fill = alpha("steelblue", 0.3),
#'           color = "steelblue") +
#'   theme_minimal()
#'
#' # Keep only complete hexagons (no cropping)
#' france_grid_complete <- grid_clip(france, grid, crop = FALSE)
#' }
grid_clip <- function(boundary, grid, crop = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required")
  }

  # Validate boundary
  if (!inherits(boundary, c("sf", "sfc"))) {
    stop("boundary must be an sf or sfc object")
  }

  g <- extract_grid(grid)

  # H3 path: fill boundary polygon using native C backend
  if (is_h3_grid(g)) {
    # Disable S2 for spatial operations
    s2_state <- sf::sf_use_s2()
    sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(s2_state), add = TRUE)

    boundary_geom <- sf::st_geometry(boundary)
    if (length(boundary_geom) > 1) {
      boundary_geom <- sf::st_union(boundary_geom)
    }
    boundary_geom <- sf::st_make_valid(boundary_geom)

    # Extract polygon rings for cpp_h3_polygonToCells
    polys <- sf::st_cast(boundary_geom, "POLYGON")
    all_cells <- character(0)
    for (p in polys) {
      rings <- unclass(p)
      outer_ring <- rings[[1]]
      hole_rings <- if (length(rings) > 1) rings[-1] else NULL
      pcells <- cpp_h3_polygonToCells(outer_ring, g@resolution, holes = hole_rings)
      all_cells <- c(all_cells, pcells)
    }
    cell_ids <- unique(all_cells)
    hex_sf <- cell_to_sf(cell_ids, g)

    if (crop) {
      hex_sf <- sf::st_make_valid(hex_sf)
      result <- tryCatch({
        suppressWarnings(sf::st_intersection(hex_sf, boundary_geom))
      }, error = function(e) {
        boundary_buf <- sf::st_buffer(boundary_geom, 0)
        hex_buf <- sf::st_buffer(hex_sf, 0)
        suppressWarnings(sf::st_intersection(hex_buf, boundary_buf))
      })
      geom_types <- sf::st_geometry_type(result)
      result <- result[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
      return(result)
    }
    return(hex_sf)
  }

  # Get bounding box of boundary
  bbox <- sf::st_bbox(boundary)

  # Disable S2 for all spatial operations (spherical geometry can cause issues)
  s2_state <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_state), add = TRUE)

  # Generate rectangular grid covering the boundary
  rect_grid <- grid_rect(bbox, g)
  rect_grid <- sf::st_make_valid(rect_grid)

  # Get boundary geometry and ensure it's valid
  boundary_geom <- sf::st_geometry(boundary)
  boundary_geom <- sf::st_make_valid(boundary_geom)
  if (inherits(boundary, "sf") || length(boundary_geom) > 1) {
    boundary_geom <- sf::st_union(boundary_geom)
    boundary_geom <- sf::st_make_valid(boundary_geom)
  }

  if (crop) {
    # Crop hexagons to boundary
    result <- tryCatch({
      suppressWarnings(sf::st_intersection(rect_grid, boundary_geom))
    }, error = function(e) {
      # If intersection fails, try with buffered geometries
      boundary_buf <- sf::st_buffer(boundary_geom, 0)
      rect_buf <- sf::st_buffer(rect_grid, 0)
      suppressWarnings(sf::st_intersection(rect_buf, boundary_buf))
    })
    # Keep only polygon geometries (filter out points/lines from edge cases)
    geom_types <- sf::st_geometry_type(result)
    result <- result[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
  } else {
    # Filter to hexagons whose centroids fall within boundary
    centroids <- suppressWarnings(sf::st_centroid(rect_grid))
    within <- suppressWarnings(
      sf::st_within(centroids, boundary_geom, sparse = FALSE)
    )
    result <- rect_grid[apply(within, 1, any), ]
  }

  result
}

# =============================================================================
# CELL AREA COMPUTATION
# =============================================================================

#' Compute per-cell area in km²
#'
#' Returns the area of each cell in square kilometers. For ISEA grids, all
#' cells have the same area (equal-area property). For H3 grids, each cell
#' has a different geodesic area depending on its location.
#'
#' @param cell_id Cell IDs to compute area for. For ISEA grids, these are
#'   numeric; for H3 grids, character strings. When \code{grid} is a HexData
#'   object and \code{cell_id} is \code{NULL}, all cell IDs from the data are
#'   used.
#' @param grid A HexGridInfo or HexData object.
#'
#' @return Named numeric vector of areas in km², one per \code{cell_id}.
#'
#' @details
#' For ISEA grids the area is constant across all cells and is read directly
#' from the grid specification.
#'
#' For H3 grids the area varies by latitude. This function computes geodesic
#' area via \code{sf::st_area()} on H3 cell polygons, with results cached in a
#' session-scoped environment so repeated calls for the same cells are fast.
#'
#' @seealso \code{\link{hex_grid}} for grid specifications,
#'   \code{\link{h3_crosswalk}} for ISEA/H3 interoperability
#'
#' @export
#' @examples
#' # ISEA: constant area
#' grid <- hex_grid(area_km2 = 1000)
#' cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
#' cell_area(cells, grid)
#'
#' # H3: area varies by location
#' \donttest{
#' h3 <- hex_grid(resolution = 5, type = "h3")
#' h3_cells <- lonlat_to_cell(c(0, 0), c(0, 80), h3)
#' cell_area(h3_cells, h3)  # equator vs polar — different areas
#' }
cell_area <- function(cell_id = NULL, grid) {

  # Handle HexData input
  if (is_hex_data(grid)) {
    if (is.null(cell_id)) {
      cell_id <- grid@cell_id
    }
    g <- grid@grid
  } else {
    g <- extract_grid(grid)
    if (is.null(cell_id)) {
      stop("cell_id required when grid is not HexData")
    }
  }

  # ISEA: constant equal-area
  if (!is_h3_grid(g)) {
    areas <- rep(g@area_km2, length(cell_id))
    if (is.numeric(cell_id)) {
      names(areas) <- as.character(as.integer(cell_id))
    } else {
      names(areas) <- as.character(cell_id)
    }
    return(areas)
  }

  # H3: per-cell area via native C backend
  cell_id <- as.character(cell_id)
  areas <- cpp_h3_cellAreaKm2(cell_id)
  names(areas) <- cell_id
  areas
}


# =============================================================================
# HIERARCHICAL INDEX HELPERS
# =============================================================================

#' Convert cell ID to hierarchical index string
#'
#' Advanced function for working with hierarchical index strings.
#' Most users don't need this - use cell IDs directly.
#'
#' @param cell_id Numeric vector of cell IDs
#' @param grid A HexGridInfo or HexData object
#'
#' @return Character vector of hierarchical index strings
#'
#' @keywords internal
#' @export
cell_to_index <- function(cell_id, grid) {
  g <- extract_grid(grid)

  # H3 cell IDs are already hierarchical index strings
  if (is_h3_grid(g)) {
    return(as.character(cell_id))
  }

  # Determine index type based on aperture
  index_type <- if (g@aperture == "3") "z3"
                else if (g@aperture == "7") "z7"
                else "zorder"


  # Convert aperture to integer for C++ functions
  aperture_int <- if (g@aperture == "4/3") 3L else as.integer(g@aperture)

  sapply(cell_id, function(id) {
    # Get quad/ij coordinates
    qij <- cpp_cell_to_quad_ij(id, g@resolution, aperture_int)
    # Encode to index
    cpp_cell_to_index(
      qij$quad, qij$i, qij$j,
      g@resolution, aperture_int, index_type
    )
  })
}

#' Get parent cell
#'
#' Returns the parent cell at a coarser resolution.
#'
#' @param cell_id Numeric vector of cell IDs
#' @param grid A HexGridInfo or HexData object
#' @param levels Number of levels up (default 1)
#'
#' @return Numeric vector of parent cell IDs
#'
#' @keywords internal
#' @export
#' @examples
#' grid <- hex_grid(resolution = 10)
#' child_cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
#' parent_cells <- get_parent(child_cells, grid)
get_parent <- function(cell_id, grid, levels = 1L) {
  g <- extract_grid(grid)

  if (g@resolution < levels) {
    stop("Cannot get parent: already at minimum resolution")
  }

  # H3 path
  if (is_h3_grid(g)) {
    parent_res <- g@resolution - as.integer(levels)
    if (parent_res < H3_MIN_RESOLUTION) {
      stop("Cannot get parent: would go below H3 minimum resolution")
    }
    return(cpp_h3_cellToParent(as.character(cell_id), parent_res))
  }

  index_type <- if (g@aperture == "3") "z3"
                else if (g@aperture == "7") "z7"
                else "zorder"

  # Convert aperture to integer for C++ functions
  aperture_int <- if (g@aperture == "4/3") 3L else as.integer(g@aperture)

  # Get index, get parent, convert back
  parent_res <- g@resolution - levels

  sapply(cell_id, function(id) {
    # Get quad/ij at current resolution
    qij <- cpp_cell_to_quad_ij(id, g@resolution, aperture_int)

    # Get index string
    idx <- cpp_cell_to_index(qij$quad, qij$i, qij$j,
                             g@resolution, aperture_int, index_type)

    # Get parent index
    parent_idx <- cpp_get_parent_index(idx, aperture_int, index_type)

    # Convert back to cell ID at parent resolution
    # cpp_index_to_cell returns face, i, j, resolution - not cell_id
    result <- cpp_index_to_cell(parent_idx, aperture_int, index_type)
    # Convert quad/ij coordinates to cell ID
    cpp_quad_ij_to_cell(result$face, result$i, result$j, result$resolution, aperture_int)
  })
}

#' Get children cells
#'
#' Returns the child cells at a finer resolution.
#'
#' @param cell_id Numeric vector of cell IDs
#' @param grid A HexGridInfo or HexData object
#' @param levels Number of levels down (default 1)
#'
#' @return List of numeric vectors containing child cell IDs
#'
#' @keywords internal
#' @export
get_children <- function(cell_id, grid, levels = 1L) {
  g <- extract_grid(grid)

  # H3 path
  if (is_h3_grid(g)) {
    child_res <- g@resolution + as.integer(levels)
    if (child_res > H3_MAX_RESOLUTION) {
      stop("Cannot get children: would exceed H3 maximum resolution (15)")
    }
    return(cpp_h3_cellToChildren(as.character(cell_id), child_res))
  }

  if (g@resolution + levels > MAX_RESOLUTION) {
    stop("Cannot get children: would exceed maximum resolution")
  }

  index_type <- if (g@aperture == "3") "z3"
                else if (g@aperture == "7") "z7"
                else "zorder"

  # Convert aperture to integer for C++ functions
  aperture_int <- if (g@aperture == "4/3") 3L else as.integer(g@aperture)

  lapply(cell_id, function(id) {
    qij <- cpp_cell_to_quad_ij(id, g@resolution, aperture_int)
    idx <- cpp_cell_to_index(qij$quad, qij$i, qij$j,
                             g@resolution, aperture_int, index_type)

    # Get children indices
    children_idx <- cpp_get_children_indices(idx, aperture_int, index_type)

    # Convert to cell IDs
    # cpp_index_to_cell returns face, i, j, resolution - not cell_id
    sapply(children_idx, function(child_idx) {
      result <- cpp_index_to_cell(child_idx, aperture_int, index_type)
      # Convert quad/ij coordinates to cell ID
      cpp_quad_ij_to_cell(result$face, result$i, result$j, result$resolution, aperture_int)
    })
  })
}
