# grid_helpers.R
# Helper functions that accept HexGridInfo or HexData objects
#
# These functions wrap the low-level coordinate conversion functions
# to accept grid specifications, eliminating the need to repeat
# aperture/resolution parameters.

#' Make a ring's longitudes continuous
#'
#' Adds or drops whole turns so that consecutive corners stay less than half a
#' turn apart, which is how far apart they are on the sphere, and centres the
#' result on the map. A ring crossing the antimeridian then runs past +/-180
#' in one piece rather than jumping the width of the map, and
#' `sf::st_wrap_dateline()` splits it at the seam.
#'
#' @param coords A 2+ column matrix whose first column is longitude, the first
#'   corner repeated last
#' @return The same matrix, with continuous longitudes
#' @noRd
unwrap_ring_longitudes <- function(coords) {
  lons <- coords[, 1]
  n <- length(lons)
  if (n < 2 || anyNA(lons)) {
    return(coords)
  }

  steps <- diff(lons)

  # A step over a pole is half a turn either way. Take the way that closes the
  # ring; the cell lies beside the pole rather than around it.
  over_pole <- abs(abs(steps) - 180) < 1e-6
  wrapped <- !over_pole & abs(steps) > 180
  steps[wrapped] <- steps[wrapped] - 360 * sign(steps[wrapped])
  if (sum(over_pole) == 1L) {
    steps[over_pole] <- -sum(steps[!over_pole])
  }

  lons <- c(lons[1], lons[1] + cumsum(steps))

  # The last corner is the first one, whole turns on: place it exactly, so that
  # accumulated rounding cannot leave the ring open.
  lons[n] <- lons[1] + 360 * round((lons[n] - lons[1]) / 360)

  center <- (max(lons) + min(lons)) / 2
  coords[, 1] <- lons - 360 * round(center / 360)
  coords
}

#' Latitude where a cell edge meets a meridian
#'
#' The edge is the great-circle arc between two corners. Defined for corners
#' that are neither on the same meridian nor half a turn apart, which is every
#' edge that meets a meridian at one point.
#'
#' @param lon1,lat1,lon2,lat2 The two corners, in degrees
#' @param lon The meridian, in degrees on the same continuous scale
#' @return The latitude in degrees
#' @noRd
edge_lat_at_lon <- function(lon1, lat1, lon2, lat2, lon) {
  rad <- pi / 180
  atan((tan(lat1 * rad) * sin((lon2 - lon) * rad) +
        tan(lat2 * rad) * sin((lon - lon1) * rad)) /
       sin((lon2 - lon1) * rad)) / rad
}

#' Close a ring that encircles a pole
#'
#' The corners of a cell holding a pole run once around the globe, so no lon/lat
#' ring closes on them. Cutting the ring at the antimeridian and carrying the two
#' ends up to the pole closes it around the cap the cell covers, within the map
#' and without the seam-crossing a whole-globe ring would otherwise need.
#'
#' @param coords A 2+ column matrix with continuous longitudes, the first corner
#'   repeated last
#' @return The same matrix when the ring closes on its own, otherwise the ring
#'   from -180 to 180 with its two polar corners
#' @noRd
close_ring_over_pole <- function(coords) {
  n <- nrow(coords)
  if (n < 4 || anyNA(coords[, 1])) {
    return(coords)
  }

  lons <- coords[, 1]
  winding <- lons[n] - lons[1]
  pole <- sign(mean(coords[, 2]))
  if (abs(abs(winding) - 360) > 1e-6 || pole == 0) {
    return(coords)
  }

  # The one antimeridian the ring spans, and the edge that meets it
  seam <- 180 + 360 * round((mean(lons) - 180) / 360)
  at <- which((lons[-n] - seam) * (lons[-1] - seam) < 0)[1]
  if (is.na(at)) {
    return(coords)
  }

  opening <- coords[at, , drop = FALSE]
  opening[, 1] <- seam
  opening[, 2] <- edge_lat_at_lon(lons[at], coords[at, 2],
                                  lons[at + 1L], coords[at + 1L, 2], seam)
  closing <- opening
  closing[, 1] <- seam + winding

  # Walk from the seam once around, so the ring runs the width of the map
  turned <- coords[seq_len(at), , drop = FALSE]
  turned[, 1] <- turned[, 1] + winding
  ring <- rbind(opening, coords[seq_len(n - 1L)[-seq_len(at)], , drop = FALSE],
                turned, closing)
  ring[, 1] <- ring[, 1] - seam - 180 * sign(winding)

  caps <- ring[c(nrow(ring), 1L), , drop = FALSE]
  caps[, 2] <- 90 * pole
  rbind(ring, caps, ring[1, , drop = FALSE])
}

#' Give a ring a corner where it meets the antimeridian
#'
#' `sf::st_wrap_dateline()` cuts a ring at +/-180 and joins the cut ends along a
#' line of constant latitude, which is not where the cell's own edge runs. A
#' corner placed on the edge itself is where the cut belongs, so the two halves
#' keep the cell's shape. An edge running over a pole meets the meridian at the
#' pole.
#'
#' @param coords A 2+ column matrix with continuous longitudes
#' @return The same matrix, with a corner inserted at each seam crossing
#' @noRd
insert_seam_corners <- function(coords) {
  lons <- coords[, 1]
  n <- nrow(coords)
  if (n < 2 || anyNA(lons) || (min(lons) >= -180 && max(lons) <= 180)) {
    return(coords)
  }

  pieces <- vector("list", 2L * n)
  at <- 0L
  for (i in seq_len(n - 1L)) {
    pieces[[at <- at + 1L]] <- coords[i, , drop = FALSE]
    over_pole <- abs(abs(lons[i + 1L] - lons[i]) - 180) < 1e-6

    for (seam in c(-180, 180)) {
      if ((lons[i] - seam) * (lons[i + 1L] - seam) >= 0) {
        next
      }
      corner <- coords[i, , drop = FALSE]
      corner[, 1] <- seam
      corner[, 2] <- if (over_pole) {
        90 * sign(coords[i, 2] + coords[i + 1L, 2])
      } else {
        edge_lat_at_lon(lons[i], coords[i, 2], lons[i + 1L], coords[i + 1L, 2], seam)
      }
      pieces[[at <- at + 1L]] <- corner
    }
  }
  pieces[[at + 1L]] <- coords[n, , drop = FALSE]

  do.call(rbind, pieces[seq_len(at + 1L)])
}

#' Prepare a spherical cell ring for a lon/lat polygon
#'
#' @param coords A 2+ column matrix whose first two columns are lon and lat
#' @return The prepared coordinate matrix
#' @noRd
lonlat_ring_coords <- function(coords) {
  insert_seam_corners(close_ring_over_pole(unwrap_ring_longitudes(coords)))
}

#' Split cells at the antimeridian
#'
#' Only cells reaching past +/-180 are handed to `sf::st_wrap_dateline()`, which
#' cuts any segment half a turn wide and would otherwise also split a cell whose
#' edge runs over a pole rather than over the seam.
#'
#' @param x An sf object of cell polygons with continuous longitudes
#' @return The same object, with seam-crossing cells split in two
#' @noRd
wrap_cells_at_dateline <- function(x) {
  crossing <- vapply(sf::st_geometry(x), function(g) {
    lons <- sf::st_coordinates(g)[, 1]
    length(lons) > 0 && (min(lons) < -180 || max(lons) > 180)
  }, logical(1))

  if (any(crossing)) {
    wrapped <- sf::st_wrap_dateline(x[crossing, ],
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)
    sf::st_geometry(x)[crossing] <- sf::st_geometry(wrapped)
  }
  x
}

#' Build hexagon polygons for ISEA cell IDs
#'
#' Rings are prepared with `lonlat_ring_coords()` so each polygon is
#' contiguous; callers that render on a flat map pass the result through
#' `sf::st_wrap_dateline()` to split them at +/-180.
#'
#' @param cell_id Numeric vector of cell IDs
#' @param resolution Grid resolution level
#' @param aperture Grid aperture: 3, 4, 7, or a mixed sequence spelling
#' @param crs CRS the polygons carry, as sf reads it
#' @return An sfc of POLYGON geometries, one per cell ID, in input order
#' @noRd
isea_cells_to_sfc <- function(cell_id, resolution, aperture, crs = 4326) {
  corners_list <- if (is_mixed_aperture(aperture)) {
    mixed_cell_corners(cell_id, resolution, aperture)
  } else {
    cpp_cell_to_corners(
      as.numeric(cell_id),
      as.integer(resolution),
      as.integer(aperture)
    )
  }

  polygons <- lapply(corners_list, function(coords) {
    sf::st_polygon(list(lonlat_ring_coords(coords)))
  })

  # A ring that crosses the antimeridian or runs over a pole carries longitudes
  # outside -180..180; st_wrap_dateline() brings those back inside the map.
  sf::st_sfc(polygons, crs = crs)
}

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

  if (is_mixed_aperture(g@aperture)) {
    cpp_lonlat_to_cell_seq(
      as.numeric(lon),
      as.numeric(lat),
      grid_ap_seq(g)
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

  if (is_mixed_aperture(g@aperture)) {
    cpp_cell_to_lonlat_seq(
      as.numeric(cell_id),
      grid_ap_seq(g)
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
#'   \code{\link[=st_as_sf.HexData]{st_as_sf}} for converting HexData to sf
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
      coords <- lonlat_ring_coords(coords)
      sf::st_polygon(list(coords))
    })
    sfc <- sf::st_sfc(polygons, crs = grid_crs(g))
    result_sf <- sf::st_sf(cell_id = as.character(cell_id), geometry = sfc)
    if (wrap_dateline) {
      result_sf <- wrap_cells_at_dateline(result_sf)
    }
    return(result_sf)
  }

  # ISEA path: generate polygons using C++ function. For globe/orthographic
  # projections, pass wrap_dateline = FALSE to keep cells intact.
  sfc <- isea_cells_to_sfc(cell_id, g@resolution, g@aperture,
                           crs = grid_crs(g))

  result_sf <- sf::st_sf(cell_id = cell_id, geometry = sfc)
  if (wrap_dateline) {
    result_sf <- wrap_cells_at_dateline(result_sf)
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
  spacing_deg <- diagonal / km_per_degree(grid_radius_km(g)) * 0.8

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
  n_cells <- aperture_n_cells(g@aperture, g@resolution)
  if (n_cells > 100000) {
    warning(sprintf(
      "This will generate approximately %.0f cells. Consider larger area_km2.",
      n_cells
    ))
  }

  # Dense sampling - use diagonal_km from grid if available
  diagonal <- if (!is.na(g@diagonal_km)) g@diagonal_km else sqrt(g@area_km2 * 2 / sqrt(3))
  spacing_deg <- diagonal / km_per_degree(grid_radius_km(g)) * 0.7

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
#' For H3 grids the area varies from cell to cell, by about a factor of 2
#' across the hexagons of a resolution. The vendored 'H3' library computes
#' each cell's spherical polygon area as a solid angle, which this function
#' reads on the grid's body, so a grid built with \code{radius_km} reports that
#' body's areas.
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
    names(areas) <- as.character(cell_id)
    return(areas)
  }

  # H3: per-cell area via native C backend, read on the grid's body
  cell_id <- as.character(cell_id)
  areas <- scale_area_to_body(cpp_h3_cellAreaKm2(cell_id), grid_radius_km(g))
  names(areas) <- cell_id
  areas
}


# =============================================================================
# HIERARCHICAL INDEX HELPERS
# =============================================================================

#' Hierarchical index strings of cells on a pure-aperture grid
#'
#' The index encoders work on (quad, i, j), so cell IDs make the round trip
#' through `cpp_cell_to_quad_ij()` first. `cell_to_index()`, `get_parent()` and
#' `get_children()` all enter the hierarchy this way.
#'
#' @param cell_id Numeric vector of cell IDs
#' @param resolution Resolution the cell IDs belong to
#' @param aperture_int Integer aperture (3, 4 or 7)
#' @param index_type One of "z3", "z7", "zorder"
#' @return Character vector of index strings
#' @noRd
isea_cells_to_index <- function(cell_id, resolution, aperture_int, index_type) {
  qij <- cpp_cell_to_quad_ij(as.numeric(cell_id), resolution, aperture_int)
  cpp_cell_to_index(qij$quad, qij$i, qij$j, resolution, aperture_int, index_type)
}

#' Cell IDs of hierarchical index strings on a pure-aperture grid
#'
#' Inverse of `isea_cells_to_index()`. Each index carries its own resolution,
#' which is the one its cell ID belongs to, so the strings are packed back into
#' cell IDs one resolution at a time.
#'
#' @param index Character vector of index strings
#' @param aperture_int Integer aperture (3, 4 or 7)
#' @param index_type One of "z3", "z7", "zorder"
#' @return Numeric vector of cell IDs
#' @noRd
isea_index_to_cells <- function(index, aperture_int, index_type) {
  cell <- cpp_index_to_cell(as.character(index), aperture_int, index_type)
  out <- rep(NA_real_, nrow(cell))

  for (resolution in unique(stats::na.omit(cell$resolution))) {
    at_res <- !is.na(cell$resolution) & cell$resolution == resolution
    out[at_res] <- cpp_quad_ij_to_cell(cell$face[at_res], cell$i[at_res],
                                       cell$j[at_res], resolution, aperture_int)
  }

  out
}

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

  # Mixed sequences use a geometric hierarchical index (see
  # R/aperture_mixed_hierarchy.R); pure apertures use the Z7/Z3/zorder encoders.
  if (is_mixed_aperture(g@aperture)) {
    return(vapply(as.numeric(cell_id),
                  function(id) mixed_cell_to_index_one(id, g@resolution, g@aperture),
                  character(1)))
  }

  # Determine index type based on aperture
  index_type <- index_type_for_aperture(g@aperture)
  aperture_int <- aperture_to_int(g@aperture)

  isea_cells_to_index(cell_id, g@resolution, aperture_int, index_type)
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

  # Mixed sequences: geometric parent (centre re-quantised at the coarser resolution).
  if (is_mixed_aperture(g@aperture)) {
    return(mixed_get_parent(as.numeric(cell_id), g@resolution, g@aperture,
                            as.integer(levels)))
  }

  index_type <- index_type_for_aperture(g@aperture)
  aperture_int <- aperture_to_int(g@aperture)
  levels <- as.integer(levels)

  # cpp_get_parent_index() strips one level off the index string, so `levels`
  # levels up is that many strips.
  idx <- isea_cells_to_index(cell_id, g@resolution, aperture_int, index_type)
  for (step in seq_len(levels)) {
    idx <- cpp_get_parent_index(idx, aperture_int, index_type)
  }
  isea_index_to_cells(idx, aperture_int, index_type)
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

  # Mixed sequences: geometric children (cells whose geometric parent is this cell).
  if (is_mixed_aperture(g@aperture)) {
    child_res <- g@resolution + as.integer(levels)
    ncc <- aperture_n_cells(g@aperture, child_res)
    return(lapply(as.numeric(cell_id), function(id)
      mixed_get_children_one(id, g@resolution, child_res, g@aperture, ncc)))
  }

  levels <- as.integer(levels)

  # One level at a time, so each level's children are read at their own
  # resolution.
  front <- lapply(as.numeric(cell_id), identity)
  for (step in seq_len(levels)) {
    front <- isea_children_one_level(front, g@resolution + step - 1L, g)
  }
  front
}

#' Children of ISEA cells, one resolution down
#'
#' A child is a cell whose parent is this one, so the children are read back
#' from `get_parent()` over a candidate set. Two things generate candidates.
#' Appending a digit to the parent's index names the children directly, which is
#' what makes the two operations inverse, but the twelve icosahedron vertices sit
#' at the corner of several quads and so have an index spelling in each: the
#' children spelled under the other quads are not reached from the one spelling
#' `cell_to_index()` returns, and a digit that names no cell at all comes back
#' out of range or NA. The ring of child cells around the parent's centre reaches
#' those, since neighbours are found on the sphere rather than in one quad's
#' lattice. The parent test then decides which candidates are children.
#'
#' @param front List of numeric vectors of cell IDs, one per requested cell
#' @param resolution Resolution the cells in `front` are at
#' @param g The grid the cells came from
#' @return A list of the same length, each element the children, sorted
#' @noRd
isea_children_one_level <- function(front, resolution, g) {
  parents <- unique(unlist(front, use.names = FALSE))
  if (length(parents) == 0) {
    return(front)
  }

  index_type <- index_type_for_aperture(g@aperture)
  aperture_int <- aperture_to_int(g@aperture)
  child_res <- resolution + 1L
  radius <- grid_radius_km(g)

  parent_grid <- hex_grid(resolution = resolution, aperture = g@aperture,
                          radius_km = radius)
  child_grid <- hex_grid(resolution = child_res, aperture = g@aperture,
                         radius_km = radius)
  n_child <- aperture_n_cells(g@aperture, child_res)

  parent_idx <- isea_cells_to_index(parents, resolution, aperture_int, index_type)
  kids <- cpp_get_children_indices(parent_idx, aperture = aperture_int,
                                   index_type = index_type)
  flat <- isea_index_to_cells(unlist(kids, use.names = FALSE), aperture_int,
                              index_type)
  expansion <- unname(split(flat, factor(rep(seq_along(kids), lengths(kids)),
                                         levels = seq_along(kids))))

  centre <- cell_to_lonlat(parents, parent_grid)
  central <- lonlat_to_cell(centre$lon_deg, centre$lat_deg, child_grid)
  ring <- get_neighbors(central, child_grid)

  candidates <- lapply(seq_along(parents), function(k) {
    cand <- unique(c(expansion[[k]], central[k], ring[[k]]))
    cand[!is.na(cand) & cand >= 1 & cand <= n_child]
  })

  pool <- unique(unlist(candidates, use.names = FALSE))
  parent_of <- stats::setNames(get_parent(pool, child_grid), as.character(pool))

  children <- lapply(seq_along(parents), function(k) {
    cand <- candidates[[k]]
    sort(cand[parent_of[as.character(cand)] == parents[k]])
  })
  names(children) <- as.character(parents)

  lapply(front, function(ids) {
    sort(unique(unlist(children[as.character(ids)], use.names = FALSE)))
  })
}
