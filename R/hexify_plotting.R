# hexify_plotting.R
# Visualization and boundary generation functions
#
# This file contains functions for generating cell boundaries and grids
# for plotting and visualization purposes.

#' @title Plotting and Visualization
#' @description Functions for generating cell boundaries and grids for plotting
#' @name hexify-plotting
NULL

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Get boundary coordinates for a single cell (internal)
#'
#' @param cell_id Cell index string
#' @param aperture Grid aperture
#' @param index_type Index encoding type
#' @return List with lon and lat vectors (closed polygon)
#' @noRd
.get_cell_boundary <- function(cell_id, aperture, index_type) {
  # Decode index to get cell parameters

cell_info <- hexify_index_to_cell(cell_id, aperture, index_type)

  # Extract values
  face_val <- as.integer(cell_info$face)
  i_val <- as.integer(cell_info$i)
  j_val <- as.integer(cell_info$j)
  res_val <- as.integer(cell_info$resolution)

  # Get corners in face plane coordinates
  corners <- if (aperture == 3) {
    cpp_hex_corners_ap3(i_val, j_val, res_val, hex_radius = 1.0)
  } else if (aperture == 4) {
    cpp_hex_corners_ap4(i_val, j_val, res_val, hex_radius = 1.0)
  } else {
    cpp_hex_corners_ap7(i_val, j_val, res_val, hex_radius = 1.0)
  }

  # Convert each corner to lon/lat
  n_corners <- length(corners$x)
  corner_lons <- numeric(n_corners)
  corner_lats <- numeric(n_corners)

  for (j in seq_len(n_corners)) {
    ll <- cpp_face_xy_to_ll(corners$x[j], corners$y[j], face_val)
    corner_lons[j] <- ll["lon"]
    corner_lats[j] <- ll["lat"]
  }

  # Close the polygon
  list(
    lon = c(corner_lons, corner_lons[1]),
    lat = c(corner_lats, corner_lats[1])
  )
}

# =============================================================================
# PUBLIC API
# =============================================================================

#' Get cell boundaries as polygon coordinates
#'
#' Generates the boundary coordinates for hexagonal cells. Returns a list
#' of polygon coordinates (lon/lat) for each cell, suitable for plotting.
#'
#' @param grid Grid specification from hexify_construct()
#' @param seqnum Cell indices (character vector)
#'
#' @return List of boundary coordinates, one element per cell. Each element
#'   contains:
#'   \item{lon}{Longitude coordinates of boundary vertices}
#'   \item{lat}{Latitude coordinates of boundary vertices}
#'
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_construct(area = 1000, aperture = 3)
#' cells <- hexify_geo_to_cell(grid, lon = c(0, 5), lat = c(0, 45))
#' boundaries <- hexify_cell_to_boundary(grid, cells$seqnum)
#'
#' # Plot first cell
#' plot(boundaries[[1]]$lon, boundaries[[1]]$lat, type = "l")
#' }
hexify_cell_to_boundary <- function(grid, seqnum) {
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_construct()")
  }

  boundaries <- vector("list", length(seqnum))
  for (i in seq_along(seqnum)) {
    if (is.na(seqnum[i])) {
      boundaries[[i]] <- NA
    } else {
      boundaries[[i]] <- .get_cell_boundary(seqnum[i], grid$aperture, grid$index_type)
    }
  }
  boundaries
}

#' Convert cell indices to spatial grid for plotting
#'
#' Converts a set of cell indices to polygon geometries suitable for
#' plotting and spatial analysis. Returns an sf object or data frame
#' depending on the return_sf parameter.
#'
#' @param dggs Grid specification from hexify_construct() or dgconstruct()
#' @param cells Character vector of cell indices
#' @param savegrid Optional file path to save as shapefile
#' @param return_sf Logical. If TRUE, returns sf object; if FALSE, returns
#'   long-format data frame (faster, more memory efficient)
#'
#' @return If return_sf = TRUE: sf object with polygon geometries
#'   If return_sf = FALSE: data frame with columns seqnum, order, long, lat
#'   If savegrid is specified: file path of saved shapefile
#'
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(sf)
#'
#' grid <- hexify_construct(area = 10000, aperture = 3)
#'
#' # Generate some cells
#' lons <- seq(-10, 10, by = 5)
#' lats <- seq(-10, 10, by = 5)
#' cells <- hexify_geo_to_cell(grid, lons, lats)$seqnum
#'
#' # Get as sf object (for ggplot2, tmap, etc.)
#' hex_sf <- dgcellstogrid(grid, cells, return_sf = TRUE)
#' plot(st_geometry(hex_sf))
#'
#' # Get as data frame (faster, for base R plotting)
#' hex_df <- dgcellstogrid(grid, cells, return_sf = FALSE)
#'
#' # Save as shapefile
#' dgcellstogrid(grid, cells, savegrid = "hexagons.shp")
#' }
dgcellstogrid <- function(dggs, cells, savegrid = NA, return_sf = TRUE) {
  if (!inherits(dggs, "hexify_grid")) {
    stop("dggs must be a hexify_grid object from hexify_construct()")
  }

  if (return_sf && !requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for return_sf = TRUE. ",
         "Install with: install.packages('sf')")
  }

  cells <- cells[!is.na(cells)]
  if (length(cells) == 0) stop("No valid cells provided")

  unique_cells <- unique(cells)

  if (return_sf) {
    polygons <- vector("list", length(unique_cells))
    for (i in seq_along(unique_cells)) {
      boundary <- .get_cell_boundary(unique_cells[i], dggs$aperture, dggs$index_type)
      poly_matrix <- cbind(boundary$lon, boundary$lat)
      polygons[[i]] <- sf::st_polygon(list(poly_matrix))
    }
    sfc <- sf::st_sfc(polygons, crs = 4326)
    result <- sf::st_sf(seqnum = unique_cells, geometry = sfc)
  } else {
    df_list <- vector("list", length(unique_cells))
    for (i in seq_along(unique_cells)) {
      boundary <- .get_cell_boundary(unique_cells[i], dggs$aperture, dggs$index_type)
      df_list[[i]] <- data.frame(
        seqnum = unique_cells[i],
        order = seq_along(boundary$lon),
        long = boundary$lon,
        lat = boundary$lat,
        stringsAsFactors = FALSE
      )
    }
    result <- do.call(rbind, df_list)
    rownames(result) <- NULL
  }

  if (!is.na(savegrid)) {
    if (!return_sf) {
      warning("Converting to sf object for shapefile export")
      result <- dgcellstogrid(dggs, cells, savegrid = NA, return_sf = TRUE)
    }
    sf::st_write(result, savegrid, delete_dsn = TRUE, quiet = TRUE)
    return(savegrid)
  }

  result
}

#' Plot hexagonal cells
#'
#' Simple plotting function for hexagonal cells. For more sophisticated
#' plotting, convert to sf object and use ggplot2, tmap, or other mapping
#' packages.
#'
#' @param grid Grid specification
#' @param cells Cell indices to plot
#' @param col Fill color for cells
#' @param border Border color for cells
#' @param add If TRUE, add to existing plot; if FALSE, create new plot
#' @param ... Additional arguments passed to polygon()
#'
#' @return NULL (creates plot as side effect)
#'
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_construct(area = 10000, aperture = 3)
#' cells <- hexify_geo_to_cell(grid, lon = seq(-10, 10, 5),
#'                                     lat = seq(-10, 10, 5))$seqnum
#'
#' # Simple plot
#' hexify_plot_cells(grid, cells, col = "lightblue", border = "blue")
#'
#' # Add to existing world map
#' library(maps)
#' map("world")
#' hexify_plot_cells(grid, cells, col = "red", add = TRUE)
#' }
hexify_plot_cells <- function(grid, cells, col = "lightgray",
                              border = "black", add = FALSE, ...) {
  boundaries <- hexify_cell_to_boundary(grid, cells)
  boundaries <- boundaries[!is.na(boundaries)]

  if (length(boundaries) == 0) {
    warning("No valid cells to plot")
    return(invisible(NULL))
  }

  if (!add) {
    all_lons <- unlist(lapply(boundaries, function(b) b$lon))
    all_lats <- unlist(lapply(boundaries, function(b) b$lat))
    plot(range(all_lons), range(all_lats), type = "n",
         xlab = "Longitude", ylab = "Latitude",
         main = sprintf("Hexify Grid (res=%d, aperture=%d)",
                       grid$resolution, grid$aperture))
  }

  for (boundary in boundaries) {
    polygon(boundary$lon, boundary$lat, col = col, border = border, ...)
  }

  invisible(NULL)
}
