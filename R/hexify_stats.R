# hexify_stats.R
# Grid statistics and utility functions
#
# This file contains functions for calculating grid statistics and
# other utility functions for working with hexagonal grids.

#' @title Grid Statistics
#' @description Functions for calculating grid statistics and utilities
#' @name hexify-stats
NULL

#' Get grid statistics for whole-body coverage
#'
#' Calculates statistics about the hexagonal grid at the current resolution,
#' including total number of cells, cell area, and cell spacing. Figures follow
#' the radius the grid carries, so a grid built with \code{radius_km} reports
#' that body's areas and distances.
#'
#' @param dggs Grid specification from hexify_grid()
#'
#' @return List with components:
#'   \item{area_km}{Total surface area of the body in km^2}
#'   \item{n_cells}{Total number of cells at this resolution}
#'   \item{cell_area_km2}{Average cell area in km^2}
#'   \item{cell_spacing_km}{Average distance between cell centers in km}
#'   \item{resolution}{Resolution level}
#'   \item{aperture}{Grid aperture}
#'
#' @family grid statistics
#' @export
#' @examples
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' stats <- dgearthstat(grid)
#'
#' print(sprintf("Resolution %d has %.0f cells",
#'               stats$resolution, stats$n_cells))
#' print(sprintf("Average cell area: %.2f km^2",
#'               stats$cell_area_km2))
#' print(sprintf("Average cell spacing: %.2f km",
#'               stats$cell_spacing_km))
dgearthstat <- function(dggs) {
  if (is_hex_grid(dggs)) {
    if (is_h3_grid(dggs)) {
      return(h3_level_stats(dggs@resolution, grid_radius_km(dggs)))
    }
    return(grid_level_stats(body_surface_km2(grid_radius_km(dggs)),
                            grid_n_cells(dggs),
                            dggs@resolution, aperture_to_int(dggs@aperture)))
  }

  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a hexify_grid or HexGridInfo object")
  }

  resolution <- get_grid_resolution(dggs, require = TRUE)
  grid_level_stats(body_surface_km2(grid_radius_km(dggs)),
                   aperture_n_cells(dggs$aperture, resolution, grid_polyhedron(dggs)),
                   resolution, dggs$aperture)
}

#' Whole-body statistics of one grid level
#'
#' @param surface_km2 Area of the body in km^2
#' @param n_cells Number of cells at this level
#' @param resolution Resolution, reported as given
#' @param aperture Aperture, reported as given
#' @param cell_area_km2 Mean cell area in km^2: the body's area over the cell
#'   count, unless the backend reports its own
#' @param ... Further fields appended to the list
#' @return The list dgearthstat() documents
#' @noRd
grid_level_stats <- function(surface_km2, n_cells, resolution, aperture,
                             cell_area_km2 = surface_km2 / n_cells, ...) {
  c(list(
    area_km = surface_km2,
    n_cells = n_cells,
    cell_area_km2 = cell_area_km2,
    cell_spacing_km = hex_spacing_km(cell_area_km2),
    cls_km = cls_km(cell_area_km2),
    resolution = resolution,
    aperture = aperture
  ), list(...))
}

#' Whole-body statistics of an H3 level
#'
#' H3 cells are not equal-area, so the mean cell area is H3's own average
#' rather than the body's area over the cell count.
#' @param resolution H3 resolution
#' @param radius_km Radius of the body in km
#' @noRd
h3_level_stats <- function(resolution, radius_km) {
  grid_level_stats(body_surface_km2(radius_km), h3_n_cells(resolution),
                   resolution, 7L,
                   cell_area_km2 = h3_avg_area_km2(resolution, radius_km),
                   grid_type = "h3")
}

#' Find closest resolution for target cell area
#'
#' Finds the grid resolution that produces cells closest to the target area.
#' This is primarily used internally by \code{\link{hexify_grid}} and
#' \code{\link{hex_grid}}. Most users should use those functions directly.
#'
#' @param dggs Grid specification (aperture and topology must be set)
#' @param area Target cell area in km^2 (if metric=TRUE)
#' @param round Rounding method ("nearest", "up", "down")
#' @param metric Whether area is in metric units
#' @param show_info Print information about chosen resolution
#'
#' @return Resolution level (integer)
#'
#' @keywords internal
#' @export
#' @examples
#' # Create a temporary grid to get aperture settings
#' temp_grid <- list(aperture = 3, topology = "HEXAGON")
#' class(temp_grid) <- "hexify_grid"
#' 
#' # Find resolution for 1000 km^2 cells
#' res <- dg_closest_res_to_area(temp_grid, area = 1000, 
#'                                metric = TRUE, show_info = TRUE)
#' print(res)
#' @family grid statistics
dg_closest_res_to_area <- function(dggs, area, round = "nearest",
                                   metric = TRUE, show_info = FALSE) {
  if (!metric) {
    # Convert from square miles to square km
    area <- area * MI2_TO_KM2
  }

  resolution <- resolve_resolution_from_area(area, dggs$aperture,
                                             grid_radius_km(dggs), round)

  if (show_info) {
    # Calculate actual area at this resolution
    temp_grid <- dggs
    temp_grid$resolution <- resolution
    temp_grid$res <- resolution
    stats <- dgearthstat(temp_grid)
    
    message(sprintf("Resolution %d:", resolution))
    message(sprintf("  Cell area: %.2f km^2", stats$cell_area_km2))
    message(sprintf("  Cell spacing: %.2f km", stats$cell_spacing_km))
    message(sprintf("  Total cells: %.0f", stats$n_cells))
  }
  
  return(resolution)
}


#' Compare grid resolutions
#'
#' Generates a table comparing different resolution levels for a given
#' grid configuration. Useful for choosing appropriate resolution.
#'
#' @param aperture Grid aperture (3, 4, or 7). Ignored for H3 grids.
#' @param res_range Range of resolutions to compare (e.g., 1:10)
#' @param type Grid type: "isea" (default) or "h3".
#' @param print If TRUE, prints a formatted table to console. If FALSE (default),
#'   returns a data frame.
#' @param radius_km Radius of the body, in kilometers, or a body name such as
#'   "mars" (default Earth). See \code{\link{hex_grid}}.
#'
#' @return If print=FALSE: data frame with columns resolution, n_cells,
#'   cell_area_km2, cell_spacing_km, cls_km.
#'   If print=TRUE: invisibly returns the data frame after printing.
#'
#' @family grid statistics
#' @export
#' @examples
#' # Get data frame of resolutions 0-10 for aperture 3
#' comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:10)
#' print(comparison)
#'
#' # Print formatted table directly
#' hexify_compare_resolutions(aperture = 3, res_range = 0:10, print = TRUE)
#'
#' # Find resolution with cells ~1000 km^2
#' subset(comparison, cell_area_km2 > 900 & cell_area_km2 < 1100)
#'
#' # Resolutions on another body
#' hexify_compare_resolutions(aperture = 3, res_range = 0:6, radius_km = "mars")
hexify_compare_resolutions <- function(aperture = 3, res_range = 0:15,
                                       type = c("isea", "h3"),
                                       print = FALSE,
                                       radius_km = EARTH_RADIUS_KM) {
  type <- match.arg(type)
  radius_km <- resolve_radius_km(radius_km)

  if (type == "h3") {
    res_range <- res_range[res_range >= H3_MIN_RESOLUTION & res_range <= H3_MAX_RESOLUTION]
    level_stats <- function(res) h3_level_stats(res, radius_km)
    title <- "Grid Resolution Comparison (H3)"
  } else {
    surface_km2 <- body_surface_km2(radius_km)
    res_range <- res_range[res_range >= MIN_RESOLUTION &
                             res_range <= isea_max_resolution(aperture)]
    level_stats <- function(res) {
      grid_level_stats(surface_km2, aperture_n_cells(aperture, res), res, aperture)
    }
    title <- sprintf("Grid Resolution Comparison (Aperture %s)", aperture)
  }

  columns <- c("resolution", "n_cells", "cell_area_km2", "cell_spacing_km", "cls_km")
  result_df <- do.call(rbind, lapply(res_range, function(res) {
    as.data.frame(level_stats(res)[columns])
  }))

  if (print) {
    if (type == "h3") {
      .print_resolution_table(result_df, title, "%-12.4f  %-12.3f  %-10.3f",
                              "Note: H3 areas are averages; actual area varies by location")
    } else {
      .print_resolution_table(result_df, title, "%-12.1f  %-12.1f  %-10.1f")
    }
    return(invisible(result_df))
  }

  result_df
}

#' Print formatted resolution table
#' @param comparison Data frame from hexify_compare_resolutions()
#' @param title Heading line
#' @param value_format sprintf format of the area, spacing and CLS columns
#' @param note Optional line printed under the table
#' @noRd
.print_resolution_table <- function(comparison, title, value_format, note = NULL) {
  cat(sprintf("\n%s\n", title))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("%-4s  %-12s  %-12s  %-12s  %-10s\n",
              "Res", "# Cells", "Area (km^2)", "Spacing (km)", "CLS (km)"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in seq_len(nrow(comparison))) {
    row <- comparison[i, ]
    cat(sprintf(paste0("%-4d  %-12s  ", value_format, "\n"),
                row$resolution, format_cell_count(row$n_cells),
                row$cell_area_km2, row$cell_spacing_km, row$cls_km))
  }

  cat(paste(rep("=", 70), collapse = ""), "\n")
  if (!is.null(note)) cat(note, "\n", sep = "")
  cat("\n")
}

#' Cell count with a T/B/M/K suffix
#' @noRd
format_cell_count <- function(n) {
  if (n > 1e12) {
    sprintf("%.1fT", n / 1e12)
  } else if (n > 1e9) {
    sprintf("%.1fB", n / 1e9)
  } else if (n > 1e6) {
    sprintf("%.1fM", n / 1e6)
  } else if (n > 1e3) {
    sprintf("%.1fK", n / 1e3)
  } else {
    sprintf("%.0f", n)
  }
}
