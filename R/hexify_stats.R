# hexify_stats.R
# Grid statistics and utility functions
#
# This file contains functions for calculating grid statistics and
# other utility functions for working with hexagonal grids.

#' @title Grid Statistics
#' @description Functions for calculating grid statistics and utilities
#' @name hexify-stats
NULL

#' Get grid statistics for Earth coverage
#'
#' Calculates statistics about the hexagonal grid at the current resolution,
#' including total number of cells, cell area, and cell spacing.
#'
#' @param dggs Grid specification from hexify_grid()
#'
#' @return List with components:
#'   \item{area_km}{Total Earth surface area in km^2}
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
  # Accept HexGridInfo objects
  if (is_hex_grid(dggs)) {
    g <- dggs
    gt <- tryCatch(g@grid_type, error = function(e) "isea")

    if (gt == "h3") {
      h3_n_cells <- 2 + 120 * 7^g@resolution
      cell_area_km2 <- H3_AVG_AREA_KM2[g@resolution + 1L]
      cell_spacing_km <- sqrt(2 * cell_area_km2 / sqrt(3))
      cls_km <- 2 * sqrt(cell_area_km2 / pi)
      return(list(
        area_km = EARTH_SURFACE_KM2,
        n_cells = h3_n_cells,
        cell_area_km2 = cell_area_km2,
        cell_spacing_km = cell_spacing_km,
        cls_km = cls_km,
        resolution = g@resolution,
        aperture = 7L,
        grid_type = "h3"
      ))
    }

    # ISEA HexGridInfo
    ap <- if (g@aperture == "4/3") {
      level <- as.integer(g@resolution / 2)
      n_cells <- 10 * (4^level) * (3^(g@resolution - level)) + 2
      3L
    } else {
      ap_int <- as.integer(g@aperture)
      n_cells <- 10 * (ap_int^g@resolution) + 2
      ap_int
    }
    cell_area_km2 <- EARTH_SURFACE_KM2 / n_cells
    cell_spacing_km <- sqrt(2 * cell_area_km2 / sqrt(3))
    cls_km <- 2 * sqrt(cell_area_km2 / pi)
    return(list(
      area_km = EARTH_SURFACE_KM2,
      n_cells = n_cells,
      cell_area_km2 = cell_area_km2,
      cell_spacing_km = cell_spacing_km,
      cls_km = cls_km,
      resolution = g@resolution,
      aperture = ap
    ))
  }

  # Legacy hexify_grid / dggs path
  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a hexify_grid or HexGridInfo object")
  }

  resolution <- get_grid_resolution(dggs)


  # DGGRID cell count formula: N = 10 * aperture^res + 2
  # This accounts for the 12 pentagon cells at icosahedron vertices
  n_cells <- 10 * (dggs$aperture ^ resolution) + 2

  # Calculate cell area
  cell_area_km2 <- EARTH_SURFACE_KM2 / n_cells

  # Approximate cell spacing (distance between cell centers)
  # For hexagons: spacing approx  sqrt(2 * area / sqrt(3))
  cell_spacing_km <- sqrt(2 * cell_area_km2 / sqrt(3))

  # Calculate characteristic length scale (CLS)
  # CLS is the diameter of a spherical cap with same area
  # CLS = 2 * sqrt(area / pi)
  cls_km <- 2 * sqrt(cell_area_km2 / pi)

  return(list(
    area_km = EARTH_SURFACE_KM2,
    n_cells = n_cells,
    cell_area_km2 = cell_area_km2,
    cell_spacing_km = cell_spacing_km,
    cls_km = cls_km,
    resolution = resolution,
    aperture = dggs$aperture
  ))
}

# NOTE: dgmaxcell() is defined in dggrid_compat.R for dggridR compatibility

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
  
  # Calculate resolution
  resolution <- calculate_resolution_for_area(area, dggs$aperture)
  
  # Apply rounding
  if (round == "up") {
    resolution <- ceiling(resolution)
  } else if (round == "down") {
    resolution <- floor(resolution)
  } else {
    resolution <- round(resolution)
  }
  
  # Ensure valid range
  resolution <- max(0, min(30, resolution))
  
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
hexify_compare_resolutions <- function(aperture = 3, res_range = 0:15,
                                       type = c("isea", "h3"),
                                       print = FALSE) {
  type <- match.arg(type)

  if (type == "h3") {
    # H3 resolution table from pre-computed area values
    res_range <- res_range[res_range >= H3_MIN_RESOLUTION & res_range <= H3_MAX_RESOLUTION]
    results <- lapply(res_range, function(res) {
      cell_area_km2 <- H3_AVG_AREA_KM2[res + 1L]
      h3_n_cells <- 2 + 120 * 7^res
      cell_spacing_km <- sqrt(2 * cell_area_km2 / sqrt(3))
      cls_km <- 2 * sqrt(cell_area_km2 / pi)
      data.frame(
        resolution = res,
        n_cells = h3_n_cells,
        cell_area_km2 = cell_area_km2,
        cell_spacing_km = cell_spacing_km,
        cls_km = cls_km
      )
    })
    result_df <- do.call(rbind, results)

    if (print) {
      cat("\nGrid Resolution Comparison (H3)\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat(sprintf("%-4s  %-12s  %-12s  %-12s  %-10s\n",
                  "Res", "# Cells", "Area (km^2)", "Spacing (km)", "CLS (km)"))
      cat(paste(rep("-", 70), collapse = ""), "\n")
      for (i in seq_len(nrow(result_df))) {
        row <- result_df[i, ]
        n_cells_str <- if (row$n_cells > 1e12) {
          sprintf("%.1fT", row$n_cells / 1e12)
        } else if (row$n_cells > 1e9) {
          sprintf("%.1fB", row$n_cells / 1e9)
        } else if (row$n_cells > 1e6) {
          sprintf("%.1fM", row$n_cells / 1e6)
        } else if (row$n_cells > 1e3) {
          sprintf("%.1fK", row$n_cells / 1e3)
        } else {
          sprintf("%.0f", row$n_cells)
        }
        cat(sprintf("%-4d  %-12s  %-12.4f  %-12.3f  %-10.3f\n",
                    row$resolution, n_cells_str,
                    row$cell_area_km2, row$cell_spacing_km, row$cls_km))
      }
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat("Note: H3 areas are averages; actual area varies by location\n\n")
      return(invisible(result_df))
    }
    return(result_df)
  }

  # ISEA path (original)
  # Create temporary grid
  temp_grid <- list(
    aperture = aperture,
    topology = "HEXAGON",
    projection = "ISEA"
  )
  class(temp_grid) <- c("hexify_grid", "dggs", "list")

  # Calculate stats for each resolution
  results <- lapply(res_range, function(res) {
    temp_grid$resolution <- res
    temp_grid$res <- res
    stats <- dgearthstat(temp_grid)

    data.frame(
      resolution = res,
      n_cells = stats$n_cells,
      cell_area_km2 = stats$cell_area_km2,
      cell_spacing_km = stats$cell_spacing_km,
      cls_km = stats$cls_km
    )
  })

  # Combine into data frame
  result_df <- do.call(rbind, results)

  if (print) {
    .print_resolution_table(result_df, aperture)
    return(invisible(result_df))
  }

  return(result_df)
}

#' Print formatted resolution table
#' @noRd
.print_resolution_table <- function(comparison, aperture) {
  cat(sprintf("\nGrid Resolution Comparison (Aperture %d)\n", aperture))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("%-4s  %-12s  %-12s  %-12s  %-10s\n",
              "Res", "# Cells", "Area (km^2)", "Spacing (km)", "CLS (km)"))
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (i in seq_len(nrow(comparison))) {
    row <- comparison[i, ]

    # Format numbers nicely
    n_cells_str <- if (row$n_cells > 1e6) {
      sprintf("%.1fM", row$n_cells / 1e6)
    } else if (row$n_cells > 1e3) {
      sprintf("%.1fK", row$n_cells / 1e3)
    } else {
      sprintf("%.0f", row$n_cells)
    }

    cat(sprintf("%-4d  %-12s  %-12.1f  %-12.1f  %-10.1f\n",
                row$resolution,
                n_cells_str,
                row$cell_area_km2,
                row$cell_spacing_km,
                row$cls_km))
  }

  cat(paste(rep("=", 70), collapse = ""), "\n\n")
}

