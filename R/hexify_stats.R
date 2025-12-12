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
#'   \item{area_km}{Total Earth surface area in km²}
#'   \item{n_cells}{Total number of cells at this resolution}
#'   \item{cell_area_km2}{Average cell area in km²}
#'   \item{cell_spacing_km}{Average distance between cell centers in km}
#'   \item{resolution}{Resolution level}
#'   \item{aperture}{Grid aperture}
#'   
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' stats <- dgearthstat(grid)
#' 
#' print(sprintf("Resolution %d has %.0f cells", 
#'               stats$resolution, stats$n_cells))
#' print(sprintf("Average cell area: %.2f km²", 
#'               stats$cell_area_km2))
#' print(sprintf("Average cell spacing: %.2f km", 
#'               stats$cell_spacing_km))
#' }
dgearthstat <- function(dggs) {
  # Validate grid
  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a hexify_grid object")
  }

  resolution <- get_grid_resolution(dggs)
  
  # Calculate number of cells at this resolution
  # For aperture A at resolution R:
  # - Number of cells = base_faces * A^R
  # - Aperture 3,4: 20 base faces (icosahedron)
  # - Aperture 7: 12 base faces (different topology)
  base_faces <- if (dggs$aperture == 7) 12 else 20
  n_cells <- base_faces * (dggs$aperture ^ resolution)
  
  # Calculate cell area
  cell_area_km2 <- EARTH_SURFACE_KM2 / n_cells

  # Approximate cell spacing (distance between cell centers)
  # For hexagons: spacing ≈ sqrt(area)
  cell_spacing_km <- sqrt(cell_area_km2)

  # Calculate characteristic length scale (CLS)
  # CLS is a measure of cell size that accounts for shape
  # For hexagons, CLS ≈ sqrt(area / (3 * sqrt(3) / 2))
  cls_km <- sqrt(cell_area_km2 / (3 * sqrt(3) / 2))

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

#' Get maximum cell index for a grid
#' 
#' Returns the total number of cells in the grid at the current resolution.
#' This is useful for validating cell indices and understanding grid size.
#' 
#' @param dggs Grid specification from hexify_grid()
#' 
#' @return Number of cells (numeric)
#' 
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' max_cells <- dgmaxcell(grid)
#' print(sprintf("This grid has %.0f cells", max_cells))
#' }
dgmaxcell <- function(dggs) {
  if (!inherits(dggs, "hexify_grid") && !inherits(dggs, "dggs")) {
    stop("dggs must be a hexify_grid object")
  }
  
  stats <- dgearthstat(dggs)
  return(stats$n_cells)
}

#' Find closest resolution for target cell area
#' 
#' Finds the grid resolution that produces cells closest to the target area.
#' This is a helper function for grid construction.
#' 
#' @param dggs Grid specification (aperture and topology must be set)
#' @param area Target cell area in km² (if metric=TRUE)
#' @param round Rounding method ("nearest", "up", "down")
#' @param metric Whether area is in metric units
#' @param show_info Print information about chosen resolution
#' 
#' @return Resolution level (integer)
#' 
#' @export
#' @examples
#' \dontrun{
#' # Create a temporary grid to get aperture settings
#' temp_grid <- list(aperture = 3, topology = "HEXAGON")
#' class(temp_grid) <- "hexify_grid"
#' 
#' # Find resolution for 1000 km² cells
#' res <- dg_closest_res_to_area(temp_grid, area = 1000, 
#'                                metric = TRUE, show_info = TRUE)
#' print(res)
#' }
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
    message(sprintf("  Cell area: %.2f km²", stats$cell_area_km2))
    message(sprintf("  Cell spacing: %.2f km", stats$cell_spacing_km))
    message(sprintf("  Total cells: %.0f", stats$n_cells))
  }
  
  return(resolution)
}

#' Find closest resolution for target cell spacing
#' 
#' Finds the grid resolution that produces cells with spacing (distance
#' between centers) closest to the target spacing.
#' 
#' @param dggs Grid specification (aperture and topology must be set)
#' @param spacing Target cell spacing in km (if metric=TRUE)
#' @param round Rounding method ("nearest", "up", "down")
#' @param metric Whether spacing is in metric units
#' @param show_info Print information about chosen resolution
#' 
#' @return Resolution level (integer)
#' 
#' @export
dg_closest_res_to_spacing <- function(dggs, spacing, round = "nearest",
                                      metric = TRUE, show_info = FALSE) {
  if (!metric) {
    # Convert from miles to km
    spacing <- spacing * MI_TO_KM
  }
  
  # Spacing ≈ sqrt(area), so area ≈ spacing²
  target_area <- spacing ^ 2
  
  return(dg_closest_res_to_area(dggs, target_area, round, 
                                metric = TRUE, show_info))
}

#' Find closest resolution for target CLS
#' 
#' Finds the grid resolution that produces cells with characteristic length
#' scale (CLS) closest to the target CLS.
#' 
#' @param dggs Grid specification (aperture and topology must be set)
#' @param cls Target CLS in km (if metric=TRUE)
#' @param round Rounding method ("nearest", "up", "down")
#' @param metric Whether CLS is in metric units
#' @param show_info Print information about chosen resolution
#' 
#' @return Resolution level (integer)
#' 
#' @export
dg_closest_res_to_cls <- function(dggs, cls, round = "nearest",
                                  metric = TRUE, show_info = FALSE) {
  if (!metric) {
    # Convert from miles to km
    cls <- cls * MI_TO_KM
  }
  
  # For hexagons: CLS = sqrt(area / (3 * sqrt(3) / 2))
  # So: area = CLS² * (3 * sqrt(3) / 2)
  target_area <- cls^2 * (3 * sqrt(3) / 2)
  
  return(dg_closest_res_to_area(dggs, target_area, round,
                                metric = TRUE, show_info))
}

#' Compare grid resolutions
#' 
#' Generates a table comparing different resolution levels for a given
#' grid configuration. Useful for choosing appropriate resolution.
#' 
#' @param aperture Grid aperture (3, 4, or 7)
#' @param res_range Range of resolutions to compare (e.g., 1:10)
#' 
#' @return Data frame with columns: resolution, n_cells, cell_area_km2,
#'   cell_spacing_km, cls_km
#'   
#' @export
#' @examples
#' \dontrun{
#' # Compare resolutions 0-10 for aperture 3
#' comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:10)
#' print(comparison)
#' 
#' # Find resolution with cells ~1000 km²
#' subset(comparison, cell_area_km2 > 900 & cell_area_km2 < 1100)
#' }
hexify_compare_resolutions <- function(aperture = 3, res_range = 0:15) {
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
  
  return(result_df)
}

#' Print resolution comparison table
#' 
#' Pretty-prints a comparison of grid resolutions with human-readable
#' formatting.
#' 
#' @param aperture Grid aperture (3, 4, or 7)
#' @param res_range Range of resolutions to display
#' 
#' @return NULL (prints to console)
#' 
#' @export
hexify_print_resolutions <- function(aperture = 3, res_range = 0:10) {
  comparison <- hexify_compare_resolutions(aperture, res_range)
  
  cat(sprintf("\nGrid Resolution Comparison (Aperture %d)\n", aperture))
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("%-4s  %-12s  %-12s  %-12s  %-10s\n",
              "Res", "# Cells", "Area (km²)", "Spacing (km)", "CLS (km)"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (i in 1:nrow(comparison)) {
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
  
  invisible(NULL)
}
