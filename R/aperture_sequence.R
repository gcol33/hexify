# =============================================================================
# Aperture sequences
# =============================================================================
#
# An ISEA grid refines by one aperture per resolution level. hex_grid() takes
# that as a family name or as a per-level vector:
#
#   "3", "4", "7"    every level refines by that aperture
#   "4/3", "4/7"     the first floor(resolution / 2) levels refine by the first
#                    aperture and the rest by the second, which is how DGGRID
#                    arranges ISEA43H
#   c(4, 4, 7, 3)    one aperture per level, in order; a list shorter than
#                    resolution recurs, as OGC Topic 21 reads its list of
#                    refinement ratios, so c(4, 3) at resolution 5 is 4,3,4,3,4
#
# A grid stores the spelling as a string, a family with "/" ("4/7") and a
# per-level sequence with "," ("4,4,7,3"), so a two-level sequence c(4, 7)
# ("4,7") stays apart from the family "4/7": the two name the same grid at
# resolution 2 but different grids at resolution 1. parse_aperture_seq() turns
# a spelling back into the sequence the C++ layer takes: entry 1 names the base
# grid and the rest are the refinement steps, so its length is resolution + 1.

#' Is this aperture spelling a mixed sequence?
#' @param aperture Character aperture spelling
#' @noRd
is_mixed_aperture <- function(aperture) {
  grepl("[/,]", as.character(aperture))
}

#' Is this aperture spelling one aperture per level?
#' @param aperture Character aperture spelling
#' @noRd
is_per_level_aperture <- function(aperture) {
  grepl(",", as.character(aperture), fixed = TRUE)
}

#' Aperture spelling to store on a grid
#'
#' A family name passes through; a per-level vector, or its comma-separated
#' spelling, is joined with "," at one aperture per level, a shorter list
#' recurring until it covers `resolution` levels.
#' @param aperture Character family name or numeric vector of apertures
#' @param resolution Integer resolution the vector spelling is given for
#' @return Single character string
#' @noRd
format_aperture <- function(aperture, resolution) {
  if (length(aperture) == 1L && is_per_level_aperture(aperture) &&
      !is.null(resolution)) {
    aperture <- suppressWarnings(
      as.integer(strsplit(as.character(aperture), ",", fixed = TRUE)[[1]]))
  }
  if (length(aperture) > 1L) {
    parts <- as.integer(aperture)
    if (anyNA(parts) || !all(parts %in% VALID_APERTURES)) {
      stop(sprintf("Aperture sequence entries must be one of: %s",
                   paste(VALID_APERTURES, collapse = ", ")))
    }
    if (length(parts) > resolution) {
      stop(sprintf(
        "An aperture sequence names one aperture per resolution level: %d given for resolution %d",
        length(parts), as.integer(resolution)
      ))
    }
    return(paste(rep(parts, length.out = resolution), collapse = ","))
  }
  as.character(aperture)
}

#' Aperture of every level of a grid
#'
#' @param aperture Character aperture spelling
#' @param resolution Integer resolution
#' @return Integer vector of length resolution + 1: the base grid followed by
#'   one aperture per refinement step
#' @noRd
parse_aperture_seq <- function(aperture, resolution) {
  resolution <- as.integer(resolution)
  per_level <- is_per_level_aperture(aperture)
  parts <- suppressWarnings(as.integer(strsplit(as.character(aperture), "[/,]")[[1]]))

  if (anyNA(parts) || !all(parts %in% VALID_APERTURES) ||
      (per_level && grepl("/", aperture, fixed = TRUE))) {
    stop(sprintf("Aperture must be one of %s, a family such as \"4/3\", or one aperture per level",
                 paste(VALID_APERTURES, collapse = ", ")))
  }

  if (length(parts) == 1L) {
    steps <- rep(parts, resolution)
  } else if (per_level) {
    if (length(parts) != resolution) {
      stop(sprintf(
        "Aperture \"%s\" names the apertures of %d levels, resolution %d needs %d",
        aperture, length(parts), resolution, resolution
      ))
    }
    steps <- parts
  } else if (length(parts) == 2L) {
    level <- as.integer(resolution / 2)
    steps <- c(rep(parts[1], level), rep(parts[2], resolution - level))
  } else {
    stop(sprintf(
      "Aperture \"%s\": a family names two apertures, such as \"4/3\"; give one aperture per level as \"%s\"",
      aperture, paste(parts, collapse = ",")
    ))
  }

  base <- if (length(steps) > 0L) steps[1] else parts[1]
  as.integer(c(base, steps))
}

#' Levels of an ISEA grid as the C++ entry points take them
#'
#' Every ISEA entry point takes `resolution, aperture, ap_seq`. A pure
#' aperture passes itself and an empty `ap_seq`; a mixed spelling passes
#' aperture 0 and its sequence read at `resolution`.
#' @param aperture Character or numeric aperture spelling
#' @param resolution Integer resolution
#' @return List with integer `resolution`, `aperture` and `ap_seq`
#' @noRd
isea_levels <- function(aperture, resolution) {
  resolution <- as.integer(resolution)
  if (is_mixed_aperture(aperture)) {
    return(list(resolution = resolution, aperture = 0L,
                ap_seq = parse_aperture_seq(aperture_at_resolution(aperture, resolution),
                                            resolution)))
  }
  list(resolution = resolution, aperture = as.integer(aperture), ap_seq = integer(0))
}

#' The same spelling read at a coarser resolution
#'
#' A family name applies at every resolution, so it passes through. A per-level
#' spelling names one aperture per level of its own resolution, so a coarser
#' grid takes the leading levels.
#' @param aperture Character aperture spelling
#' @param resolution Integer resolution to read it at
#' @noRd
aperture_at_resolution <- function(aperture, resolution) {
  if (!is_per_level_aperture(aperture)) return(as.character(aperture))
  parts <- strsplit(as.character(aperture), ",", fixed = TRUE)[[1]]
  if (resolution >= length(parts)) return(as.character(aperture))
  paste(parts[seq_len(max(resolution, 1L))], collapse = ",")
}

#' Cell count of an aperture spelling at a resolution, exact
#'
#' N = d * (product of the refinement apertures) + 2, d the solid's diamond
#' quads, as quad_frame() in src/rcpp_cell.cpp counts the cells it numbers;
#' the largest cell ID. `NA` where N passes 2^63 - 1.
#' @param aperture Character aperture spelling
#' @param resolution Integer resolution
#' @param polyhedron The solid the grid is built on
#' @return integer64
#' @noRd
isea_cell_count <- function(aperture, resolution, polyhedron = "icosahedron") {
  lv <- isea_levels(aperture, resolution)
  cpp_grid_n_cells(standard_icosa(polyhedron), lv$resolution, lv$aperture,
                   lv$ap_seq)
}

#' Cell count of an aperture spelling at a resolution, as a double
#' @param aperture Character aperture spelling
#' @param resolution Integer resolution
#' @param polyhedron The solid the grid is built on
#' @noRd
aperture_n_cells <- function(aperture, resolution, polyhedron = "icosahedron") {
  # Above 2^53 the count rounds to the nearest double, which is all a cell
  # area or a spacing reads from it
  suppressWarnings(as.numeric(isea_cell_count(aperture, resolution, polyhedron)))
}

#' Calculate resolution for target area
#'
#' Uses the cell count formula N = d * aperture^res + 2, d the solid's diamond
#' quads; on the icosahedron (d = 10) this is the 'ISEA3H'/'ISEA4H'/'ISEA7H'
#' count, which matches 'dggridR' resolution numbering exactly.
#'
#' @param target_area_km2 Target area in square kilometers
#' @param aperture Aperture (3, 4, or 7)
#' @param radius_km Radius of the body, in kilometers
#' @param polyhedron The solid the grid is built on
#' @return Resolution level, not rounded. A target larger than the cells of
#'   resolution 0 gives -Inf.
#' @keywords internal
calculate_resolution_for_area <- function(target_area_km2, aperture = 3,
                                          radius_km = EARTH_RADIUS_KM,
                                          polyhedron = "icosahedron") {
  n_cells <- body_surface_km2(radius_km) / target_area_km2
  # Hex9 has 12 * 9^res cells and no vertex cells
  if (is_hex9_aperture(aperture)) {
    return(log(pmax(n_cells / HEX9_BASE_CELLS, 0)) / log(HEX9_APERTURE))
  }
  d <- polyhedron_diamonds(polyhedron)

  # Solving N = d * aperture^res + 2 for res:
  # res = log((surface / area - 2) / d) / log(aperture)
  log(pmax((n_cells - 2) / d, 0)) / log(aperture)
}

#' Resolution whose cells have a target area, for a mixed sequence
#'
#' Cell area is the sphere over the cell count, and for a sequence that count is
#' a product rather than a power, so the resolution comes from the two levels
#' the target falls between, interpolated on the log-area scale that a pure
#' aperture's closed form gives exactly.
#' @param area_km2 Target cell area
#' @param aperture Character aperture spelling
#' @param radius_km Radius of the body, in kilometers
#' @param polyhedron The solid the grid is built on
#' @return Numeric resolution, not rounded
#' @noRd
calculate_resolution_for_area_mixed <- function(area_km2, aperture,
                                                radius_km = EARTH_RADIUS_KM,
                                                polyhedron = "icosahedron") {
  res <- seq.int(MIN_RESOLUTION, isea_max_resolution(aperture, polyhedron))
  log_area <- vapply(res, function(r) {
    log(mean_cell_area_km2(aperture, r, radius_km, polyhedron))
  }, numeric(1))
  target <- log(area_km2)

  if (target >= log_area[1]) return(as.numeric(MIN_RESOLUTION))
  if (target <= log_area[length(log_area)]) return(as.numeric(res[length(res)]))

  k <- max(which(log_area > target))
  frac <- (log_area[k] - target) / (log_area[k] - log_area[k + 1])
  res[k] + frac
}

#' ISEA resolution for a target cell area
#'
#' The resolution whose mean cell area is closest to the target in the
#' direction `round` asks for, clamped to the resolutions whose cell IDs fit
#' in 64 bits (isea_max_resolution()).
#' @param area_km2 Target cell area in km^2
#' @param aperture Aperture spelling, pure or mixed
#' @param radius_km Radius of the body, in kilometers
#' @param round "nearest", "up" (finer cells) or "down" (coarser cells)
#' @param polyhedron The solid the grid is built on
#' @return Resolution, a whole number
#' @noRd
resolve_resolution_from_area <- function(area_km2, aperture,
                                         radius_km = EARTH_RADIUS_KM,
                                         round = "nearest",
                                         polyhedron = "icosahedron") {
  res_exact <- if (is_mixed_aperture(aperture)) {
    calculate_resolution_for_area_mixed(area_km2, aperture, radius_km, polyhedron)
  } else {
    calculate_resolution_for_area(area_km2, as.integer(aperture), radius_km,
                                  polyhedron)
  }

  resolution <- switch(round,
    "nearest" = base::round(res_exact),
    "up" = ceiling(res_exact),
    "down" = floor(res_exact),
    stop("resround must be 'nearest', 'up', or 'down'")
  )

  max(MIN_RESOLUTION, min(isea_max_resolution(aperture, polyhedron), resolution))
}
