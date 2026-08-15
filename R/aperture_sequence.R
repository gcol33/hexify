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
#   c(4, 4, 7, 3)    one aperture per level, in order, length = resolution
#
# A grid stores the spelling as a string ("4/7", "4/4/7/3") and
# parse_aperture_seq() turns it back into the sequence the C++ layer takes:
# entry 1 names the base grid and the rest are the refinement steps, so its
# length is resolution + 1.

#' Is this aperture spelling a mixed sequence?
#' @param aperture Character aperture spelling
#' @noRd
is_mixed_aperture <- function(aperture) {
  grepl("/", as.character(aperture), fixed = TRUE)
}

#' Aperture spelling to store on a grid
#'
#' A family name passes through; a per-level vector is joined with "/".
#' @param aperture Character family name or numeric vector of apertures
#' @param resolution Integer resolution the vector spelling is given for
#' @return Single character string
#' @noRd
format_aperture <- function(aperture, resolution) {
  if (length(aperture) > 1L) {
    parts <- as.integer(aperture)
    if (anyNA(parts) || !all(parts %in% VALID_APERTURES)) {
      stop(sprintf("Aperture sequence entries must be one of: %s",
                   paste(VALID_APERTURES, collapse = ", ")))
    }
    if (length(parts) != resolution) {
      stop(sprintf(
        "An aperture sequence needs one aperture per resolution level: %d given for resolution %d",
        length(parts), as.integer(resolution)
      ))
    }
    return(paste(parts, collapse = "/"))
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
  parts <- as.integer(strsplit(as.character(aperture), "/", fixed = TRUE)[[1]])

  if (anyNA(parts) || !all(parts %in% VALID_APERTURES)) {
    stop(sprintf("Aperture must be one of %s, a family such as \"4/3\", or one aperture per level",
                 paste(VALID_APERTURES, collapse = ", ")))
  }

  if (length(parts) == 1L) {
    steps <- rep(parts, resolution)
  } else if (length(parts) == resolution) {
    steps <- parts
  } else if (length(parts) == 2L) {
    level <- as.integer(resolution / 2)
    steps <- c(rep(parts[1], level), rep(parts[2], resolution - level))
  } else {
    stop(sprintf(
      "Aperture \"%s\" names %d levels, resolution %d needs %d",
      aperture, length(parts), resolution, resolution
    ))
  }

  base <- if (length(steps) > 0L) steps[1] else parts[1]
  as.integer(c(base, steps))
}

#' Aperture sequence of a grid object
#' @param g HexGridInfo object
#' @noRd
grid_ap_seq <- function(g) {
  parse_aperture_seq(g@aperture, g@resolution)
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
  parts <- strsplit(as.character(aperture), "/", fixed = TRUE)[[1]]
  if (length(parts) <= 2L || resolution >= length(parts)) return(as.character(aperture))
  paste(parts[seq_len(max(resolution, 1L))], collapse = "/")
}

#' Cell count of an aperture sequence
#'
#' N = 10 * (product of the refinement apertures) + 2, the formula
#' calc_grid_params_mixed() in src/rcpp_cell.cpp packs cell IDs for.
#' @param ap_seq Integer aperture sequence
#' @noRd
ap_seq_n_cells <- function(ap_seq) {
  10 * prod(as.numeric(ap_seq[-1])) + 2
}

#' Cell count of an aperture spelling at a resolution
#' @param aperture Character aperture spelling
#' @param resolution Integer resolution
#' @noRd
aperture_n_cells <- function(aperture, resolution) {
  if (is_mixed_aperture(aperture)) {
    ap_seq_n_cells(parse_aperture_seq(aperture, resolution))
  } else {
    max_cell_id(resolution, as.integer(aperture))
  }
}

#' Resolution whose cells have a target area, for a mixed sequence
#'
#' Cell area is the sphere over the cell count, and for a sequence that count is
#' a product rather than a power, so the resolution comes from the two levels
#' the target falls between, interpolated on the log-area scale that a pure
#' aperture's closed form gives exactly.
#' @param area_km2 Target cell area
#' @param aperture Character aperture spelling
#' @return Numeric resolution, not rounded
#' @noRd
calculate_resolution_for_area_mixed <- function(area_km2, aperture) {
  res <- seq.int(MIN_RESOLUTION, MAX_RESOLUTION)
  log_area <- vapply(res, function(r) {
    log(EARTH_SURFACE_KM2 / aperture_n_cells(aperture, r))
  }, numeric(1))
  target <- log(area_km2)

  if (target >= log_area[1]) return(as.numeric(MIN_RESOLUTION))
  if (target <= log_area[length(log_area)]) return(as.numeric(MAX_RESOLUTION))

  k <- max(which(log_area > target))
  frac <- (log_area[k] - target) / (log_area[k] - log_area[k + 1])
  res[k] + frac
}
