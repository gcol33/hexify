# crosswalk.R
# Crosswalk between ISEA and H3 grid systems
#
# Maps cell IDs between hexify's native ISEA grids and the H3 hierarchical
# hexagonal system via cell-center coordinate lookup.

#' Crosswalk Between ISEA and H3 Cell IDs
#'
#' Maps cell IDs between ISEA (equal-area) and H3 grid systems by looking up
#' each cell's center coordinate in the target grid. This enables workflows
#' where analysis is done in ISEA (exact equal-area) and reporting in H3
#' (industry-standard).
#'
#' @param cell_id Cell IDs to translate. Numeric for ISEA, character for H3.
#'   When \code{grid} is a HexData object and \code{cell_id} is \code{NULL},
#'   all cell IDs from the data are used.
#' @param grid A HexGridInfo or HexData object. For \code{direction =
#'   "isea_to_h3"}, this must be an ISEA grid. For \code{direction =
#'   "h3_to_isea"}, this must be an H3 grid.
#' @param h3_resolution Target H3 resolution for \code{"isea_to_h3"}, or the
#'   source H3 resolution for \code{"h3_to_isea"}. When \code{NULL} (default),
#'   the closest H3 resolution matching the ISEA cell area is selected
#'   automatically.
#' @param isea_grid A HexGridInfo for the target ISEA grid. Required when
#'   \code{direction = "h3_to_isea"}.
#' @param direction One of \code{"isea_to_h3"} (default) or \code{"h3_to_isea"}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{isea_cell_id}{ISEA cell ID (numeric)}
#'   \item{h3_cell_id}{H3 cell ID (character)}
#'   \item{isea_area_km2}{Area of the ISEA cell in km2}
#'   \item{h3_area_km2}{Geodesic area of the H3 cell in km2}
#'   \item{area_ratio}{Ratio of ISEA area to H3 area}
#' }
#'
#' @details
#' The crosswalk works by computing the center coordinate of each source cell,
#' then finding which cell in the target grid contains that center. This is a
#' many-to-one mapping: multiple ISEA cells may map to the same H3 cell (or
#' vice versa) depending on the relative resolutions.
#'
#' When \code{h3_resolution} is \code{NULL} and \code{direction = "isea_to_h3"},
#' the H3 resolution whose average cell area is closest to the ISEA cell area
#' is chosen automatically. This gives the best 1:1 correspondence.
#'
#' @seealso \code{\link{cell_area}} for per-cell area computation,
#'   \code{\link{hex_grid}} for creating grids
#'
#' @export
#' @examples
#' \donttest{
#' # ISEA -> H3
#' grid <- hex_grid(area_km2 = 1000)
#' cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
#' xwalk <- h3_crosswalk(cells, grid)
#' head(xwalk)
#'
#' # H3 -> ISEA
#' h3 <- hex_grid(resolution = 5, type = "h3")
#' h3_cells <- lonlat_to_cell(c(0, 10), c(45, 50), h3)
#' xwalk2 <- h3_crosswalk(h3_cells, h3, isea_grid = grid, direction = "h3_to_isea")
#' }
h3_crosswalk <- function(cell_id = NULL,
                         grid,
                         h3_resolution = NULL,
                         isea_grid = NULL,
                         direction = c("isea_to_h3", "h3_to_isea")) {

  direction <- match.arg(direction)

  # -------------------------------------------------------------------------
  # Extract grid and cell IDs
  # -------------------------------------------------------------------------
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

  # -------------------------------------------------------------------------
  # Validate direction vs grid type
  # -------------------------------------------------------------------------
  if (direction == "isea_to_h3") {
    if (is_h3_grid(g)) {
      stop("direction = 'isea_to_h3' requires an ISEA grid, got H3")
    }
  } else {
    if (!is_h3_grid(g)) {
      stop("direction = 'h3_to_isea' requires an H3 grid, got ISEA")
    }
    if (is.null(isea_grid)) {
      stop("isea_grid is required when direction = 'h3_to_isea'")
    }
    isea_g <- extract_grid(isea_grid)
    if (is_h3_grid(isea_g)) {
      stop("isea_grid must be an ISEA grid, not H3")
    }
  }

  # -------------------------------------------------------------------------
  # Deduplicate cell IDs
  # -------------------------------------------------------------------------
  unique_ids <- unique(cell_id)

  # -------------------------------------------------------------------------
  # ISEA -> H3
  # -------------------------------------------------------------------------
  if (direction == "isea_to_h3") {
    # Auto-select H3 resolution if not provided
    if (is.null(h3_resolution)) {
      h3_resolution <- closest_h3_resolution(g@area_km2)
    } else {
      h3_resolution <- as.integer(h3_resolution)
      if (h3_resolution < H3_MIN_RESOLUTION || h3_resolution > H3_MAX_RESOLUTION) {
        stop(sprintf("h3_resolution must be between %d and %d",
                     H3_MIN_RESOLUTION, H3_MAX_RESOLUTION))
      }
    }

    # Get ISEA cell centers
    coords <- cell_to_lonlat(unique_ids, g)

    # Convert centers to H3 cell IDs
    h3_ids <- cpp_h3_latLngToCell(coords$lon_deg, coords$lat_deg, h3_resolution)

    # Compute areas
    isea_areas <- rep(g@area_km2, length(unique_ids))
    h3_areas <- cpp_h3_cellAreaKm2(h3_ids)

    data.frame(
      isea_cell_id = unique_ids,
      h3_cell_id = h3_ids,
      isea_area_km2 = isea_areas,
      h3_area_km2 = h3_areas,
      area_ratio = isea_areas / h3_areas,
      stringsAsFactors = FALSE
    )

  } else {
    # -----------------------------------------------------------------------
    # H3 -> ISEA
    # -----------------------------------------------------------------------

    # Get H3 cell centers
    center_df <- cpp_h3_cellToLatLng(as.character(unique_ids))

    # Convert centers to ISEA cell IDs
    isea_ids <- lonlat_to_cell(center_df$lon, center_df$lat, isea_g)

    # Compute areas
    h3_areas <- cpp_h3_cellAreaKm2(as.character(unique_ids))
    isea_areas <- rep(isea_g@area_km2, length(unique_ids))

    data.frame(
      isea_cell_id = isea_ids,
      h3_cell_id = as.character(unique_ids),
      isea_area_km2 = isea_areas,
      h3_area_km2 = h3_areas,
      area_ratio = isea_areas / h3_areas,
      stringsAsFactors = FALSE
    )
  }
}
