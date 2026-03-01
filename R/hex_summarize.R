# R/hex_summarize.R
# Cell-level aggregation for HexData objects

#' Summarize Data by Hex Cell
#'
#' Aggregates data within each hexagonal cell, similar to
#' `dplyr::group_by(cell_id) |> summarize(...)`. Returns a data.frame
#' with one row per unique cell, including cell center coordinates and area.
#'
#' @param hex_data A HexData object (from [hexify()]).
#' @param ... Named summary expressions (tidyeval). Each expression is
#'   evaluated per cell group. Examples: `mean_temp = mean(temperature)`,
#'   `n = dplyr::n()`, `max_elev = max(elevation, na.rm = TRUE)`.
#' @param .fns Optional named list of functions for formula-style aggregation.
#'   Example: `.fns = list(mean_temp = ~mean(temperature))`.
#' @param geometry Logical. If `TRUE`, attach cell center points as an
#'   sf geometry column (requires sf). Default `FALSE`.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{cell_id}{Unique cell identifier}
#'     \item{cell_cen_lon, cell_cen_lat}{Cell center coordinates}
#'     \item{cell_area_km2}{Cell area in km^2}
#'     \item{n_points}{Number of data points in this cell}
#'     \item{...}{User-defined summary columns}
#'   }
#'   If `geometry = TRUE`, returns an sf object with POINT geometry.
#'
#' @details
#' The function works entirely in R (no C++ needed). It groups by `cell_id`
#' and evaluates the summary expressions within each group.
#'
#' If no summary expressions are provided, returns cell counts only.
#'
#' @seealso [hexify()] for creating HexData objects,
#'   [get_neighbors()] for finding neighboring cells
#'
#' @export
#' @examples
#' \donttest{
#' df <- data.frame(
#'   lon = runif(100, -10, 10),
#'   lat = runif(100, 40, 55),
#'   temperature = rnorm(100, 15, 5),
#'   species = sample(letters[1:5], 100, replace = TRUE)
#' )
#' hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 500)
#'
#' # Count points per cell
#' hex_summarize(hd)
#'
#' # Custom summaries
#' hex_summarize(hd, mean_temp = mean(temperature),
#'                   n_species = length(unique(species)))
#' }
hex_summarize <- function(hex_data, ..., .fns = NULL, geometry = FALSE) {
  if (!is_hex_data(hex_data)) {
    stop("hex_data must be a HexData object")
  }

  # Build the combined data.frame
  df <- as.data.frame(hex_data)

  # Capture user expressions
  dots <- rlang::enquos(...)

  # Group by cell_id
  cell_ids <- unique(df$cell_id)
  split_df <- split(df, df$cell_id)

  if (length(dots) == 0 && is.null(.fns)) {
    # No summary expressions: just count
    result <- data.frame(
      cell_id = cell_ids,
      n_points = vapply(split_df, nrow, integer(1)),
      stringsAsFactors = FALSE
    )
  } else if (length(dots) > 0) {
    # Tidyeval path
    summaries <- lapply(split_df, function(grp) {
      env <- rlang::as_data_mask(grp)
      vals <- lapply(dots, function(q) rlang::eval_tidy(q, data = env))
      as.data.frame(vals, stringsAsFactors = FALSE)
    })

    summary_df <- do.call(rbind, summaries)
    result <- data.frame(
      cell_id = cell_ids,
      n_points = vapply(split_df, nrow, integer(1)),
      summary_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    # Formula-style (.fns)
    if (!is.list(.fns) || is.null(names(.fns))) {
      stop(".fns must be a named list of functions or formulas")
    }

    # Convert formulas to functions
    fns <- lapply(.fns, function(f) {
      if (inherits(f, "formula")) rlang::as_function(f) else f
    })

    summaries <- lapply(split_df, function(grp) {
      vals <- lapply(fns, function(fn) fn(grp))
      as.data.frame(vals, stringsAsFactors = FALSE)
    })

    summary_df <- do.call(rbind, summaries)
    result <- data.frame(
      cell_id = cell_ids,
      n_points = vapply(split_df, nrow, integer(1)),
      summary_df,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  # Attach cell center coordinates (use first occurrence per cell)
  first_idx <- match(cell_ids, df$cell_id)
  result$cell_cen_lon <- df$cell_cen_lon[first_idx]
  result$cell_cen_lat <- df$cell_cen_lat[first_idx]
  result$cell_area_km2 <- df$cell_area_km2[first_idx]

  # Reorder columns: cell info first, then summaries
  info_cols <- c("cell_id", "cell_cen_lon", "cell_cen_lat", "cell_area_km2",
                  "n_points")
  other_cols <- setdiff(names(result), info_cols)
  result <- result[, c(info_cols, other_cols), drop = FALSE]
  rownames(result) <- NULL


  if (geometry) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required when geometry = TRUE")
    }
    result <- sf::st_as_sf(result,
                            coords = c("cell_cen_lon", "cell_cen_lat"),
                            crs = hex_data@grid@crs)
  }

  result
}
