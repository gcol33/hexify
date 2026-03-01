# R/h3_compat.R
# H3 interoperability layer
#
# import_h3() — ingest external H3 cell IDs into hexify

# =============================================================================
# IMPORT FROM EXTERNAL H3
# =============================================================================

#' Import External H3 Cell IDs into hexify
#'
#' Ingests H3 cell IDs from an external source (another H3 library, a database,
#' or a CSV file) into hexify. Validates cell IDs, infers the H3 resolution,
#' and optionally attaches data to build a HexData object.
#'
#' For converting between grid specs, use
#' \code{\link{hex_grid}(type = "h3")} directly. For cell-level ISEA/H3
#' mapping, use \code{\link{h3_crosswalk}}.
#'
#' @param cell_ids Character vector of H3 cell ID strings
#' @param data Optional data frame to attach. Must have the same number of rows
#'   as \code{cell_ids}. If \code{NULL}, returns a HexGridInfo object.
#' @param validate If \code{TRUE} (default), checks that all cell IDs are valid
#'   H3 cells at the same resolution before proceeding. Set to \code{FALSE} to
#'   skip validation when cell IDs are known to be correct.
#'
#' @return If \code{data = NULL}, a HexGridInfo object for the inferred H3
#'   resolution. If \code{data} is provided, a HexData object with data
#'   attached at the specified cells.
#'
#' @details
#' H3 cell IDs encode their resolution in the index itself, so no resolution
#' argument is needed. The resolution is inferred automatically. All cell IDs
#' must share the same resolution; mixed resolutions produce an error.
#'
#' @seealso \code{\link{hex_grid}} for creating grids directly,
#'   \code{\link{h3_crosswalk}} for cell-level ISEA/H3 mapping
#'
#' @export
#' @examples
#' \donttest{
#' # Import external H3 cell IDs (grid spec only)
#' h3_ids <- c("8528342bfffffff", "85283473fffffff", "85283447fffffff")
#' grid <- import_h3(h3_ids)
#' grid
#'
#' # Import with data attached
#' df <- data.frame(species = c("oak", "pine", "birch"), count = c(10, 5, 3))
#' hd <- import_h3(h3_ids, data = df)
#' hd
#' }
import_h3 <- function(cell_ids, data = NULL, validate = TRUE) {

  if (!is.character(cell_ids)) {
    stop("cell_ids must be a character vector of H3 cell ID strings")
  }
  if (length(cell_ids) == 0) {
    stop("cell_ids must not be empty")
  }

  # Track non-NA positions
  non_na <- !is.na(cell_ids)
  valid_ids <- cell_ids[non_na]

  if (length(valid_ids) == 0) {
    stop("cell_ids contains only NA values")
  }

  # Validate cell IDs
  if (validate) {
    is_valid <- cpp_h3_isValidCell(valid_ids)
    if (any(!is_valid)) {
      n_invalid <- sum(!is_valid)
      first_bad <- valid_ids[which(!is_valid)[1]]
      stop(sprintf(
        "%d invalid H3 cell ID(s). First invalid: '%s'",
        n_invalid, first_bad
      ))
    }
  }

  # Infer resolution (always needed, even when validate = FALSE)
  resolutions <- cpp_h3_getResolution(valid_ids)
  unique_res <- unique(resolutions)

  if (length(unique_res) > 1) {
    stop(sprintf(
      "All cell_ids must share the same H3 resolution. Found resolutions: %s",
      paste(sort(unique_res), collapse = ", ")
    ))
  }

  res <- unique_res[1L]

  # Build HexGridInfo
  grid <- hex_grid(resolution = res, type = "h3")

  # If no data, return grid only
  if (is.null(data)) {
    return(grid)
  }

  # Validate data dimensions
  if (!inherits(data, "data.frame")) {
    stop(sprintf("data must be a data.frame, got %s", class(data)[1]))
  }
  if (nrow(data) != length(cell_ids)) {
    stop(sprintf(
      "data has %d rows but cell_ids has %d elements; they must match",
      nrow(data), length(cell_ids)
    ))
  }

  # Get cell centers for HexData
  centers <- cpp_h3_cellToLatLng(cell_ids)
  cell_center <- cbind(lon = centers$lon, lat = centers$lat)

  new_hex_data(
    data = data,
    grid = grid,
    cell_id = cell_ids,
    cell_center = cell_center
  )
}
