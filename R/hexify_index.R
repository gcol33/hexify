# hexify_index.R
# Unified cell indexing: space-filling curves (Z3, Z7, Z-order) and sequential numbers
#
# This module provides all cell identification and indexing operations.
# Each cell can be identified by:
#   1. Coordinates: (face, i, j, resolution)
#   2. Index string: hierarchical encoding (Z3, Z7, Z-order)
#   3. Seqnum: integer sequential number
#
# @name hexify-indexing

# =============================================================================
# INDEX STRING OPERATIONS
# =============================================================================

#' Convert cell coordinates to index string
#'
#' Converts DGGRID cell coordinates (face, i, j) to a hierarchical index string.
#' The index type is automatically selected based on aperture unless specified.
#'
#' @param face Face/quad number (0-19)
#' @param i I coordinate
#' @param j J coordinate
#' @param resolution Resolution level
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding: "auto" (default), "z3", "z7", or "zorder"
#'
#' @return Index string (e.g., "051223")
#'
#' @details
#' Default index types by aperture:
#' - Aperture 3: Z3 (optimized digit selection)
#' - Aperture 4: Z-order (Morton curve)
#' - Aperture 7: Z7 (hierarchical with Class III handling)
#'
#' @family hierarchical index
#' @keywords internal
#' @export
#' @examples
#' idx <- hexify_cell_to_index(5, 10, 15, resolution = 3, aperture = 3)
hexify_cell_to_index <- function(face, i, j, resolution, aperture = 3L,
                                  index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_cell_to_index(
    as.integer(face), as.numeric(i), as.numeric(j),
    as.integer(resolution), as.integer(aperture), index_type
  )
}

#' Convert index string to cell coordinates
#'
#' Decodes a hierarchical index string back to cell coordinates.
#'
#' @param index Index string
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding used. Default "auto" infers from aperture.
#'
#' @return List with face, i, j, and resolution
#'
#' @family hierarchical index
#' @keywords internal
#' @export
#' @examples
#' cell <- hexify_index_to_cell("0012012", aperture = 3)
hexify_index_to_cell <- function(index, aperture = 3L,
                                  index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_index_to_cell(as.character(index), as.integer(aperture), index_type)
}

#' Convert longitude/latitude to index string
#'
#' Main entry point for geocoding points to grid cells.
#'
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#' @param resolution Resolution level
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding: "auto" (default), "z3", "z7", or "zorder"
#'
#' @return Index string
#'
#' @family hierarchical index
#' @keywords internal
#' @export
#' @examples
#' idx <- hexify_lonlat_to_index(16.37, 48.21, resolution = 5, aperture = 3)
hexify_lonlat_to_index <- function(lon, lat, resolution, aperture = 3L,
                                    index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_lonlat_to_index(lon, lat, as.integer(resolution),
                      as.integer(aperture), index_type)
}

#' Convert index string to longitude/latitude
#'
#' Returns the cell center coordinates for a given index.
#'
#' @param index Index string
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding. Default "auto" infers from aperture.
#'
#' @return Named numeric vector with lon and lat in degrees
#'
#' @family hierarchical index
#' @keywords internal
#' @export
#' @examples
#' coords <- hexify_index_to_lonlat("0012012", aperture = 3)
hexify_index_to_lonlat <- function(index, aperture = 3L,
                                    index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_index_to_lonlat(as.character(index), as.integer(aperture), index_type)
}

# =============================================================================
# INDEX HIERARCHY OPERATIONS
# =============================================================================

#' Get parent index
#'
#' Returns the parent index (one resolution coarser).
#'
#' @param index Index string
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding. Default "auto".
#'
#' @return Parent index string
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_get_parent <- function(index, aperture = 3L,
                               index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_get_parent_index(as.character(index), as.integer(aperture), index_type)
}

#' Get children indices
#'
#' Returns all children indices (one resolution finer).
#'
#' @param index Index string
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding. Default "auto".
#'
#' @return Character vector of child indices
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_get_children <- function(index, aperture = 3L,
                                 index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_get_children_indices(as.character(index), as.integer(aperture), index_type)
}

#' Get index resolution
#'
#' Returns the resolution level encoded in an index string.
#'
#' @param index Index string
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index encoding. Default "auto".
#'
#' @return Integer resolution level
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_get_resolution <- function(index, aperture = 3L,
                                   index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_get_index_resolution(as.character(index), as.integer(aperture), index_type)
}

#' Compare two indices
#'
#' Lexicographic comparison of two index strings.
#'
#' @param idx1 First index string
#' @param idx2 Second index string
#'
#' @return Integer: -1 if idx1 < idx2, 0 if equal, 1 if idx1 > idx2
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_compare_indices <- function(idx1, idx2) {
  cpp_compare_indices(as.character(idx1), as.character(idx2))
}

# =============================================================================
# INTEGER CELL IDs
# =============================================================================

#' Convert longitude/latitude to cell ID
#'
#' Converts geographic coordinates to DGGRID-compatible cell identifiers.
#' This is the primary function for geocoding points to grid cells.
#'
#' @param lon Numeric vector of longitudes in degrees
#' @param lat Numeric vector of latitudes in degrees
#' @param resolution Grid resolution (integer >= 0)
#' @param aperture Grid aperture (3, 4, or 7)
#'
#' @return Numeric vector of cell IDs (1-based)
#'
#' @details
#' Returns DGGRID-compatible cell identifiers. The cell ID
#' uniquely identifies each hexagonal cell in the global grid.
#'
#' For a grid-object interface, use \code{\link{lonlat_to_cell}}.
#'
#' @family coordinate conversion
#' @seealso \code{\link{lonlat_to_cell}} for the recommended S4 interface,
#'   \code{\link{hexify_cell_to_lonlat}} for the inverse operation
#'
#' @keywords internal
#' @export
#' @examples
#' cell_id <- hexify_lonlat_to_cell(0, 45, resolution = 5, aperture = 3)
hexify_lonlat_to_cell <- function(lon, lat, resolution, aperture) {
  validate_lon(lon)
  validate_lat(lat)
  validate_resolution(resolution)
  validate_aperture(aperture)
  cpp_lonlat_to_cell(lon, lat, resolution, aperture)
}

#' Convert cell ID to longitude/latitude
#'
#' Converts cell identifiers back to cell center coordinates.
#' This is the inverse of \code{\link{hexify_lonlat_to_cell}}.
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution (integer >= 0)
#' @param aperture Grid aperture (3, 4, or 7)
#'
#' @return Data frame with lon_deg and lat_deg columns
#'
#' @family coordinate conversion
#' @seealso \code{\link{cell_to_lonlat}} for the recommended S4 interface,
#'   \code{\link{hexify_lonlat_to_cell}} for the forward operation
#'
#' @keywords internal
#' @export
#' @examples
#' coords <- hexify_cell_to_lonlat(1702, resolution = 5, aperture = 3)
hexify_cell_to_lonlat <- function(cell_id, resolution, aperture) {
  validate_resolution(resolution)
  validate_aperture(aperture)
  validate_cell_id(cell_id, resolution, aperture)
  cpp_cell_to_lonlat(cell_id, resolution, aperture)
}

#' Get cell info from cell ID
#'
#' Converts cell ID to cell components (quad, i, j).
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution (integer >= 0)
#' @param aperture Grid aperture (3, 4, or 7)
#'
#' @return Data frame with quad, i, j columns
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' info <- hexify_cell_to_quad_ij(1702, resolution = 5, aperture = 3)
hexify_cell_id_to_quad_ij <- function(cell_id, resolution, aperture) {
  validate_resolution(resolution)
  validate_aperture(aperture)
  validate_cell_id(cell_id, resolution, aperture)
  result <- cpp_cell_to_quad_ij(cell_id, resolution, aperture)
  as.data.frame(result)
}

# =============================================================================
# Z7 SPECIAL OPERATIONS
# =============================================================================

#' Get canonical form of Z7 index
#'
#' For Z7 indices that form cycles during decode/encode, returns the
#' lexicographically smallest index in the cycle. Provides stable
#' unique identifiers for aperture 7 grids.
#'
#' @param index Z7 index string
#' @param max_iterations Maximum iterations for cycle detection (default 128)
#'
#' @return Canonical form (lexicographically smallest in cycle)
#'
#' @family hierarchical index
#' @keywords internal
#' @export
#' @examples
#' # These all return the same canonical form
#' hexify_z7_canonical("110001")
#' hexify_z7_canonical("110002")
hexify_z7_canonical <- function(index, max_iterations = 128L) {
  cpp_z7_canonical_form(as.character(index), as.integer(max_iterations))
}

# =============================================================================
# INDEX TYPE UTILITIES
# =============================================================================

#' Check if index type is valid for aperture
#'
#' @param aperture Aperture (3, 4, or 7)
#' @param index_type Index type to check
#'
#' @return Logical: TRUE if valid combination
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_is_valid_index_type <- function(aperture,
                                        index_type = c("auto", "z3", "z7", "zorder")) {
  index_type <- match.arg(index_type)
  cpp_is_valid_index_type(as.integer(aperture), index_type)
}

#' Get default index type for aperture
#'
#' Returns the recommended index type:
#' - Aperture 3: "z3"
#' - Aperture 4: "zorder"
#' - Aperture 7: "z7"
#'
#' @param aperture Aperture (3, 4, or 7)
#'
#' @return String: "z3", "z7", or "zorder"
#'
#' @family hierarchical index
#' @keywords internal
#' @export
hexify_default_index_type <- function(aperture) {
  cpp_get_default_index_type(as.integer(aperture))
}

# =============================================================================
# RESOLUTION CONVERSION (ISEA3H-specific)
# =============================================================================

#' Convert effective resolution to area
#'
#' Calculates approximate cell area for ISEA3H (aperture 3) grids.
#' Based on calibration: eff_res 10 = 863.8006 km^2.
#'
#' @param eff_res Effective resolution (0.5 * resolution for aperture 3)
#'
#' @return Area in km^2
#'
#' @family grid statistics
#' @keywords internal
#' @export
hexify_eff_res_to_area <- function(eff_res) {
  ISEA3H_RES10_AREA_KM2 * (3^(10 - eff_res))
}

#' Convert area to effective resolution
#'
#' Finds the effective resolution that gives approximately the target area.
#'
#' @param area_km2 Target area in km^2
#'
#' @return Effective resolution (may be fractional)
#'
#' @family grid statistics
#' @keywords internal
#' @export
hexify_area_to_eff_res <- function(area_km2) {
  10 - log(area_km2 / ISEA3H_RES10_AREA_KM2) / log(3)
}

#' Convert effective resolution to index resolution
#'
#' For aperture 3: res_index = 2 * eff_res - 1 for odd index resolutions.
#'
#' @param eff_res Effective resolution
#'
#' @return Index resolution (integer)
#'
#' @family grid statistics
#' @keywords internal
#' @export
hexify_eff_res_to_resolution <- function(eff_res) {

  as.integer(2 * eff_res - 1)
}

#' Convert index resolution to effective resolution
#'
#' @param resolution Index resolution
#'
#' @return Effective resolution
#'
#' @family grid statistics
#' @keywords internal
#' @export
hexify_resolution_to_eff_res <- function(resolution) {
  (resolution + 1) / 2
}
