#' hexify
#'
#' Core icosahedron and 'Snyder' projection helpers.
#'
#' @useDynLib hexify, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods new setClass setMethod setGeneric setValidity is as
#' @importFrom grDevices adjustcolor
#' @importFrom graphics points polygon
#' @importFrom utils head
#' @importFrom rlang .data
"_PACKAGE"

# Global variables to avoid R CMD check notes
utils::globalVariables(c(
  "hexify_world", "cpp_decode_z3", "cpp_decode_zorder",
  "cpp_h3_latLngToCell", "cpp_h3_cellToLatLng", "cpp_h3_isValidCell",
  "cpp_h3_cellToParent", "cpp_h3_cellToChildren", "cpp_h3_cellToBoundary",
  "cpp_h3_polygonToCells", "cpp_h3_cellAreaKm2", "cpp_h3_getResolution",
  "cpp_h3_gridDisk", "cpp_h3_gridDiskDistances", "cpp_h3_gridRingUnsafe",
  "cpp_h3_compactCells", "cpp_h3_uncompactCells",
  "cpp_h3_isPentagon", "cpp_h3_gridDistance",
  "cpp_get_neighbors_isea", "cpp_get_neighbors_z7"
))

#' @importFrom sf st_crs
#' @export
sf::st_crs

#' @importFrom sf st_as_sf
#' @export
sf::st_as_sf
