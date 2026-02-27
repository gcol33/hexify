# R/dggrid.R - DGGRID helper functions
#
# Note: C++ bindings are auto-generated in RcppExports.R.
# This file contains only pure R helper functions.

#' Create DGGRID 43H aperture sequence
#'
#' Create an aperture sequence following DGGRID's 43H pattern:
#' first num_ap4 resolutions use aperture 4, then aperture 3.
#'
#' @param num_ap4 Number of aperture-4 resolutions
#' @param num_ap3 Number of aperture-3 resolutions
#' @return Integer vector of aperture sequence
#'
#' @family 'dggridR' compatibility
#' @keywords internal
#' @export
#' @examples
#' # DGGRID 43H with 2 ap4 resolutions, then 3 ap3 resolutions
#' seq <- dggrid_43h_sequence(2, 3)  # c(4, 4, 3, 3, 3)
dggrid_43h_sequence <- function(num_ap4, num_ap3) {
  c(rep(4L, num_ap4), rep(3L, num_ap3))
}
