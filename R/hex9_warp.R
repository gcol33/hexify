# The trained warp of Hex9's own projection (projection = "akw")

#' Where libhex9 publishes its warp field: the file at the commit hexify's
#' Hex9 agreement was checked against (libhex9 2.3.1)
#' @noRd
HEX9_WARP_URL <- paste0(
  "https://raw.githubusercontent.com/MrBenGriffin/libhex9/",
  "1eb0be7d43255b30ecb0c46611fd56b10bf7d773/core/Sphere_l6_fund.f64g.h9warp"
)

#' Name of the warp field file
#' @noRd
HEX9_WARP_FILE <- "Sphere_l6_fund.f64g.h9warp"

#' Where hexify reads the warp field: the option \code{hexify.hex9_warp}, or
#' hexify's user data directory
#' @noRd
hex9_warp_path <- function() {
  getOption("hexify.hex9_warp",
            file.path(tools::R_user_dir("hexify", "data"), HEX9_WARP_FILE))
}

#' Load the warp field once per session, before a grid on projection "akw"
#' reaches the C++ layer
#' @noRd
ensure_hex9_warp <- function() {
  if (cpp_hex9_warp_ready()) return(invisible(TRUE))
  path <- hex9_warp_path()
  if (!file.exists(path)) {
    stop("projection = \"akw\" reads the trained warp field of 'libhex9' ",
         "(19 MB), which hexify does not ship. Fetch it once with ",
         "hex9_warp_download(), or set options(hexify.hex9_warp = ) to a copy of ",
         "libhex9's core/", HEX9_WARP_FILE, ".", call. = FALSE)
  }
  cpp_hex9_warp_load(readBin(path, "raw", file.info(path)$size))
  invisible(TRUE)
}

#' Download Hex9's warp field
#'
#' Hex9's own projection (\code{hex_grid(aperture = 9, projection = "akw")})
#' is Kaseorg's octahedral projection with an area-correcting warp that Griffin
#' trained by optimal transport and publishes with 'libhex9', his reference
#' implementation (Griffin 2026, Sec. 11b). The field is 19 MB of stored
#' displacements and gradients, which hexify reads but does not ship; this
#' function downloads it once from the 'libhex9' repository, checks it, and
#' saves it where hexify looks for it.
#'
#' @param path Where to save the field. The default is hexify's user data
#'   directory (\code{tools::R_user_dir("hexify", "data")}), or the option
#'   \code{hexify.hex9_warp} when set; hexify reads it from the same place.
#' @param quiet Passed to \code{utils::download.file()}.
#'
#' @return The path of the saved file, invisibly.
#'
#' @details
#' The file is 'libhex9''s \code{core/Sphere_l6_fund.f64g.h9warp} at the
#' commit of version 2.3.1, which the Hex9 agreement checks in hexify were run
#' against. 'libhex9' is Apache-2.0 licensed. Its field is defined on the unit
#' sphere: 'libhex9''s WGS84 functions first turn geodetic latitude into
#' authalic latitude, while hexify reads latitude as spherical latitude, as
#' 'libhex9''s \code{_sphere} functions do.
#'
#' @references
#' Griffin, B. (2026). Hex9: A Quasi-Authalic, Quasi-Continuous Hexagonal
#' DGGS on the Reference Ellipsoid. arXiv:2608.00022.
#'
#' @export
#' @examples
#' if (interactive()) {
#'   hex9_warp_download()
#'   g <- hex_grid(resolution = 5, aperture = 9, projection = "akw")
#'   lonlat_to_cell(16.37, 48.21, g)
#' }
hex9_warp_download <- function(path = hex9_warp_path(), quiet = FALSE) {
  tmp <- tempfile(fileext = ".h9warp")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(HEX9_WARP_URL, tmp, mode = "wb", quiet = quiet)
  cpp_hex9_warp_load(readBin(tmp, "raw", file.info(tmp)$size))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(tmp, path, overwrite = TRUE)) {
    stop("could not write the warp field to ", path, call. = FALSE)
  }
  invisible(path)
}
