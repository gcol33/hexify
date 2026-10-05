# sphere_geometry.R
# Points on the unit sphere: lon/lat <-> unit vectors, cross products, arcs

#' Unit vectors of lon/lat points
#'
#' @param lon,lat Longitudes and latitudes in degrees
#' @return Matrix with one row (x, y, z) per point
#' @noRd
unit_vec <- function(lon, lat) {
  lon <- lon * pi / 180
  lat <- lat * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

#' Longitude and latitude of unit vectors
#'
#' The z coordinate is clamped to `[-1, 1]` so rounding just past a pole still
#' gives its latitude.
#' @param P Matrix with one row (x, y, z) per point, or one vector of length 3
#' @return Matrix with columns lon and lat in degrees, one row per point
#' @noRd
vec_lonlat <- function(P) {
  P <- matrix(P, ncol = 3L)
  cbind(atan2(P[, 2], P[, 1]), asin(pmax(-1, pmin(1, P[, 3])))) * 180 / pi
}

#' Cross products of matching rows
#'
#' @param a,b Matrices with three columns, or two vectors of length 3
#' @return The row-wise cross products a x b, a vector when both inputs are
#'   vectors
#' @noRd
cross3 <- function(a, b) {
  vectors <- is.null(dim(a)) && is.null(dim(b))
  a <- matrix(a, ncol = 3L)
  b <- matrix(b, ncol = 3L)
  out <- cbind(a[, 2] * b[, 3] - a[, 3] * b[, 2],
               a[, 3] * b[, 1] - a[, 1] * b[, 3],
               a[, 1] * b[, 2] - a[, 2] * b[, 1])
  if (vectors) drop(out) else out
}

#' Points along the great-circle arc a -> b, from a and short of b
#' @noRd
slerp <- function(a, b, max_angle) {
  w <- acos(max(-1, min(1, sum(a * b))))
  n <- max(1L, ceiling(w / max_angle))
  if (w < 1e-12) return(rbind(a))
  s <- (seq_len(n) - 1L) / n
  outer(sin((1 - s) * w), a) / sin(w) + outer(sin(s * w), b) / sin(w)
}
