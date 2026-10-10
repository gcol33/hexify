# construction.R
# Lambert's construction of Snyder's projection, step by step

#' Points through the steps of Snyder's projection
#'
#' Snyder's equal-area projection of a face is Lambert's azimuthal equal-area
#' projection about the face's centre, adjusted to fill the face's plane
#' triangle. Lambert's projection has a construction with compass and ruler,
#' and this function returns each point at every step of it, on one face.
#'
#' With \eqn{T} the face centre on the sphere and \eqn{S} a point at arc
#' \eqn{z} from it and azimuth \eqn{Az}:
#' \describe{
#'   \item{sphere}{\eqn{S} itself.}
#'   \item{lambert}{\eqn{S} swings about \eqn{T}, in the plane of \eqn{S} and
#'     the normal at \eqn{T}, down onto the tangent plane at \eqn{T},
#'     keeping its distance \eqn{|TS| = 2 \sin(z / 2)}: the Lambert point
#'     \eqn{L}, at azimuth \eqn{Az}. This map is equal-area and takes the
#'     spherical face to a triangle with curved edges.}
#'   \item{nudge}{Within the tangent plane, Snyder's adjustment moves
#'     \eqn{L} to azimuth \eqn{Az'} and radius \eqn{2 f \sin(z / 2)}, which
#'     fill each of the face's six right triangles in the plane with area in
#'     proportion (Snyder 1992, eqs. 6-11).}
#'   \item{snyder}{The tangent plane scales by \eqn{R'} about the sphere's
#'     centre onto the face plane, at distance \eqn{R'} from it; the point
#'     lands on Snyder's point \eqn{P}, at radius \eqn{2 R' f \sin(z / 2)},
#'     and the face on the plane triangle of the face's area.}
#' }
#'
#' Plane coordinates (\code{u}, \code{v}) run from the face centre along the
#' x and y axes of the face's plane triangle, the axes of \code{tx} and
#' \code{ty}, in units of the sphere's radius; an azimuth \eqn{a} lies along
#' \eqn{(\sin a, \cos a)}. Seen from the sphere's centre the planes are
#' parallel, and the point's height along the face's normal is
#' \eqn{\cos z} for \eqn{S}, 1 on the tangent plane and \eqn{R'} on the
#' face plane. The sphere point's \code{u} and \code{v} are its projection
#' onto the planes along the normal.
#'
#' The scale factors at every step come from
#' \code{\link{projection_distortion}} with its \code{stage} argument; they
#' and these points are computed by the same code as the forward projection,
#' and \code{tx}, \code{ty} are those \code{\link{hexify_forward}} returns.
#'
#' @param x A HexGridInfo object from \code{\link{hex_grid}} with Snyder's
#'   projection (\code{projection = "isea"}).
#' @param lon,lat Longitudes and latitudes in degrees (geodetic on the
#'   grid's ellipsoid, if it has one: the construction starts from the
#'   point's place on the authalic sphere).
#' @param face \code{NULL} to take each point onto the face it lies on, or
#'   face numbers from 0, one or one per point; the face's formulas are then
#'   applied as they stand to points off it.
#'
#' @return A data frame with one row per point: \code{lon}, \code{lat},
#'   \code{face} (from 0); \code{arc}, \eqn{z} in degrees; \code{az}
#'   (\eqn{Az}) and \code{az_prime} (\eqn{Az'}) in degrees from the face's
#'   first vertex, clockwise seen from outside the sphere; Snyder's factor
#'   \code{f}; the plane coordinates \code{sphere_u}, \code{sphere_v},
#'   \code{lambert_u}, \code{lambert_v}, \code{nudge_u}, \code{nudge_v},
#'   \code{snyder_u}, \code{snyder_v}; and \code{tx}, \code{ty}, Snyder's
#'   point in triangle coordinates (a face edge is 1). The attribute
#'   \code{"r1"} holds \eqn{R'}.
#'
#' @references
#' Snyder, J. P. (1992). An equal-area map projection for polyhedral globes.
#' \emph{Cartographica} 29(1), 10-21. \doi{10.3138/27H7-8K88-4882-1752}
#'
#' @seealso \code{\link{projection_distortion}};
#'   \code{\link{hex_globe}} with \code{lambert} to follow the construction
#'   on a globe
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' s <- projection_stages(g, c(0, 10, 16.37), c(40, 45, 48.21))
#' s
#'
#' # The Lambert point keeps its distance from the face centre
#' cbind(lambert = sqrt(s$lambert_u^2 + s$lambert_v^2),
#'       chord = 2 * sin(s$arc * pi / 360))
projection_stages <- function(x, lon, lat, face = NULL) {
  g <- face_projection_grid(x, "projection_stages()")
  check_lambert_grid(g)
  check_lonlat_vectors(lon, lat)
  n <- length(lon)
  icosa <- icosa_arg(g)
  if (is.null(face)) {
    face <- rep(NA_integer_, n)
  } else {
    n_faces <- nrow(icosa_solid(icosa)$faces)
    if (!is.numeric(face) || !length(face) %in% c(1L, n) || anyNA(face) ||
        any(face != round(face)) || any(face < 0 | face >= n_faces)) {
      stop("face must be NULL or face numbers from 0 to ", n_faces - 1L,
           ", one or one per point", call. = FALSE)
    }
    face <- rep(as.integer(face), length.out = n)
  }
  s <- cpp_lonlat_construction(icosa, as.numeric(lon), as.numeric(lat), face)
  deg <- 180 / pi
  out <- data.frame(
    lon = as.numeric(lon), lat = as.numeric(lat), face = s$face,
    arc = s$z * deg, az = s$az * deg, az_prime = s$az_prime * deg, f = s$f,
    sphere_u = sin(s$z) * sin(s$az), sphere_v = sin(s$z) * cos(s$az),
    lambert_u = s$lambert_u, lambert_v = s$lambert_v,
    nudge_u = s$plane_u / s$r1, nudge_v = s$plane_v / s$r1,
    snyder_u = s$plane_u, snyder_v = s$plane_v,
    tx = s$tx, ty = s$ty
  )
  attr(out, "r1") <- s$r1
  out
}

# =============================================================================
# DRAWING ONE FACE AT A STEP
# =============================================================================

#' The steps of Lambert's construction whose plane coordinates
#' projection_stages() returns, and "face", the grid's own face projection
#' (on any ISEA grid)
#' @noRd
STAGE_NAMES <- c("sphere", "lambert", "nudge", "snyder", "face")

#' Triangle coordinates of a face as plane coordinates (u, v) from its
#' centre, in sphere radii on the plane triangle with the face's area
#' @noRd
face_plane_uv <- function(g, tx, ty) {
  edge <- face_plane_edge(icosa_arg(g))
  cbind((tx - 0.5) * edge, (ty - 1 / (2 * sqrt(3))) * edge)
}

#' Points given in lon/lat at one step on one face: an n x 2 matrix of the
#' step's plane coordinates (u, v)
#' @noRd
stage_uv <- function(g, face, stage, lon, lat) {
  s <- projection_stages(g, lon, lat, face = face)
  cbind(s[[paste0(stage, "_u")]], s[[paste0(stage, "_v")]])
}

#' Unit vectors (rows) as lon/lat at one step on one face, NA rows kept
#' @noRd
stage_xyz <- function(g, face, stage, P) {
  out <- matrix(NA_real_, nrow(P), 2L)
  ok <- stats::complete.cases(P)
  if (any(ok)) {
    ll <- vec_lonlat(P[ok, , drop = FALSE], icosa_arg(g))
    out[ok, ] <- stage_uv(g, face, stage, ll[, 1], ll[, 2])
  }
  out
}

#' The walls of a grid's cells on one face at one step of Lambert's
#' construction: each cell's run of boundary points on the face, `step`
#' radians apart on the sphere, as plane coordinates (u, v), the runs
#' separated by NA rows for lines() and polygon()
#' @noRd
stage_walls <- function(g, face, stage = STAGE_NAMES, cells = NULL, step = 0.004) {
  stage <- match.arg(stage)
  if (stage != "face") check_lambert_grid(g)
  p <- grid_surface_paths(g, grid_cells(g, cells), step)
  at <- which(p[, "face"] == face)
  if (length(at) == 0L) return(matrix(numeric(0), 0L, 2L))
  p <- p[at, , drop = FALSE]
  n <- nrow(p)
  breaks <- c(FALSE, p[-1L, "cell"] != p[-n, "cell"] | diff(at) > 1L)
  run <- cumsum(breaks)
  rows <- unlist(lapply(split(seq_len(n), run), function(i) c(i, NA_integer_)),
                 use.names = FALSE)
  if (stage == "face") return(face_plane_uv(g, p[rows, "tx"], p[rows, "ty"]))
  P <- p[, c("sphere_x", "sphere_y", "sphere_z"), drop = FALSE]
  stage_xyz(g, face, stage, P[rows, , drop = FALSE])
}

#' The outline of one face at one step: its three edges as great-circle arcs
#' of `n` points each, as plane coordinates (u, v), a closed ring
#' @noRd
stage_outline <- function(g, face, stage = STAGE_NAMES, n = 200L) {
  stage <- match.arg(stage)
  if (stage == "face") return(face_plane_uv(g, FACE_TRIANGLE[, 1], FACE_TRIANGLE[, 2]))
  solid <- icosa_solid(icosa_arg(g))
  V <- solid$vertices / sqrt(rowSums(solid$vertices^2))
  v <- solid$faces[face + 1L, ]
  t <- seq(0, 1, length.out = n)
  ring <- do.call(rbind, lapply(1:3, function(k) {
    a <- V[v[k], ]
    b <- V[v[k %% 3L + 1L], ]
    P <- outer(1 - t, a) + outer(t, b)
    P / sqrt(rowSums(P^2))
  }))
  stage_xyz(g, face, stage, ring)
}

#' Tissot's ellipses on one face at one step (or under the grid's face
#' projection, "face"): at the points of a triangular
#' lattice `steps` to the edge of the plane triangle (corners and edges left
#' out), each the image of a circle of `radius` sphere radii, as plane
#' coordinates (u, v) with an NA row after each ellipse; with the points'
#' angular deformation in degrees as attribute "angular"
#' @noRd
stage_ellipses <- function(g, face, stage = c("lambert", "nudge", "snyder", "face"),
                           steps = 7L, radius = 0.025, n = 49L) {
  stage <- match.arg(stage)
  if (stage != "face") check_lambert_grid(g)
  icosa <- icosa_arg(g)
  ij <- expand.grid(i = seq_len(steps - 1L), j = seq_len(steps - 1L))
  ij <- ij[ij$i + ij$j < steps, ]
  tx <- (ij$i + 0.5 * ij$j) / steps
  ty <- ij$j * sqrt(3) / 2 / steps
  ll <- t(vapply(seq_along(tx), function(k) {
    cpp_face_xy_to_ll(icosa, tx[k], ty[k], as.integer(face))
  }, numeric(2)))
  centre <- if (stage == "face") {
    face_plane_uv(g, tx, ty)
  } else {
    stage_uv(g, face, stage, ll[, 1], ll[, 2])
  }
  code <- CONSTRUCTION_STAGES[[if (stage == "snyder") "face" else stage]]
  s <- cpp_lonlat_tissot(icosa, ll[, 1], ll[, 2], rep(as.integer(face), nrow(ll)), code)
  th <- seq(0, 2 * pi, length.out = n)
  rings <- lapply(seq_len(nrow(centre)), function(k) {
    x <- radius * s$a[k] * cos(th)
    y <- radius * s$b[k] * sin(th)
    rbind(cbind(centre[k, 1] + cos(s$angle[k]) * x - sin(s$angle[k]) * y,
                centre[k, 2] + sin(s$angle[k]) * x + cos(s$angle[k]) * y),
          c(NA, NA))
  })
  out <- do.call(rbind, rings)
  attr(out, "angular") <- angular_deformation(s$a, s$b)
  out
}
