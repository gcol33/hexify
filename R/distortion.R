# distortion.R
# Distortion of the face projection: Tissot's indicatrix of the map from the
# sphere onto the faces of the solid

#' Distortion of a grid's face projection
#'
#' Tissot's indicatrix of the projection that maps the sphere onto the faces
#' of the grid's solid (Snyder's or the vertex-oriented equal-area projection,
#' or Fuller's): a small
#' circle on the sphere maps to an ellipse on the face, with semi-axes
#' \code{a} and \code{b} times the circle's radius.
#'
#' The scale factors are the singular values of the projection's derivative,
#' computed exactly by forward-mode automatic differentiation of the
#' projection code. Face-plane lengths are measured on the plane triangle
#' that has the face's area, so the areal scale \code{a * b} is 1 everywhere
#' under Snyder's projection and IVEA, and its mean over a face is 1 under
#' Fuller's.
#'
#' Snyder's projection bends along the arcs from each face's centre to its
#' corners, where its derivative jumps, and IVEA along those and the arcs to
#' the edge midpoints; there and on face edges the result is the derivative
#' on one side. Its angular deformation is largest at the
#' face centre, approached along those arcs: 17.27 degrees, with scale
#' factors 1.163 and 0.860 (Snyder 1992, Table 1). At a face centre the
#' projection has a different derivative along each direction, and the
#' result is the limit along the azimuth the point's coordinates give.
#'
#' \code{stage} reads the indicatrix at an earlier step of Lambert's
#' construction of Snyder's projection (see \code{\link{projection_stages}}):
#' \code{"lambert"}, Lambert's azimuthal equal-area projection about the face
#' centre onto the tangent plane, with scale factors
#' \eqn{1 / \cos(z / 2)} across the radius and \eqn{\cos(z / 2)} along it at
#' arc \eqn{z} from the centre; \code{"nudge"}, Snyder's point on the tangent
#' plane before the scaling onto the face plane, whose scale factors are those
#' of \code{"face"} over \eqn{R'}, so its areal scale is \eqn{1 / R'^2}
#' everywhere. Both read lengths on the tangent plane in units of the
#' sphere's radius.
#'
#' @param x A HexGridInfo object from \code{\link{hex_grid}} with an ISEA
#'   grid. H3 is built on a gnomonic projection of its own and is not
#'   covered.
#' @param lon,lat Longitudes and latitudes in degrees.
#' @param stage \code{"face"}, the grid's face projection, or a step of
#'   Lambert's construction, \code{"lambert"} or \code{"nudge"}, on a grid
#'   with Snyder's projection (\code{projection = "isea"}).
#'
#' @return A data frame with one row per point: \code{lon}, \code{lat},
#'   \code{face} (from 0), the scale factors \code{a} (largest) and \code{b}
#'   (smallest), the maximum angular deformation \code{angular}
#'   \eqn{2 \arcsin((a - b) / (a + b))} in degrees, the areal scale
#'   \code{areal} \eqn{a b}, and \code{angle}, the direction of the longer
#'   axis on the face plane in degrees anticlockwise from the face's x axis.
#'
#' @references
#' Snyder, J. P. (1992). An equal-area map projection for polyhedral globes.
#' \emph{Cartographica} 29(1), 10-21. \doi{10.3138/27H7-8K88-4882-1752}
#'
#' van Leeuwen, D. and Strebe, D. (2006). A "slice-and-dice" approach to area
#' equivalence in polyhedral map projections. \emph{Cartography and
#' Geographic Information Science} 33(4), 269-286.
#' \doi{10.1559/152304006779500687}
#'
#' @seealso \code{\link[=plot,HexGridInfo,missing-method]{plot}} with
#'   \code{distortion} and \code{tissot} to map it;
#'   \code{\link{projection_stages}} for the points at each step of Lambert's
#'   construction
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' projection_distortion(g, c(0, 16.37), c(0, 48.21))
#' projection_distortion(g, c(0, 16.37), c(0, 48.21), stage = "lambert")
#'
#' gf <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
#' d <- projection_distortion(gf, runif(1000, -180, 180),
#'                            asin(runif(1000, -1, 1)) * 180 / pi)
#' range(d$areal)
projection_distortion <- function(x, lon, lat, stage = c("face", "lambert", "nudge")) {
  stage <- match.arg(stage)
  g <- face_projection_grid(x, "projection_distortion()")
  check_lonlat_vectors(lon, lat)
  if (stage != "face") check_lambert_grid(g)
  s <- cpp_lonlat_tissot(icosa_arg(g), as.numeric(lon), as.numeric(lat),
                         rep(NA_integer_, length(lon)), CONSTRUCTION_STAGES[[stage]])
  tissot_table(data.frame(lon = as.numeric(lon), lat = as.numeric(lat)), s)
}

#' The steps of Lambert's construction the C++ layer reads
#' (\code{hexify::ConstructionStage}): the Lambert point, the nudged point,
#' the face projection
#' @noRd
CONSTRUCTION_STAGES <- c(lambert = 0L, nudge = 1L, face = 2L)

#' The grid of an object whose ISEA face projection a function reads; H3 is
#' refused
#' @noRd
face_projection_grid <- function(x, what) {
  g <- extract_grid(x)
  if (is_h3_grid(g)) {
    stop("an H3 grid is built on a gnomonic projection of its own; ",
         what, " reads the ISEA face projections", call. = FALSE)
  }
  g
}

#' Stops unless the grid's face projection is Lambert's construction
#' @noRd
check_lambert_grid <- function(g) {
  if (!is_lambert_projection(grid_projection(g))) {
    stop("Lambert's construction is the one behind Snyder's projection; ",
         "the grid needs projection = \"isea\"", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stops unless lon and lat are finite numbers, as many of each
#' @noRd
check_lonlat_vectors <- function(lon, lat) {
  if (!is.numeric(lon) || !is.numeric(lat) || length(lon) != length(lat)) {
    stop("lon and lat must be numeric vectors of the same length", call. = FALSE)
  }
  if (any(!is.finite(lon) | !is.finite(lat))) {
    stop("lon and lat must be finite", call. = FALSE)
  }
  invisible(TRUE)
}

#' Points with their indicatrix: the scale factors from cpp_*_tissot() and
#' the measures read from them
#' @noRd
tissot_table <- function(points, s) {
  cbind(points, data.frame(face = s$face, a = s$a, b = s$b,
                           angular = angular_deformation(s$a, s$b),
                           areal = s$a * s$b, angle = s$angle * 180 / pi))
}

#' Maximum angular deformation in degrees from the scale factors
#' @noRd
angular_deformation <- function(a, b) 2 * asin((a - b) / (a + b)) * 180 / pi

#' The edge of a face's plane triangle in units of the unit sphere: the
#' triangle with the face's share of the sphere's area
#' @noRd
face_plane_edge <- function(icosa) {
  n_faces <- nrow(icosa_solid(icosa)$faces)
  sqrt(16 * pi / (sqrt(3) * n_faces))
}

# =============================================================================
# DRAWING
# =============================================================================

#' Check the distortion and tissot arguments of the plot method
#' @noRd
resolve_distortion <- function(distortion, tissot, surface, g) {
  if ((distortion != "none" || !isFALSE(tissot)) && is_h3_grid(g)) {
    stop("distortion and tissot map the ISEA face projection; an H3 grid ",
         "is built on a projection of its own", call. = FALSE)
  }
  if (!isFALSE(tissot) && surface == "sphere") {
    stop("Tissot's indicatrix is drawn on the face plane; use surface = ",
         "\"net\" or \"solid\"", call. = FALSE)
  }
  if (isFALSE(tissot)) return(NULL)
  if (isTRUE(tissot)) return(30)
  if (!is.numeric(tissot) || length(tissot) != 1L || !is.finite(tissot) ||
      tissot <= 0 || tissot > 90) {
    stop("tissot must be TRUE, FALSE, or degrees between ellipses (up to 90)",
         call. = FALSE)
  }
  tissot
}

#' The faces as a mesh of small triangles, each coloured by the distortion
#' at its centroid
#'
#' The mesh is cpp_globe_faces(): each face cut into its sectors (three, six
#' under IVEA), so no triangle straddles a crease of the projection. Returns the mesh with
#' `value` and `col` per triangle and the `scale` (breaks and colours) they
#' were drawn from.
#' @noRd
distortion_mesh <- function(icosa, measure, max_len = 0.025) {
  m <- cpp_globe_faces(icosa, max_len)
  idx <- matrix(m$index + 1L, ncol = 3L, byrow = TRUE)
  tri <- matrix(m$tri, ncol = 2L, byrow = TRUE)
  face <- m$item[idx[, 1]] - 1L
  cx <- (tri[idx[, 1], 1] + tri[idx[, 2], 1] + tri[idx[, 3], 1]) / 3
  cy <- (tri[idx[, 1], 2] + tri[idx[, 2], 2] + tri[idx[, 3], 2]) / 3
  value <- numeric(nrow(idx))
  for (f in unique(face)) {
    at <- which(face == f)
    s <- cpp_face_tri_tissot(icosa, f, cx[at], cy[at])
    value[at] <- if (measure == "angular") angular_deformation(s$a, s$b) else s$a * s$b
  }
  scale <- distortion_scale(measure, value)
  list(mesh = m, idx = idx, tri = tri, face = face,
       centroid = cbind(cx, cy), value = value,
       col = scale$col[findInterval(value, scale$breaks, all.inside = TRUE)],
       scale = scale)
}

#' Classes of a distortion map: angular deformation in one-degree steps from
#' 0, areal scale in steps of 0.01 around 1 on a diverging ramp
#' @noRd
distortion_scale <- function(measure, value) {
  if (measure == "angular") {
    breaks <- seq(0, max(1, ceiling(max(value))), by = 1)
    col <- grDevices::hcl.colors(length(breaks) - 1L, "YlOrRd", rev = TRUE)
  } else {
    reach <- max(0.01, ceiling(round(max(abs(value - 1)) * 100, 6)) / 100)
    breaks <- seq(1 - reach, 1 + reach, by = 0.01)
    col <- grDevices::hcl.colors(length(breaks) - 1L, "Blue-Red 3")
  }
  list(measure = measure, breaks = breaks, col = col)
}

#' Fill triangles, given as an n x 2 matrix per corner, in their colours
#' @noRd
fill_triangles <- function(p1, p2, p3, col) {
  n <- nrow(p1)
  if (n == 0L) return(invisible())
  x <- rbind(p1[, 1], p2[, 1], p3[, 1], NA)
  y <- rbind(p1[, 2], p2[, 2], p3[, 2], NA)
  graphics::polygon(as.vector(x), as.vector(y), col = col, border = col, lwd = 0.25)
}

#' The distortion map on the sphere or the solid, seen in `view`
#' @noRd
draw_distortion_surface <- function(d, view, surface, icosa) {
  key <- if (surface == "sphere") "sphere" else "solid"
  P <- matrix(d$mesh[[key]], ncol = 3L, byrow = TRUE)
  if (surface == "sphere") {
    c3 <- (P[d$idx[, 1], ] + P[d$idx[, 2], ] + P[d$idx[, 3], ]) / 3
    seen <- faces_camera(c3 / sqrt(rowSums(c3^2)), view)
  } else {
    seen <- icosa_front(icosa_solid(icosa), view)[d$face + 1L]
  }
  i <- d$idx[seen, , drop = FALSE]
  fill_triangles(project(P[i[, 1], , drop = FALSE], view),
                 project(P[i[, 2], , drop = FALSE], view),
                 project(P[i[, 3], , drop = FALSE], view), d$col[seen])
}

#' The distortion map on a net: each mesh triangle placed by every piece
#' that holds it
#' @noRd
draw_distortion_net <- function(d, net) {
  for (p in net$pieces) {
    at <- which(d$face == p$face)
    at <- at[in_triangle(d$centroid[at, , drop = FALSE], p$region)]
    if (length(at) == 0L) next
    i <- d$idx[at, , drop = FALSE]
    fill_triangles(place_points(d$tri[i[, 1], , drop = FALSE], p),
                   place_points(d$tri[i[, 2], , drop = FALSE], p),
                   place_points(d$tri[i[, 3], , drop = FALSE], p), d$col[at])
  }
}

#' A colour bar for a distortion map along the bottom of the plot region
#' @noRd
draw_distortion_legend <- function(scale, xlim, ylim) {
  n <- length(scale$col)
  w <- diff(xlim)
  x0 <- xlim[1] + 0.2 * w
  x1 <- xlim[1] + 0.8 * w
  y1 <- ylim[1] + 0.045 * diff(ylim)
  y0 <- ylim[1] + 0.025 * diff(ylim)
  xs <- seq(x0, x1, length.out = n + 1L)
  graphics::rect(xs[-(n + 1L)], y0, xs[-1L], y1, col = scale$col, border = NA)
  graphics::rect(x0, y0, x1, y1, border = "#5A636B", lwd = 0.5)
  ticks <- pretty(range(scale$breaks), n = 5)
  ticks <- ticks[ticks >= min(scale$breaks) & ticks <= max(scale$breaks)]
  tx <- x0 + (ticks - min(scale$breaks)) / diff(range(scale$breaks)) * (x1 - x0)
  graphics::text(tx, y0, labels = format(ticks), pos = 1, cex = 0.7, col = "#2E3439")
  label <- if (scale$measure == "angular") "maximum angular deformation (degrees)" else
    "areal scale"
  graphics::text(x0, y1, labels = label, pos = 3, offset = 0.3, cex = 0.7,
                 adj = 0, col = "#2E3439")
}

#' Tissot ellipses at points every `spacing` degrees, in the triangle
#' coordinates of the face each centre lies on: one closed ring per point,
#' with the face, as a matrix with columns cell (the point), face, tx, ty
#' @noRd
tissot_paths <- function(spacing, icosa, n = 49L) {
  lats <- seq(-90 + spacing, 90 - spacing / 2, by = spacing)
  pts <- do.call(rbind, lapply(lats, function(la) {
    k <- max(1L, round(360 / spacing * cos(la * pi / 180)))
    cbind(lon = seq(-180, 180, length.out = k + 1L)[-1L] - 180 / k, lat = la)
  }))
  s <- cpp_lonlat_tissot(icosa, pts[, 1], pts[, 2], rep(NA_integer_, nrow(pts)))
  # The circle's radius on the sphere: a fifth of the spacing.
  r <- spacing / 5 * pi / 180 / face_plane_edge(icosa)
  t <- seq(0, 2 * pi, length.out = n)
  rows <- lapply(seq_len(nrow(pts)), function(k) {
    c0 <- cpp_lonlat_to_face_solid(icosa, s$face[k], pts[k, 1], pts[k, 2])
    u <- r * s$a[k] * cos(t)
    v <- r * s$b[k] * sin(t)
    cbind(cell = k, face = s$face[k],
          tx = c0[1, "tx"] + cos(s$angle[k]) * u - sin(s$angle[k]) * v,
          ty = c0[1, "ty"] + sin(s$angle[k]) * u + cos(s$angle[k]) * v)
  })
  do.call(rbind, rows)
}

#' The segments of each Tissot ellipse on one face, clipped to `region`
#' (triangle coordinates): a matrix with columns x0, y0, x1, y1
#' @noRd
tissot_segments <- function(ell, face, region) {
  e <- ell[ell[, "face"] == face, , drop = FALSE]
  n <- nrow(e)
  if (n < 2L) return(matrix(numeric(0), 0L, 4L))
  s <- which(e[-1L, "cell"] == e[-n, "cell"])
  clip_segments(e[s, c("tx", "ty"), drop = FALSE],
                e[s + 1L, c("tx", "ty"), drop = FALSE], region)
}

#' Tissot ellipses on the solid, each on the face of its centre and cut at
#' the face's edges
#' @noRd
draw_tissot_solid <- function(ell, view, icosa, col, lwd) {
  front <- icosa_front(icosa_solid(icosa), view)
  for (f in which(front) - 1L) {
    seg <- tissot_segments(ell, f, FACE_TRIANGLE)
    if (nrow(seg) == 0L) next
    a <- project(cpp_face_tri_to_solid(icosa, f, seg[, 1], seg[, 2]), view)
    b <- project(cpp_face_tri_to_solid(icosa, f, seg[, 3], seg[, 4]), view)
    graphics::segments(a[, 1], a[, 2], b[, 1], b[, 2], col = col, lwd = lwd)
  }
}

#' Tissot ellipses on a net, cut to each piece of their face
#' @noRd
draw_tissot_net <- function(ell, net, col, lwd) {
  for (p in net$pieces) {
    seg <- tissot_segments(ell, p$face, p$region)
    if (nrow(seg) == 0L) next
    draw_net_segments(cbind(place_points(seg[, 1:2, drop = FALSE], p),
                            place_points(seg[, 3:4, drop = FALSE], p)), col, lwd)
  }
}
