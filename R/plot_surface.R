# plot_surface.R
# A grid drawn on the sphere and on the icosahedron it is built on

#' Globe center presets
#'
#' Named list of lon/lat coordinates for common globe views, accepted by the
#' \code{center} argument of the grid \code{\link[=plot,HexGridInfo,missing-method]{plot}}
#' method.
#'
#' @format Named list with elements:
#' \describe{
#'   \item{europe}{c(10, 50) - Western/Central Europe}
#'   \item{north_america}{c(-100, 45) - USA and Canada}
#'   \item{south_america}{c(-60, -15) - Full continent}
#'   \item{africa}{c(20, 5) - Central Africa}
#'   \item{asia}{c(100, 35) - China, SE Asia, Japan}
#'   \item{oceania}{c(135, -25) - Australia, NZ, Indonesia}
#'   \item{middle_east}{c(45, 25) - Arabian Peninsula, Iran, Turkey}
#'   \item{south_asia}{c(80, 20) - India, Pakistan, Bangladesh}
#'   \item{pacific}{c(-160, -10) - Polynesia, Pacific islands}
#'   \item{caribbean}{c(-70, 18) - Caribbean islands}
#'   \item{arctic}{c(0, 90) - North pole view}
#'   \item{antarctic}{c(0, -90) - South pole view}
#' }
#'
#' @export
#' @examples
#' globe_centers$europe
#' globe_centers$oceania
globe_centers <- list(
  europe        = c(lon = 10,   lat = 50),
  north_america = c(lon = -100, lat = 45),
  south_america = c(lon = -60,  lat = -15),
  africa        = c(lon = 20,   lat = 5),
  asia          = c(lon = 100,  lat = 35),
  oceania       = c(lon = 135,  lat = -25),
  middle_east   = c(lon = 45,   lat = 25),
  south_asia    = c(lon = 80,   lat = 20),
  pacific       = c(lon = -160, lat = -10),
  caribbean     = c(lon = -70,  lat = 18),
  arctic        = c(lon = 0,    lat = 90),
  antarctic     = c(lon = 0,    lat = -90)
)

#' Plot a grid on the sphere or on its icosahedron
#'
#' Draws the cells of a grid in 3D, seen from above a chosen point: on the
#' sphere, or on the flat faces of the icosahedron the grid is built on. Both
#' surfaces take their cell boundaries from the same points, so a cell on the
#' icosahedron is the cell on the sphere folded flat.
#'
#' @param x A HexGridInfo object from \code{\link{hex_grid}}
#' @param y Ignored
#' @param surface \code{"sphere"} or \code{"icosahedron"}. The icosahedron
#'   needs an ISEA grid; an H3 grid is drawn on the sphere.
#' @param center Point the view looks down on: a preset name from
#'   \code{\link{globe_centers}} or \code{c(lon, lat)}.
#' @param cells Cell IDs to draw. \code{NULL} draws every cell of the grid.
#' @param land \code{TRUE} for the built-in \code{\link{hexify_world}},
#'   \code{FALSE} for none, or an sf/sfc object of polygons.
#' @param ocean_fill Fill colour of the surface.
#' @param land_fill Fill colour of land; \code{NA} leaves land unfilled.
#' @param land_border Colour of country outlines; \code{NA} draws none.
#' @param land_lwd Line width of country outlines.
#' @param grid_border Colour of cell boundaries.
#' @param grid_lwd Line width of cell boundaries.
#' @param face_edges Draw the edges of the icosahedron's faces. \code{NULL}
#'   draws them for an ISEA grid; H3 is built on an icosahedron of its own,
#'   so its grid takes none.
#' @param edge_col Colour of face edges.
#' @param edge_lwd Line width of face edges.
#' @param step Spacing of the points along a cell boundary, as a fraction of
#'   a face edge.
#' @param main Plot title.
#' @param ... Ignored
#'
#' @return The grid, invisibly
#'
#' @seealso \code{\link{grid_global}} for the cells as sf polygons
#'
#' @export
#' @examples
#' grid <- hex_grid(resolution = 3, aperture = 3)
#' plot(grid)
#' plot(grid, surface = "icosahedron")
#' plot(grid, surface = "icosahedron", land = FALSE, center = "pacific")
setMethod("plot", signature(x = "HexGridInfo", y = "missing"),
  function(x, y,
           surface = c("sphere", "icosahedron"),
           center = c(lon = 15, lat = 32),
           cells = NULL,
           land = TRUE,
           ocean_fill = "#F4F6F8",
           land_fill = "#C9D0D6",
           land_border = "#8C969F",
           land_lwd = 0.4,
           grid_border = "#0072B2",
           grid_lwd = 0.8,
           face_edges = NULL,
           edge_col = "black",
           edge_lwd = 1.1,
           step = 0.01,
           main = NULL,
           ...) {
    surface <- match.arg(surface)
    g <- extract_grid(x)
    if (surface == "icosahedron" && is_h3_grid(g)) {
      stop("surface = \"icosahedron\" needs an ISEA grid; H3 cells are drawn ",
           "on the sphere", call. = FALSE)
    }

    if (is.null(face_edges)) face_edges <- !is_h3_grid(g)
    if (face_edges && is_h3_grid(g)) {
      stop("face_edges = TRUE draws the ISEA icosahedron, which an H3 grid is ",
           "not built on", call. = FALSE)
    }

    view <- surface_view(resolve_center(center))
    paths <- grid_surface_paths(g, cells, step)
    land <- surface_land(land)
    style <- list(ocean_fill = ocean_fill, land_fill = land_fill,
                  land_border = land_border, land_lwd = land_lwd,
                  grid_border = grid_border, grid_lwd = grid_lwd,
                  face_edges = face_edges, edge_col = edge_col,
                  edge_lwd = edge_lwd)

    graphics::plot.new()
    old <- graphics::par(mar = c(0, 0, if (is.null(main)) 0 else 2, 0))
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot.window(c(-1.02, 1.02), c(-1.02, 1.02), asp = 1)
    if (surface == "sphere") {
      draw_sphere(paths, land, view, style)
    } else {
      draw_icosahedron(paths, land, view, style)
    }
    if (!is.null(main)) graphics::title(main = main)
    invisible(x)
  }
)

# =============================================================================
# VIEW
# =============================================================================

#' Resolve center preset or coordinates
#' @noRd
resolve_center <- function(center) {
  if (is.character(center)) {
    if (!center %in% names(globe_centers)) {
      valid <- paste(names(globe_centers), collapse = ", ")
      stop(sprintf("Unknown center preset '%s'. Valid options: %s", center, valid))
    }
    return(globe_centers[[center]])
  }

  if (is.numeric(center) && length(center) == 2) {
    if (is.null(names(center))) {
      return(c(lon = center[[1]], lat = center[[2]]))
    }
    if (!setequal(names(center), c("lon", "lat"))) {
      stop(sprintf(
        "center must be named c(lon = ..., lat = ...); got names: %s",
        paste(names(center), collapse = ", ")
      ))
    }
    return(center[c("lon", "lat")])
  }

  stop("center must be a preset name or numeric c(lon, lat)")
}

#' Unit vectors of lon/lat points
#' @noRd
unit_vec <- function(lon, lat) {
  lon <- lon * pi / 180
  lat <- lat * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

#' The orthographic view down onto a point of the sphere
#'
#' Columns of `basis` are screen right, screen up and the direction towards
#' the viewer, so `P %*% basis` gives screen x, y and depth. `light` is the
#' direction light comes from, fixed above the viewer's left shoulder.
#' @noRd
surface_view <- function(center) {
  lon <- center[["lon"]] * pi / 180
  e3 <- drop(unit_vec(center[["lon"]], center[["lat"]]))
  e1 <- c(-sin(lon), cos(lon), 0)
  if (sum(e1^2) < 1e-12) e1 <- c(0, 1, 0)
  e2 <- c(e3[2] * e1[3] - e3[3] * e1[2],
          e3[3] * e1[1] - e3[1] * e1[3],
          e3[1] * e1[2] - e3[2] * e1[1])
  basis <- cbind(e1, e2, e3)
  light <- c(-0.45, 0.55, 0.70)
  list(basis = basis, light = drop(basis %*% (light / sqrt(sum(light^2)))))
}

#' Draw the visible pieces of polylines given in 3D
#'
#' `P` holds points of the scene, `brk` marks rows that start a new line and
#' `keep` the rows whose segment to the next row is drawn.
#' @noRd
draw_segments <- function(P, view, brk, keep, col, lwd) {
  n <- nrow(P)
  if (n < 2L) return(invisible())
  xy <- P %*% view$basis
  s <- which(!brk[-1L] & keep[-n])
  graphics::segments(xy[s, 1], xy[s, 2], xy[s + 1L, 1], xy[s + 1L, 2],
                     col = col, lwd = lwd, lend = "round")
}

#' Fill rings given in 3D as one even-odd path
#' @noRd
fill_rings <- function(rings, view, col) {
  if (length(rings) == 0L || is.na(col)) return(invisible())
  xy <- do.call(rbind, lapply(rings, function(r) rbind((r %*% view$basis)[, 1:2, drop = FALSE], NA)))
  xy <- xy[-nrow(xy), , drop = FALSE]
  graphics::polypath(xy[, 1], xy[, 2], col = col, border = NA, rule = "evenodd")
}

#' A colour darkened by a shading factor
#' @noRd
shade_col <- function(col, s) {
  if (is.na(col)) return(col)
  v <- grDevices::col2rgb(col)[, 1] / 255 * s
  grDevices::rgb(v[1], v[2], v[3])
}

# =============================================================================
# GEOMETRY
# =============================================================================

#' The icosahedron of the ISEA grids
#'
#' Vertices on the unit sphere, the vertex indices of each face, each face's
#' outward normal, and its 30 edges as vertex pairs with the two faces that
#' share them.
#' @noRd
icosa_solid <- function() {
  s <- cpp_icosa_solid()
  V <- s$vertices
  Fv <- s$faces
  N <- t(apply(Fv, 1, function(f) {
    c <- colMeans(V[f, ])
    c / sqrt(sum(c^2))
  }))
  pairs <- do.call(rbind, lapply(seq_len(nrow(Fv)), function(f) {
    v <- Fv[f, ]
    cbind(pmin(v, v[c(2, 3, 1)]), pmax(v, v[c(2, 3, 1)]), f)
  }))
  key <- paste(pairs[, 1], pairs[, 2])
  edges <- t(vapply(split(seq_len(nrow(pairs)), key), function(i) {
    c(pairs[i[1], 1:2], pairs[i, 3])
  }, numeric(4)))
  colnames(edges) <- c("v1", "v2", "f1", "f2")
  list(vertices = V, faces = Fv, normals = N, edges = edges)
}

#' Cell boundaries of a grid on the sphere and on the icosahedron
#'
#' One closed path per cell, as a matrix with columns cell, face, solid_x/y/z
#' and sphere_x/y/z. An ISEA grid reads them from its faces; an H3 cell edge is
#' a great-circle arc between corners, so H3 paths carry sphere positions only.
#' @noRd
grid_surface_paths <- function(g, cells, step) {
  if (is_h3_grid(g)) {
    if (is.null(cells)) cells <- h3_all_cells(g@resolution)
    return(h3_sphere_paths(as.character(cells), step))
  }
  if (is.null(cells)) cells <- seq_len(grid_n_cells(g))
  mixed <- is_mixed_aperture(g@aperture)
  cpp_cell_surface_paths(
    as.numeric(cells), g@resolution,
    if (mixed) 0L else aperture_to_int(g@aperture),
    if (mixed) grid_ap_seq(g) else integer(0),
    step
  )
}

#' H3 cell boundaries as great-circle arcs on the sphere
#' @noRd
h3_sphere_paths <- function(cells, step) {
  arc <- step * atan(2)
  boundaries <- cpp_h3_cellToBoundary(cells)
  rows <- lapply(seq_along(cells), function(k) {
    b <- lonlat_ring_coords(boundaries[[k]])
    V <- unit_vec(b[, 1], b[, 2])
    if (any(V[1, ] != V[nrow(V), ])) V <- rbind(V, V[1, ])
    pts <- do.call(rbind, lapply(seq_len(nrow(V) - 1L), function(i) {
      slerp(V[i, ], V[i + 1L, ], arc)
    }))
    pts <- rbind(pts, V[1, ])
    cbind(cell = k, face = NA_real_, solid_x = NA_real_, solid_y = NA_real_,
          solid_z = NA_real_, sphere_x = pts[, 1], sphere_y = pts[, 2],
          sphere_z = pts[, 3])
  })
  do.call(rbind, rows)
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

#' The land a plot draws
#' @noRd
surface_land <- function(land) {
  if (isFALSE(land) || is.null(land)) return(NULL)
  if (isTRUE(land)) land <- hexify::hexify_world
  if (inherits(land, "sf")) land <- sf::st_geometry(land)
  if (!inherits(land, "sfc")) {
    stop("land must be TRUE, FALSE, or an sf/sfc object", call. = FALSE)
  }
  land <- sf::st_transform(land, 4326)
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  land <- suppressMessages(sf::st_make_valid(land))
  land[!sf::st_is_empty(land)]
}

#' Polygons of land within a spherical cap, as lists of lon/lat rings
#' @noRd
land_in_cap <- function(land, lon, lat, radius_deg) {
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  cap <- sf::st_buffer(sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326),
                       radius_deg * pi / 180 * 6371008.8)
  part <- suppressMessages(sf::st_intersection(land, cap))
  part <- part[!sf::st_is_empty(part)]
  if (length(part) == 0L) return(list())
  part <- polygons_only(part)
  unlist(lapply(part, function(p) {
    if (inherits(p, "MULTIPOLYGON")) unclass(p) else list(unclass(p))
  }), recursive = FALSE)
}

#' The polygons of a geometry set, dropping the points and lines an
#' intersection can leave
#' @noRd
polygons_only <- function(x) {
  if (any(sf::st_is(x, "GEOMETRYCOLLECTION"))) {
    x <- sf::st_collection_extract(x, "POLYGON")
  }
  x[sf::st_is(x, c("POLYGON", "MULTIPOLYGON"))]
}

#' Country outlines as one set of 3D polylines on the unit sphere
#' @noRd
land_outline_points <- function(land) {
  rings <- unlist(lapply(land, function(p) {
    if (inherits(p, "MULTIPOLYGON")) unlist(unclass(p), recursive = FALSE) else unclass(p)
  }), recursive = FALSE)
  P <- do.call(rbind, lapply(rings, function(r) unit_vec(r[, 1], r[, 2])))
  brk <- unlist(lapply(rings, function(r) c(TRUE, rep(FALSE, nrow(r) - 1L))))
  list(P = P, brk = brk)
}

# =============================================================================
# SURFACES
# =============================================================================

#' Draw the grid on the sphere
#' @noRd
draw_sphere <- function(paths, land, view, style) {
  ring <- seq(0, 2 * pi, length.out = 361)
  graphics::polygon(cos(ring), sin(ring), col = style$ocean_fill, border = NA)
  e3 <- view$basis[, 3]
  center_ll <- c(atan2(e3[2], e3[1]), asin(e3[3])) * 180 / pi

  if (!is.null(land)) {
    for (p in land_in_cap(land, center_ll[1], center_ll[2], 89.5)) {
      fill_rings(lapply(p, function(r) unit_vec(r[, 1], r[, 2])), view, style$land_fill)
    }
    if (!is.na(style$land_border)) {
      o <- land_outline_points(land)
      draw_segments(o$P, view, o$brk, drop(o$P %*% e3) > 0 & c(drop(o$P %*% e3)[-1], 0) > 0,
                    style$land_border, style$land_lwd)
    }
  }

  S <- paths[, c("sphere_x", "sphere_y", "sphere_z"), drop = FALSE]
  depth <- drop(S %*% e3)
  front <- depth > 0
  draw_segments(S, view, c(TRUE, diff(paths[, "cell"]) != 0),
                front & c(front[-1], FALSE), style$grid_border, style$grid_lwd)

  if (style$face_edges) {
    solid <- icosa_solid()
    V <- solid$vertices
    arcs <- lapply(seq_len(nrow(solid$edges)), function(j) {
      a <- V[solid$edges[j, "v1"], ]
      b <- V[solid$edges[j, "v2"], ]
      rbind(slerp(a, b, 0.5 * pi / 180), b)
    })
    P <- do.call(rbind, arcs)
    brk <- unlist(lapply(arcs, function(a) c(TRUE, rep(FALSE, nrow(a) - 1L))))
    vis <- drop(P %*% e3) > 0
    draw_segments(P, view, brk, vis & c(vis[-1], FALSE), style$edge_col, style$edge_lwd)
  }
  graphics::lines(cos(ring), sin(ring), col = "#9AA3AB", lwd = 0.8)
}

#' Draw the grid on the icosahedron
#' @noRd
draw_icosahedron <- function(paths, land, view, style) {
  solid <- icosa_solid()
  V <- solid$vertices
  e3 <- view$basis[, 3]
  front <- drop(solid$normals %*% e3) > 0
  shade <- 0.80 + 0.20 * pmax(0, drop(solid$normals %*% view$light))

  for (f in which(front)) {
    tri <- V[solid$faces[f, ], ]
    graphics::polygon((tri %*% view$basis)[, 1:2],
                      col = shade_col(style$ocean_fill, shade[f]), border = NA)
    if (!is.null(land)) draw_face_land(f - 1L, tri, land, view, style, shade[f])
  }

  face <- paths[, "face"] + 1L
  Q <- paths[, c("solid_x", "solid_y", "solid_z"), drop = FALSE]
  draw_segments(Q, view, c(TRUE, diff(paths[, "cell"]) != 0), front[face],
                style$grid_border, style$grid_lwd)

  if (style$face_edges) {
    E <- solid$edges
    vis <- front[E[, "f1"]] | front[E[, "f2"]]
    for (j in which(vis)) {
      graphics::lines((V[E[j, c("v1", "v2")], ] %*% view$basis)[, 1:2],
                      col = style$edge_col, lwd = style$edge_lwd)
    }
  }
}

#' Land and country outlines on one flat face
#'
#' Land within a cap around the face's centre is read in the face's triangle
#' coordinates, cut to the triangle there, and placed on the face. The cap is
#' wider than the face's circumradius, so its rim never reaches the triangle.
#' @noRd
draw_face_land <- function(face, tri, land, view, style, shade) {
  c3 <- colMeans(tri)
  c_ll <- c(atan2(c3[2], c3[1]), asin(c3[3] / sqrt(sum(c3^2)))) * 180 / pi
  polys <- land_in_cap(land, c_ll[1], c_ll[2], 40)
  if (length(polys) == 0L) return(invisible())

  to_tri <- function(r) cpp_lonlat_to_face_solid(face, r[, 1], r[, 2])[, c("tx", "ty"), drop = FALSE]
  tri_ll <- cbind(atan2(tri[, 2], tri[, 1]), asin(tri[, 3])) * 180 / pi
  tri_t <- to_tri(tri_ll)
  triangle <- sf::st_polygon(list(rbind(tri_t, tri_t[1, ])))

  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  flat <- sf::st_sfc(lapply(polys, function(p) {
    sf::st_polygon(lapply(p, function(r) {
      t_r <- to_tri(r)
      rbind(t_r[-nrow(t_r), , drop = FALSE], t_r[1, ])
    }))
  }))
  flat <- suppressWarnings(sf::st_make_valid(flat))
  to_solid <- function(m) cpp_face_tri_to_solid(face, m[, 1], m[, 2])

  if (!is.na(style$land_fill)) {
    filled <- suppressWarnings(sf::st_intersection(flat, sf::st_sfc(triangle)))
    filled <- filled[!sf::st_is_empty(filled)]
    if (length(filled) > 0L) {
      filled <- polygons_only(filled)
      rings <- unlist(lapply(filled, function(p) {
        if (inherits(p, "MULTIPOLYGON")) unlist(unclass(p), recursive = FALSE) else unclass(p)
      }), recursive = FALSE)
      fill_rings(lapply(rings, to_solid), view, shade_col(style$land_fill, shade))
    }
  }

  if (!is.na(style$land_border)) {
    lines_t <- suppressWarnings(sf::st_intersection(sf::st_boundary(flat), sf::st_sfc(triangle)))
    lines_t <- lines_t[!sf::st_is_empty(lines_t)]
    if (length(lines_t) > 0L) {
      lines_t <- lines_t[sf::st_dimension(lines_t) == 1L]
      if (any(sf::st_is(lines_t, "GEOMETRYCOLLECTION"))) {
        lines_t <- sf::st_collection_extract(lines_t, "LINESTRING")
      }
      parts <- unlist(lapply(lines_t, function(l) {
        if (inherits(l, "MULTILINESTRING")) unclass(l) else list(unclass(l))
      }), recursive = FALSE)
      P <- do.call(rbind, lapply(parts, to_solid))
      brk <- unlist(lapply(parts, function(l) c(TRUE, rep(FALSE, nrow(l) - 1L))))
      draw_segments(P, view, brk, rep(TRUE, nrow(P)), style$land_border, style$land_lwd)
    }
  }
}
