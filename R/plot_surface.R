# plot_surface.R
# A grid drawn on the sphere and on the solid it is built on

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

#' Plot a grid on the sphere, on its solid, or on the unfolded net
#'
#' Draws the cells of a grid in 3D, seen from above a chosen point: on the
#' sphere, or on the flat faces of the solid the grid is built on. All
#' surfaces take their cell boundaries from the same points, so a cell on the
#' solid is the cell on the sphere folded flat.
#'
#' \code{surface = "net"} lays the 20 faces out flat in the PLANE layout of
#' DGGRID (\code{\link{hexify_cell_to_plane}} gives cell centres in the same
#' coordinates). A cell on a face edge that the net cuts is drawn in two
#' parts, one on each face. The net is drawn without a camera, so
#' \code{center}, \code{projection}, \code{distance}, \code{tilt},
#' \code{rotation} and \code{fov} apply to the other two surfaces.
#'
#' The view is an orthographic projection by default: parallel lines of sight,
#' so the whole near half of the sphere shows. \code{projection =
#' "perspective"} places a camera \code{distance} sphere radii from the
#' sphere's centre, above \code{center}; it sees a smaller cap, and nearer
#' cells appear larger. \code{tilt} swings that camera about the point at
#' \code{center}, which it keeps looking at, so the surface is seen at a slant
#' with the horizon towards the top; by default the plot frames the whole
#' visible sphere. \code{fov} sets how wide
#' the camera sees, and \code{rotation} turns the view about the line of sight
#' in either projection.
#'
#' @param x A HexGridInfo object from \code{\link{hex_grid}}
#' @param y Ignored
#' @param surface \code{"sphere"}, \code{"solid"} or \code{"net"}.
#'   \code{"solid"} draws the flat faces of the grid's polyhedron (the
#'   icosahedron, or the octahedron of a grid built on it). The solid and the
#'   net need an ISEA grid; an H3 grid is drawn on the
#'   sphere.
#' @param center Point the view looks down on: a preset name from
#'   \code{\link{globe_centers}} or \code{c(lon, lat)}.
#' @param projection \code{"orthographic"} or \code{"perspective"}.
#' @param distance Distance of the untilted perspective camera from the
#'   sphere's centre, in sphere radii; greater than 1. \code{NULL} uses 3. A
#'   geostationary satellite sits at about 6.6 Earth radii.
#' @param tilt Degrees the perspective camera swings about the point at
#'   \code{center}, between -90 and 90, keeping its range of
#'   \code{distance - 1} radii to that point. Positive values move the camera
#'   towards the bottom of the view, so it looks towards the top.
#' @param rotation Degrees the view turns about the line of sight; positive
#'   values turn north clockwise.
#' @param fov Field of view of the perspective camera in degrees, across the
#'   square the plot shows, centred on \code{center}; below 170. \code{NULL}
#'   frames the whole visible sphere instead.
#' @param cells Cell IDs to draw. \code{NULL} draws every cell of the grid.
#' @param land \code{TRUE} for the built-in \code{\link{hexify_world}},
#'   \code{FALSE} for none, or an sf/sfc object of polygons.
#' @param ocean_fill Fill colour of the surface.
#' @param land_fill Fill colour of land; \code{NA} leaves land unfilled.
#' @param land_border Colour of country outlines; \code{NA} draws none.
#' @param land_lwd Line width of country outlines.
#' @param grid_border Colour of cell boundaries.
#' @param grid_lwd Line width of cell boundaries.
#' @param face_edges Draw the edges of the solid's faces. \code{NULL}
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
#' plot(grid, surface = "solid")
#' plot(grid, surface = "solid", land = FALSE, center = "pacific")
#' plot(grid, projection = "perspective", distance = 2, center = "europe")
#' plot(grid, surface = "solid", projection = "perspective",
#'      distance = 2.5, tilt = 25, rotation = 30)
#' plot(grid, surface = "net")
setMethod("plot", signature(x = "HexGridInfo", y = "missing"),
  function(x, y,
           surface = c("sphere", "solid", "net"),
           center = c(lon = 15, lat = 32),
           projection = c("orthographic", "perspective"),
           distance = NULL,
           tilt = 0,
           rotation = 0,
           fov = NULL,
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
    face_edges <- resolve_surface(surface, face_edges, g)

    projection <- match.arg(projection)
    if (surface == "net" && (!missing(center) || projection != "orthographic" ||
                             !is.null(distance) || !is.null(fov) ||
                             tilt != 0 || rotation != 0)) {
      stop("the net is drawn flat; center, projection, distance, tilt, ",
           "rotation and fov apply to the sphere and the solid",
           call. = FALSE)
    }
    camera <- resolve_camera(projection, distance, tilt, rotation, fov)
    view <- surface_view(resolve_center(center), camera$distance, tilt, rotation)
    paths <- grid_surface_paths(g, cells, step)
    land <- surface_land(land)
    icosa <- icosa_arg(g)
    style <- list(ocean_fill = ocean_fill, land_fill = land_fill,
                  land_border = land_border, land_lwd = land_lwd,
                  grid_border = grid_border, grid_lwd = grid_lwd,
                  face_edges = face_edges, edge_col = edge_col,
                  edge_lwd = edge_lwd)

    old <- graphics::par(mar = c(0, 0, if (is.null(main)) 0 else 2, 0))
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot.new()
    if (surface == "net") {
      tris <- net_triangles(icosa)
      xy <- do.call(rbind, lapply(tris, `[[`, "plane"))
      pad <- 0.02 * diff(range(xy[, 1]))
      xlim <- range(xy[, 1]) + c(-1, 1) * pad
      ylim <- range(xy[, 2]) + c(-1, 1) * pad
    } else {
      frame <- view_frame(view, camera$fov, surface_outline(surface, view, icosa))
      xlim <- frame[1] + c(-1, 1) * frame[3]
      ylim <- frame[2] + c(-1, 1) * frame[3]
    }
    graphics::plot.window(xlim, ylim, asp = 1, xaxs = "i", yaxs = "i")
    graphics::clip(xlim[1], xlim[2], ylim[1], ylim[2])
    switch(surface,
      sphere = draw_sphere(paths, land, view, style, icosa),
      solid = draw_solid(paths, land, view, style, icosa),
      net = draw_net(paths, land, tris, style, icosa)
    )
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

#' Check the surface a grid is drawn on, and whether its face edges are drawn
#'
#' Returns `face_edges`, \code{NULL} becoming TRUE for an ISEA grid and FALSE for H3.
#' @noRd
resolve_surface <- function(surface, face_edges, g) {
  if (surface != "sphere" && is_h3_grid(g)) {
    stop(sprintf("surface = \"%s\" needs an ISEA grid; H3 cells are drawn ",
                 surface), "on the sphere", call. = FALSE)
  }
  if (is.null(face_edges)) face_edges <- !is_h3_grid(g)
  if (face_edges && is_h3_grid(g)) {
    stop("face_edges = TRUE draws the ISEA solid, which an H3 grid is ",
         "not built on", call. = FALSE)
  }
  face_edges
}

#' Check the projection arguments of the plot method
#'
#' Returns the camera distance (\code{Inf} for the orthographic view) and the
#' field of view in degrees (\code{NA} for the orthographic view).
#' @noRd
resolve_camera <- function(projection, distance, tilt, rotation, fov = NULL) {
  for (a in list(list("tilt", tilt), list("rotation", rotation))) {
    if (!is.numeric(a[[2]]) || length(a[[2]]) != 1L || !is.finite(a[[2]])) {
      stop(a[[1]], " must be a single number of degrees", call. = FALSE)
    }
  }
  if (projection == "orthographic") {
    if (!is.null(distance)) {
      stop("distance places a perspective camera; set projection = \"perspective\"",
           call. = FALSE)
    }
    if (tilt != 0) {
      stop("tilt pitches a perspective camera; set projection = \"perspective\"",
           call. = FALSE)
    }
    if (!is.null(fov)) {
      stop("fov frames a perspective camera; set projection = \"perspective\"",
           call. = FALSE)
    }
    return(list(distance = Inf, fov = NA_real_))
  }
  if (is.null(distance)) distance <- 3
  if (!is.numeric(distance) || length(distance) != 1L || !is.finite(distance) ||
      distance <= 1) {
    stop("distance must be a single number greater than 1 (sphere radii from ",
         "the centre)", call. = FALSE)
  }
  if (abs(tilt) >= 90) {
    stop("tilt must lie strictly between -90 and 90 degrees", call. = FALSE)
  }
  if (is.null(fov)) return(list(distance = distance, fov = NA_real_))
  if (!is.numeric(fov) || length(fov) != 1L || !is.finite(fov) ||
      fov <= 0 || fov >= 170) {
    stop("fov must be a single number of degrees between 0 and 170",
         call. = FALSE)
  }
  list(distance = distance, fov = fov)
}

#' The square of screen a plot shows, as its centre and half its width
#'
#' Orthographic: the unit disc of the sphere. Perspective with a field of
#' view: the square that field of view spans around the middle of the view.
#' Perspective with `fov = NA`: the smallest square that holds `outline`, the
#' points in the scene that bound what is drawn. For the sphere that is its
#' rim, and an untilted camera frames the unit disc again. The outline is
#' first held to the square 120 degrees wide around the middle of the view: a
#' camera close and tilted sees the nearby surface spread far wider, and
#' outline points nearing the camera's plane run off towards infinity before
#' they pass behind it, so held there they keep the frame at that edge.
#' @noRd
view_frame <- function(view, fov, outline = horizon_ring(view, 721L)) {
  if (is.null(view$eye)) return(c(0, 0, 1.02))
  if (!is.na(fov)) return(c(0, 0, 1.02 * view$scale * tan(fov / 2 * pi / 180)))
  widest <- view$scale * tan(60 * pi / 180)
  xy <- project(outline, view)
  xy <- pmin(pmax(xy[stats::complete.cases(xy), , drop = FALSE], -widest), widest)
  lo <- apply(xy, 2, min)
  hi <- apply(xy, 2, max)
  c((lo + hi) / 2, 1.02 * max(hi - lo) / 2)
}

#' Points that bound the drawing of a surface
#'
#' The sphere's rim, or the vertices of the solid's faces that face the
#' camera: the solid's outline runs along edges between a face that
#' faces the camera and one that does not, so its corners are among these.
#' @noRd
surface_outline <- function(surface, view, icosa) {
  if (surface == "sphere") return(horizon_ring(view, 721L))
  solid <- icosa_solid(icosa)
  front <- icosa_front(solid, view)
  solid$vertices[unique(as.vector(solid$faces[front, ])), , drop = FALSE]
}

#' Icosahedron faces the camera sees
#'
#' A face is seen when the camera lies outside its plane, which sits at
#' distance h from the centre along the face normal.
#' @noRd
icosa_front <- function(solid, view) {
  h <- rowSums(solid$normals * solid$vertices[solid$faces[, 1], ])
  drop(solid$normals %*% view$dir) > h * view$horizon
}

#' The view of the sphere from above a point
#'
#' The camera looks at the point `center` on the sphere, turned by `rotation`
#' degrees about the line of sight (north turns clockwise). Untilted it sits
#' straight above that point, `distance` sphere radii from the sphere's centre
#' (`Inf` for the orthographic view). `tilt` swings it by that many degrees
#' about the point, towards the bottom of the view, at the same range
#' `distance - 1` from the point, so the point stays in the middle of the view.
#'
#' Columns of `cam` are the camera's right, up and backward axes in the scene.
#' `dir` points from the sphere's centre to the camera, and a point P of the
#' unit sphere faces the camera when `P . dir > horizon`, with `horizon` the
#' cosine of the visible cap's angular radius: 0 seen from infinitely far,
#' one over the camera's distance from the centre otherwise. `u` and `v`
#' complete `dir` to an orthonormal basis. `light` is the direction light
#' comes from, fixed above the camera's left shoulder.
#' @noRd
surface_view <- function(center, distance = Inf, tilt = 0, rotation = 0) {
  lon <- center[["lon"]] * pi / 180
  e3 <- drop(unit_vec(center[["lon"]], center[["lat"]]))
  e1 <- c(-sin(lon), cos(lon), 0)
  e2 <- cross3(e3, e1)
  r <- rotation * pi / 180
  u <- cos(r) * e1 + sin(r) * e2
  v <- -sin(r) * e1 + cos(r) * e2
  t <- tilt * pi / 180
  back <- cos(t) * e3 - sin(t) * v
  cam <- unname(cbind(u, cos(t) * v + sin(t) * e3, back))
  light <- c(-0.45, 0.55, 0.70)
  light <- drop(cam %*% (light / sqrt(sum(light^2))))
  if (!is.finite(distance)) {
    return(list(cam = cam, eye = NULL, scale = 1, near = NA_real_,
                dir = e3, horizon = 0, u = u, v = v, light = light))
  }
  eye <- e3 + (distance - 1) * back
  reach <- sqrt(sum(eye^2))
  dir <- eye / reach
  list(cam = cam, eye = eye, scale = sqrt(distance^2 - 1),
       near = 0.01 * (reach - 1), dir = dir, horizon = 1 / reach, u = u,
       v = cross3(dir, u), light = light)
}

#' Points in camera coordinates: right, up, and depth in front of the camera
#' @noRd
camera_coords <- function(P, view) {
  q <- sweep(P, 2, view$eye) %*% view$cam
  q[, 3] <- -q[, 3]
  q
}

#' Screen coordinates of points in the scene
#'
#' Orthographic: the points' coordinates along the camera's right and up
#' axes. Perspective: the points seen through a pinhole at the camera, scaled
#' so that the untilted sphere fills the unit disc. Points closer to the
#' camera's plane than `view$near` are NA. No point of the sphere lies within
#' `distance - 1` of the camera, so a point that close to the plane is at least
#' 100 times further sideways than ahead and falls outside any field of view.
#' @noRd
project <- function(P, view) {
  if (is.null(view$eye)) return(P %*% view$cam[, 1:2])
  q <- camera_coords(P, view)
  z <- q[, 3]
  z[z < view$near] <- NA
  view$scale * q[, 1:2, drop = FALSE] / z
}

#' Screen coordinates of a closed ring, cut to the part in front of the camera
#'
#' The ring is clipped against the plane `view$near` ahead of the camera
#' (Sutherland-Hodgman against one plane) before it is projected, so a filled
#' shape that reaches behind the camera keeps the part the camera sees.
#' Returns NULL when nothing is left.
#' @noRd
project_ring <- function(R, view) {
  if (is.null(view$eye)) return(R %*% view$cam[, 1:2])
  q <- camera_coords(R, view)
  inside <- q[, 3] >= view$near
  if (!any(inside)) return(NULL)
  if (!all(inside)) {
    n <- nrow(q)
    nxt <- c(seq_len(n)[-1L], 1L)
    cross <- which(inside != inside[nxt])
    a <- q[cross, , drop = FALSE]
    b <- q[nxt[cross], , drop = FALSE]
    s <- (view$near - a[, 3]) / (b[, 3] - a[, 3])
    q <- rbind(q[inside, , drop = FALSE], a + s * (b - a))
    q <- q[order(c(which(inside), cross + 0.5)), , drop = FALSE]
  }
  view$scale * q[, 1:2, drop = FALSE] / q[, 3]
}

#' Points of the unit sphere on the side that faces the camera
#' @noRd
faces_camera <- function(P, view) {
  drop(P %*% view$dir) > view$horizon
}

#' The rim of the sphere as the camera sees it
#' @noRd
horizon_ring <- function(view, n = 361L) {
  a <- seq(0, 2 * pi, length.out = n)
  outer(rep(view$horizon, n), view$dir) +
    sqrt(1 - view$horizon^2) * (outer(cos(a), view$u) + outer(sin(a), view$v))
}

#' Draw the visible pieces of polylines given in 3D
#'
#' `P` holds points of the scene, `brk` marks rows that start a new line and
#' `keep` the rows whose segment to the next row is drawn.
#' @noRd
draw_segments <- function(P, view, brk, keep, col, lwd) {
  n <- nrow(P)
  if (n < 2L) return(invisible())
  xy <- project(P, view)
  s <- which(!brk[-1L] & keep[-n])
  graphics::segments(xy[s, 1], xy[s, 2], xy[s + 1L, 1], xy[s + 1L, 2],
                     col = col, lwd = lwd, lend = "round")
}

#' Fill polygons given in 3D as one path
#'
#' Each polygon is a list of rings, its exterior first and then its holes.
#' On screen every exterior is turned anticlockwise and every hole
#' clockwise, and the path is filled by the nonzero winding rule, so holes
#' stay open while polygons that overlap, such as neighbouring countries whose
#' simplified borders cross, stay filled where they overlap.
#' @noRd
fill_polygons <- function(polys, view, col) {
  if (length(polys) == 0L || is.na(col)) return(invisible())
  xy <- do.call(rbind, lapply(polys, function(rings) {
    do.call(rbind, lapply(seq_along(rings), function(k) {
      p <- project_ring(rings[[k]], view)
      if (is.null(p) || nrow(p) < 3L) return(NULL)
      rbind(orient_ring(p, anticlockwise = k == 1L), NA)
    }))
  }))
  if (is.null(xy)) return(invisible())
  xy <- xy[-nrow(xy), , drop = FALSE]
  graphics::polypath(xy[, 1], xy[, 2], col = col, border = NA, rule = "winding")
}

#' A ring of screen points turned anticlockwise or clockwise
#' @noRd
orient_ring <- function(p, anticlockwise) {
  n <- nrow(p)
  twice_area <- sum(p[, 1] * p[c(2:n, 1), 2] - p[c(2:n, 1), 1] * p[, 2])
  if ((twice_area > 0) != anticlockwise) p[n:1, , drop = FALSE] else p
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

#' The solid of an ISEA grid
#'
#' Vertices on the unit sphere, the vertex indices of each face, each face's
#' outward normal, and its edges (30 on the icosahedron) as vertex pairs with
#' the two faces that share them.
#' @noRd
icosa_solid <- function(icosa) {
  s <- cpp_icosa_solid(icosa)
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

#' Cell boundaries of a grid on the sphere and on its solid
#'
#' One closed path per cell, as a matrix with columns cell, face, solid_x/y/z,
#' sphere_x/y/z and plane_x/y. An ISEA grid reads them from its faces; an H3
#' cell edge is a great-circle arc between corners, so H3 paths carry sphere
#' positions only.
#' @noRd
grid_surface_paths <- function(g, cells, step) {
  cells <- surface_cells(g, cells)
  if (is_h3_grid(g)) return(h3_sphere_paths(as.character(cells), step))
  lv <- isea_levels(g@aperture, g@resolution)
  cpp_cell_surface_paths(icosa_arg(g), as.numeric(cells), lv$resolution,
                         lv$aperture, lv$ap_seq, step)
}

#' The cells to draw: those given, or every cell of the grid
#' @noRd
surface_cells <- function(g, cells) {
  if (!is.null(cells)) return(cells)
  if (is_h3_grid(g)) h3_all_cells(g@resolution) else seq_len(grid_n_cells(g))
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
          sphere_z = pts[, 3], plane_x = NA_real_, plane_y = NA_real_)
  })
  do.call(rbind, rows)
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
                       radius_deg * pi / 180 * EARTH_RADIUS_KM * 1000)
  part <- suppressMessages(sf::st_intersection(land, cap))
  part <- part[!sf::st_is_empty(part)]
  if (length(part) == 0L) return(list())
  sfc_polygons(polygons_only(part))
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

#' The polygons of a set of polygon geometries, each a list of lon/lat rings
#' (exterior first, then holes)
#' @noRd
sfc_polygons <- function(x) {
  unlist(lapply(x, function(p) {
    if (inherits(p, "MULTIPOLYGON")) unclass(p) else list(unclass(p))
  }), recursive = FALSE)
}

#' Country outlines as one set of 3D polylines on the unit sphere
#' @noRd
land_outline_points <- function(land) {
  rings <- unlist(sfc_polygons(land), recursive = FALSE)
  P <- do.call(rbind, lapply(rings, function(r) unit_vec(r[, 1], r[, 2])))
  brk <- unlist(lapply(rings, function(r) c(TRUE, rep(FALSE, nrow(r) - 1L))))
  list(P = P, brk = brk)
}

# =============================================================================
# SURFACES
# =============================================================================

#' Draw the grid on the sphere
#' @noRd
draw_sphere <- function(paths, land, view, style, icosa) {
  rim <- horizon_ring(view)
  disc <- project_ring(rim, view)
  if (!is.null(disc)) graphics::polygon(disc, col = style$ocean_fill, border = NA)
  d <- view$dir
  center_ll <- vec_lonlat(d)
  cap_deg <- acos(view$horizon) * 180 / pi - 0.5

  if (!is.null(land)) {
    polys <- lapply(land_in_cap(land, center_ll[1], center_ll[2], cap_deg), function(p) {
      lapply(p, function(r) unit_vec(r[, 1], r[, 2]))
    })
    fill_polygons(polys, view, style$land_fill)
    if (!is.na(style$land_border)) {
      o <- land_outline_points(land)
      near <- faces_camera(o$P, view)
      draw_segments(o$P, view, o$brk, near & c(near[-1], FALSE),
                    style$land_border, style$land_lwd)
    }
  }

  S <- paths[, c("sphere_x", "sphere_y", "sphere_z"), drop = FALSE]
  front <- faces_camera(S, view)
  draw_segments(S, view, c(TRUE, diff(paths[, "cell"]) != 0),
                front & c(front[-1], FALSE), style$grid_border, style$grid_lwd)

  if (style$face_edges) {
    solid <- icosa_solid(icosa)
    V <- solid$vertices
    arcs <- lapply(seq_len(nrow(solid$edges)), function(j) {
      a <- V[solid$edges[j, "v1"], ]
      b <- V[solid$edges[j, "v2"], ]
      rbind(slerp(a, b, 0.5 * pi / 180), b)
    })
    P <- do.call(rbind, arcs)
    brk <- unlist(lapply(arcs, function(a) c(TRUE, rep(FALSE, nrow(a) - 1L))))
    vis <- faces_camera(P, view)
    draw_segments(P, view, brk, vis & c(vis[-1], FALSE), style$edge_col, style$edge_lwd)
  }
  graphics::lines(project(rim, view), col = "#9AA3AB", lwd = 0.8)
}

#' Draw the grid on its solid
#' @noRd
draw_solid <- function(paths, land, view, style, icosa) {
  solid <- icosa_solid(icosa)
  V <- solid$vertices
  front <- icosa_front(solid, view)
  shade <- 0.80 + 0.20 * pmax(0, drop(solid$normals %*% view$light))

  for (f in which(front)) {
    tri <- V[solid$faces[f, ], ]
    flat <- project_ring(tri, view)
    if (!is.null(flat)) {
      graphics::polygon(flat, col = shade_col(style$ocean_fill, shade[f]), border = NA)
    }
    if (!is.null(land)) draw_face_land(f - 1L, tri, land, view, style, shade[f], icosa)
  }

  face <- paths[, "face"] + 1L
  Q <- paths[, c("solid_x", "solid_y", "solid_z"), drop = FALSE]
  draw_segments(Q, view, c(TRUE, diff(paths[, "cell"]) != 0), front[face],
                style$grid_border, style$grid_lwd)

  if (style$face_edges) {
    E <- solid$edges
    vis <- front[E[, "f1"]] | front[E[, "f2"]]
    for (j in which(vis)) {
      graphics::lines(project(V[E[j, c("v1", "v2")], ], view),
                      col = style$edge_col, lwd = style$edge_lwd)
    }
  }
}

#' Land and country outlines on one flat face, in its triangle coordinates
#'
#' Land within a cap around the face's centre is read in the face's triangle
#' coordinates and cut to the triangle there. The cap is wider than the face's
#' circumradius, so its rim never reaches the triangle. Returns the filled
#' polygons (lists of rings) and the outlines (matrices), each NULL when the
#' face holds none.
#' @noRd
face_land <- function(face, tri, land, icosa) {
  c3 <- colMeans(tri)
  c_ll <- vec_lonlat(c3 / sqrt(sum(c3^2)))
  polys <- land_in_cap(land, c_ll[1], c_ll[2], 40)
  if (length(polys) == 0L) return(list(fill = NULL, lines = NULL))

  to_tri <- function(r) {
    cpp_lonlat_to_face_solid(icosa, face, r[, 1], r[, 2])[, c("tx", "ty"), drop = FALSE]
  }
  tri_t <- to_tri(vec_lonlat(tri))
  triangle <- sf::st_sfc(sf::st_polygon(list(rbind(tri_t, tri_t[1, ]))))

  old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  flat <- sf::st_sfc(lapply(polys, function(p) {
    sf::st_polygon(lapply(p, function(r) {
      t_r <- to_tri(r)
      rbind(t_r[-nrow(t_r), , drop = FALSE], t_r[1, ])
    }))
  }))
  flat <- suppressWarnings(sf::st_make_valid(flat))

  filled <- suppressWarnings(sf::st_intersection(flat, triangle))
  filled <- filled[!sf::st_is_empty(filled)]
  fill <- if (length(filled) > 0L) sfc_polygons(polygons_only(filled))

  lines_t <- suppressWarnings(sf::st_intersection(sf::st_boundary(flat), triangle))
  lines_t <- lines_t[!sf::st_is_empty(lines_t)]
  lines <- NULL
  if (length(lines_t) > 0L) {
    lines_t <- lines_t[sf::st_dimension(lines_t) == 1L]
    if (any(sf::st_is(lines_t, "GEOMETRYCOLLECTION"))) {
      lines_t <- sf::st_collection_extract(lines_t, "LINESTRING")
    }
    lines <- unlist(lapply(lines_t, function(l) {
      if (inherits(l, "MULTILINESTRING")) unclass(l) else list(unclass(l))
    }), recursive = FALSE)
  }
  list(fill = fill, lines = lines)
}

#' Land and country outlines on one face of the solid
#' @noRd
draw_face_land <- function(face, tri, land, view, style, shade, icosa) {
  part <- face_land(face, tri, land, icosa)
  to_solid <- function(m) cpp_face_tri_to_solid(icosa, face, m[, 1], m[, 2])
  if (!is.na(style$land_fill) && !is.null(part$fill)) {
    fill_polygons(lapply(part$fill, function(p) lapply(p, to_solid)), view,
                  shade_col(style$land_fill, shade))
  }
  if (!is.na(style$land_border) && !is.null(part$lines)) {
    P <- do.call(rbind, lapply(part$lines, to_solid))
    brk <- unlist(lapply(part$lines, function(l) c(TRUE, rep(FALSE, nrow(l) - 1L))))
    draw_segments(P, view, brk, rep(TRUE, nrow(P)), style$land_border, style$land_lwd)
  }
}

# =============================================================================
# NET
# =============================================================================

#' The faces of the solid laid out flat
#'
#' One entry per face (from 0): its vertices on the unit sphere and its
#' triangle in the PLANE layout.
#' @noRd
net_triangles <- function(icosa) {
  solid <- icosa_solid(icosa)
  lapply(seq_len(nrow(solid$faces)), function(f) {
    tri <- solid$vertices[solid$faces[f, ], ]
    ll <- vec_lonlat(tri)
    t <-cpp_lonlat_to_face_solid(icosa, f - 1L, ll[, 1], ll[, 2])
    list(tri = tri, plane = to_plane(f - 1L, t[, c("tx", "ty"), drop = FALSE], icosa))
  })
}

#' Triangle coordinates of one face in the PLANE layout
#' @noRd
to_plane <- function(face, m, icosa) {
  p <- cpp_icosa_tri_to_plane(icosa, rep(as.integer(face), nrow(m)), m[, 1], m[, 2])
  cbind(p$plane_x, p$plane_y)
}

#' Draw the grid on the unfolded solid
#' @noRd
draw_net <- function(paths, land, tris, style, icosa) {
  for (t in tris) graphics::polygon(t$plane, col = style$ocean_fill, border = NA)

  if (!is.null(land)) {
    for (f in seq_along(tris)) {
      part <- face_land(f - 1L, tris[[f]]$tri, land, icosa)
      if (!is.na(style$land_fill) && !is.null(part$fill)) {
        xy <- do.call(rbind, lapply(part$fill, function(p) {
          do.call(rbind, lapply(seq_along(p), function(k) {
            rbind(orient_ring(to_plane(f - 1L, p[[k]], icosa), anticlockwise = k == 1L), NA)
          }))
        }))
        xy <- xy[-nrow(xy), , drop = FALSE]
        graphics::polypath(xy[, 1], xy[, 2], col = style$land_fill, border = NA,
                           rule = "winding")
      }
      if (!is.na(style$land_border) && !is.null(part$lines)) {
        for (l in part$lines) {
          graphics::lines(to_plane(f - 1L, l, icosa), col = style$land_border,
                          lwd = style$land_lwd)
        }
      }
    }
  }

  n <- nrow(paths)
  if (n >= 2L) {
    s <- which(paths[-1L, "cell"] == paths[-n, "cell"] &
               paths[-1L, "face"] == paths[-n, "face"])
    graphics::segments(paths[s, "plane_x"], paths[s, "plane_y"],
                       paths[s + 1L, "plane_x"], paths[s + 1L, "plane_y"],
                       col = style$grid_border, lwd = style$grid_lwd,
                       lend = "round")
  }

  if (style$face_edges) {
    for (t in tris) {
      graphics::polygon(t$plane, col = NA, border = style$edge_col,
                        lwd = style$edge_lwd)
    }
  }
}
