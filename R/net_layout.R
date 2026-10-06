# net_layout.R
# Layouts of the unfolded solid: the pieces of its faces a flat map shows and
# where each piece lies in the plane.
#
# A piece is a triangle of one face, given in that face's triangle coordinates
# (the plane of Snyder's or Fuller's face projection): the whole face, or the
# third of it between the face centre and one edge. A plane isometry places
# it: the triangle point t lies at A t + b. Piece corners carry labels, "v<k>"
# for vertex k of the solid and "c<f>" for the centre of face f, so two piece
# edges are the same edge on the sphere when their labels agree, and the map
# is continuous across them when their corners also land on the same points.
# The face projection is applied first and the layout moves its pieces
# rigidly, so an equal-area face projection gives an equal-area map. A piece
# may be placed more than once.

#' Lay out the faces of a grid's solid in the plane
#'
#' Builds a flat map of the sphere from the faces of the solid an ISEA grid is
#' built on: the face projection of the grid (Snyder's equal-area projection,
#' or Fuller's) maps the sphere onto the faces, and the layout cuts the faces
#' into pieces and places each piece in the plane by a rotation and a
#' translation (the rhombic layout also shears). Draw a grid on it with
#' \code{\link[=plot,HexGridInfo,missing-method]{plot}(grid, surface = "net",
#' layout = )}, and place points on it with \code{\link{net_project}}.
#'
#' The layouts are:
#' \describe{
#'   \item{\code{"plane"}}{The faces in the PLANE layout of DGGRID, the layout
#'     \code{\link{hexify_cell_to_plane}} reads: five strips of four triangles
#'     on the icosahedron, four diamonds side by side on the octahedron.}
#'   \item{\code{"gosper"}}{Van de Sande's Gosper World on the octahedron: four
#'     regular hexagons, each one face of the octahedron and a third of each of
#'     its three neighbours, the three neighbouring faces cut along the lines
#'     from their centres to their corners. The tile holding \code{centre} lies
#'     in the middle and the other three are attached on alternate sides.}
#'   \item{\code{"gosper_flower"}}{The tile holding \code{centre} surrounded by
#'     the six tiles that border it, so each of the other three tiles appears
#'     twice (Van de Sande 2024, Fig. 5). The map is continuous across the
#'     sides of the middle tile; the three seams run out from its corners that
#'     are vertices of the octahedron.}
#'   \item{\code{"gosper_land"}}{The four Gosper tiles joined along the three
#'     hexagon sides, out of the twelve, that cross the most land of
#'     \code{land}, so the continents stay whole where they can, as in Van de
#'     Sande's Fig. 1 (2024).}
#'   \item{\code{"rhombic"}}{DGGAL's 5 x 6 rhombic space on the icosahedron
#'     (Jacovella-St-Louis et al. 2025): each face sheared onto half of a unit
#'     square, the ten diamonds a staircase of squares. Areas keep one
#'     constant. DGGAL's y axis points down the page, so a point at DGGAL's
#'     (x, y) is placed at (x, -y).}
#'   \item{\code{"land"}}{A net cut along the solid's edges that cross the
#'     least land of \code{land}, so the faces stay joined across land where
#'     they can, in the spirit of Fuller's Dymaxion map, which also splits two
#'     faces. With \code{hex_grid(orientation = "dymaxion", projection =
#'     "fuller")} it is a Dymaxion-style map of whole faces.}
#' }
#'
#' On the octahedron with Snyder's projection the Gosper World is equal-area,
#' where Van de Sande built it from a gnomonic projection onto a rhombic
#' dodecahedron. His poles lie at midpoints of octahedron edges, which is
#' \code{hex_grid(polyhedron = "octahedron", orientation = "gosper")}; Rus's
#' four-hexagon map puts them at octahedron vertices, which is the standard
#' orientation of the octahedron.
#'
#' @param x A HexGridInfo object from \code{\link{hex_grid}}, built on the
#'   icosahedron or the octahedron.
#' @param layout Name of the layout; see Details.
#' @param centre For the Gosper layouts, the point \code{c(lon, lat)} whose
#'   tile lies in the middle, at the origin. The map is turned so the
#'   hexagons have horizontal top and bottom sides, with north at that point
#'   within 30 degrees of up the page. \code{NULL} uses the north pole.
#' @param mirror For the Gosper layouts, \code{TRUE} builds the hexagons on
#'   the other four faces of the octahedron, the mirror image of the default
#'   tiling.
#' @param land For \code{layout = "land"} and \code{"gosper_land"}, the land
#'   the net keeps joined: \code{TRUE} for the built-in
#'   \code{\link{hexify_world}}, or an sf/sfc object of polygons.
#'
#' @return A \code{hexify_net} object: a list with the layout's \code{name},
#'   the grid's solid in \code{polyhedron}, the \code{icosa} argument of the
#'   grid, and \code{pieces}, one list per placed piece with its \code{face}
#'   (from 0), its \code{region} (three corners in the face's triangle
#'   coordinates), their \code{labels}, the affine map \code{A} and \code{b}
#'   placing a point t at \code{A t + b}, and
#'   the \code{tile} it belongs to.
#'
#' @references
#' Van de Sande, A. (2024). Gosper World: A Hexagonal Map Using Gosper
#' Fractals. \emph{Bridges 2024 Conference Proceedings}, 507-510.
#'
#' Rus, J. (2017). Flowsnake Earth. \emph{Bridges 2017 Conference
#' Proceedings}, 237-244.
#'
#' Jacovella-St-Louis, J., Padilla Ruiz, M., Moreira de Sousa, L.,
#' Peterson, P., Stefanakis, E. (2025). Interoperable global and local
#' indexing of discrete global grid systems based on the ISEA projection for
#' efficient storage, processing and transmission. \emph{Abstracts of the
#' International Cartographic Association} 10, 126.
#' \doi{10.5194/ica-abs-10-126-2025}
#'
#' @seealso \code{\link{net_project}} to place points on a layout
#'
#' @export
#' @examples
#' grid <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
#'                  orientation = "gosper")
#' world <- net_layout(grid, "gosper")
#' plot(grid, surface = "net", layout = world, seams = TRUE)
#' plot(grid, surface = "net", layout = "gosper_flower", land = TRUE)
net_layout <- function(x, layout = c("plane", "gosper", "gosper_flower",
                                     "gosper_land", "rhombic", "land"),
                       centre = NULL, mirror = FALSE, land = TRUE) {
  layout <- match.arg(layout)
  g <- extract_grid(x)
  if (is_h3_grid(g)) {
    stop("a net lays out the solid of an ISEA grid; H3 cells are drawn on ",
         "the sphere", call. = FALSE)
  }
  spec <- NET_LAYOUTS[[layout]]
  polyhedron <- grid_polyhedron(g)
  if (!polyhedron %in% spec$polyhedra) {
    stop(sprintf("layout = \"%s\" lays out the %s; this grid is built on the %s",
                 layout, paste(spec$polyhedra, collapse = " or "), polyhedron),
         call. = FALSE)
  }
  if (!is.null(centre) && !isTRUE(spec$centred)) {
    stop("centre places the middle tile of a Gosper layout", call. = FALSE)
  }
  if (!isFALSE(mirror) && !isTRUE(spec$centred)) {
    stop("mirror picks the hexagons of a Gosper layout", call. = FALSE)
  }
  if (!isTRUE(land) && !isTRUE(spec$joins_land)) {
    stop("land picks what the \"land\" and \"gosper_land\" layouts keep joined",
         call. = FALSE)
  }
  if (!is.null(centre)) centre <- check_lonlat(centre, "centre")
  if (!is.logical(mirror) || length(mirror) != 1L || is.na(mirror)) {
    stop("mirror must be TRUE or FALSE", call. = FALSE)
  }

  icosa <- icosa_arg(g)
  faces <- net_faces(icosa)
  pieces <- spec$build(faces = faces, icosa = icosa, centre = centre,
                       mirror = mirror, land = land)
  structure(list(name = layout, polyhedron = polyhedron, icosa = icosa,
                 pieces = pieces),
            class = "hexify_net")
}

#' @export
print.hexify_net <- function(x, ...) {
  tiles <- length(unique(vapply(x$pieces, `[[`, numeric(1), "tile")))
  cat(sprintf("<hexify_net> \"%s\" layout of the %s: %d pieces in %d tiles\n",
              x$name, x$polyhedron, length(x$pieces), tiles))
  invisible(x)
}

#' Place points on a net layout
#'
#' Projects longitude/latitude points onto their face of the grid's solid
#' and places them where the layout puts that part of the face. A point on a
#' piece the layout shows more than once, or on a cut between two pieces,
#' has one row per place.
#'
#' @param layout A \code{hexify_net} object from \code{\link{net_layout}}
#' @param lon,lat Longitudes and latitudes in degrees
#'
#' @return A data frame with columns \code{point} (position in \code{lon}),
#'   \code{piece} (position in \code{layout$pieces}), \code{x} and \code{y},
#'   in units of a face edge.
#'
#' @export
#' @examples
#' grid <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
#'                  orientation = "gosper")
#' world <- net_layout(grid, "gosper")
#' net_project(world, c(16.37, -74.0), c(48.21, 40.71))
net_project <- function(layout, lon, lat) {
  if (!inherits(layout, "hexify_net")) {
    stop("layout must be a hexify_net object from net_layout()", call. = FALSE)
  }
  if (length(lon) != length(lat)) stop("lon and lat must have the same length", call. = FALSE)
  P <- unit_vec(lon, lat)
  face <- point_faces(P, layout$icosa)
  out <- lapply(seq_along(layout$pieces), function(k) {
    p <- layout$pieces[[k]]
    i <- which(face == p$face)
    if (length(i) == 0L) return(NULL)
    t <- cpp_lonlat_to_face_solid(layout$icosa, p$face, lon[i], lat[i])[, c("tx", "ty"), drop = FALSE]
    inside <- in_triangle(t, p$region)
    if (!any(inside)) return(NULL)
    xy <- place_points(t[inside, , drop = FALSE], p)
    data.frame(point = i[inside], piece = k, x = xy[, 1], y = xy[, 2])
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    return(data.frame(point = integer(0), piece = integer(0), x = numeric(0),
                      y = numeric(0)))
  }
  out <- out[order(out$point, out$piece), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# =============================================================================
# REGISTRY
# =============================================================================

#' The layouts net_layout() builds: the solids each lays out, whether it takes
#' a centre tile, and its builder. A builder takes the faces from net_faces(),
#' the icosa argument, `centre`, `mirror` and `land`, and returns the pieces.
#' @noRd
NET_LAYOUTS <- list(
  plane = list(polyhedra = c("icosahedron", "octahedron"),
               build = function(faces, icosa, ...) layout_plane(faces, icosa)),
  gosper = list(polyhedra = "octahedron", centred = TRUE,
                build = function(faces, icosa, centre, mirror, ...) {
                  layout_gosper(faces, icosa, centre, mirror, "star")
                }),
  gosper_land = list(polyhedra = "octahedron", centred = TRUE, joins_land = TRUE,
                     build = function(faces, icosa, centre, mirror, land) {
                       layout_gosper(faces, icosa, centre, mirror, "land", land)
                     }),
  rhombic = list(polyhedra = "icosahedron",
                 build = function(faces, ...) layout_rhombic(faces)),
  gosper_flower = list(polyhedra = "octahedron", centred = TRUE,
                       build = function(faces, icosa, centre, mirror, ...) {
                         layout_gosper(faces, icosa, centre, mirror, "flower")
                       }),
  land = list(polyhedra = c("icosahedron", "octahedron"), joins_land = TRUE,
              build = function(faces, icosa, land, ...) layout_land(faces, icosa, land))
)

# =============================================================================
# FACES AND PIECES
# =============================================================================

#' The faces of a solid in their triangle coordinates
#'
#' One entry per face (from 0): its vertex indices (from 1), its corners in
#' its triangle coordinates, and its edge neighbours, `nb[k]` being the face
#' (from 0) across the edge from corner k to corner k + 1.
#' @noRd
net_faces <- function(icosa) {
  solid <- icosa_solid(icosa)
  Fv <- solid$faces
  faces <- lapply(seq_len(nrow(Fv)), function(f) {
    v <- Fv[f, ]
    ll <- vec_lonlat(solid$vertices[v, ])
    t <- cpp_lonlat_to_face_solid(icosa, f - 1L, ll[, 1], ll[, 2])
    list(face = f - 1L, verts = v, tri = unname(t[, c("tx", "ty"), drop = FALSE]))
  })
  for (f in seq_along(faces)) {
    v <- faces[[f]]$verts
    faces[[f]]$nb <- vapply(1:3, function(k) {
      e <- c(v[k], v[k %% 3L + 1L])
      hit <- which(vapply(faces, function(o) sum(e %in% o$verts) == 2L, logical(1)))
      hit[hit != f] - 1L
    }, numeric(1))
  }
  faces
}

#' A piece: the whole of face `fc`, or with `edge` = k the third of it between
#' its centre and the edge from corner k to corner k + 1, placed by `iso`
#' @noRd
net_piece <- function(fc, iso, edge = NULL, tile = 0L) {
  if (is.null(edge)) {
    region <- fc$tri
    labels <- paste0("v", fc$verts)
  } else {
    k2 <- edge %% 3L + 1L
    region <- rbind(colMeans(fc$tri), fc$tri[edge, ], fc$tri[k2, ])
    labels <- c(paste0("c", fc$face), paste0("v", fc$verts[c(edge, k2)]))
  }
  list(face = fc$face, region = unname(region), labels = labels,
       A = iso$A, b = iso$b, tile = tile)
}

#' The isometry carrying face `to` across its edge with face `from`, placed
#' by `iso`, so the shared corners land where `from` puts them: the two faces
#' unfold about their common edge
#' @noRd
unfold_iso <- function(from, iso, to) {
  shared <- intersect(from$verts, to$verts)
  d <- place_points(from$tri[match(shared, from$verts), , drop = FALSE], iso)
  s <- to$tri[match(shared, to$verts), , drop = FALSE]
  isometry_of(s[1, ], s[2, ], d[1, ], d[2, ])
}

#' The rotation and translation carrying s1 to d1 and the direction s1 -> s2
#' along d1 -> d2
#' @noRd
isometry_of <- function(s1, s2, d1, d2) {
  a <- atan2(d2[2] - d1[2], d2[1] - d1[1]) - atan2(s2[2] - s1[2], s2[1] - s1[1])
  A <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2L)
  list(A = A, b = d1 - drop(A %*% s1))
}

#' The identity placement
#' @noRd
identity_iso <- function() list(A = diag(2), b = c(0, 0))

#' Points (rows) placed by an isometry or a piece
#' @noRd
place_points <- function(m, iso) {
  out <- m %*% t(iso$A)
  out[, 1] <- out[, 1] + iso$b[1]
  out[, 2] <- out[, 2] + iso$b[2]
  out
}

#' Isometry `outer` applied after `inner`
#' @noRd
compose_iso <- function(outer, inner) {
  list(A = outer$A %*% inner$A, b = drop(outer$A %*% inner$b) + outer$b)
}

#' Which points (rows of `m`) lie in a triangle, with slack `tol`
#' @noRd
in_triangle <- function(m, tri, tol = 1e-9) {
  side <- function(a, b) {
    (b[1] - a[1]) * (m[, 2] - a[2]) - (b[2] - a[2]) * (m[, 1] - a[1])
  }
  s <- sign((tri[2, 1] - tri[1, 1]) * (tri[3, 2] - tri[1, 2]) -
              (tri[2, 2] - tri[1, 2]) * (tri[3, 1] - tri[1, 1]))
  s * side(tri[1, ], tri[2, ]) >= -tol & s * side(tri[2, ], tri[3, ]) >= -tol &
    s * side(tri[3, ], tri[1, ]) >= -tol
}

#' The face (from 0) each unit vector (row) lies on
#'
#' The faces of a regular solid are the Voronoi cells of their centres, so a
#' point lies on the face whose centre is nearest.
#' @noRd
point_faces <- function(P, icosa) {
  N <- icosa_solid(icosa)$normals
  max.col(P %*% t(N), ties.method = "first") - 1L
}

#' Turn and shift a layout so (lon, lat) sits at the origin and north there
#' points up; with `snap`, the turn goes on to the nearest one that makes
#' the first piece's edges run at 30 degrees to the horizontal, so a hexagon
#' around that piece has horizontal top and bottom sides
#' @noRd
orient_layout <- function(pieces, icosa, lon, lat, snap = FALSE) {
  probe <- rbind(c(lon, lat), c(lon, lat + if (lat > 89) -0.01 else 0.01))
  layout <- list(icosa = icosa, pieces = pieces)
  class(layout) <- "hexify_net"
  xy <- net_project(layout, probe[, 1], probe[, 2])
  a <- xy[xy$point == 1L, , drop = FALSE][1, ]
  b <- xy[xy$point == 2L & xy$piece == a$piece, , drop = FALSE]
  if (nrow(b) == 0L) b <- xy[xy$point == 2L, , drop = FALSE]
  heading <- atan2(b$y[1] - a$y, b$x[1] - a$x)
  # At the north pole "north" is the meridian of `lon` continued over it.
  turn <- pi / 2 - heading + if (lat > 89) pi else 0
  if (snap) {
    edge <- place_points(pieces[[1]]$region[1:2, , drop = FALSE], pieces[[1]])
    now <- atan2(edge[2, 2] - edge[1, 2], edge[2, 1] - edge[1, 1]) + turn
    turn <- turn + ((pi / 6 - now + pi / 6) %% (pi / 3)) - pi / 6
  }
  rot <- list(A = matrix(c(cos(turn), sin(turn), -sin(turn), cos(turn)), 2L),
              b = c(0, 0))
  shift <- list(A = diag(2), b = -drop(rot$A %*% c(a$x, a$y)))
  move <- compose_iso(shift, rot)
  lapply(pieces, function(p) {
    m <- compose_iso(move, p)
    p$A <- m$A
    p$b <- m$b
    p
  })
}

# =============================================================================
# LAYOUTS
# =============================================================================

#' DGGRID's PLANE layout: each whole face where face_tri_to_plane() puts it
#' @noRd
layout_plane <- function(faces, icosa) {
  lapply(faces, function(fc) {
    p <- cpp_icosa_tri_to_plane(icosa, rep(as.integer(fc$face), 3L), c(0, 1, 0),
                                c(0, 0, 1))
    A <- cbind(c(p$plane_x[2] - p$plane_x[1], p$plane_y[2] - p$plane_y[1]),
               c(p$plane_x[3] - p$plane_x[1], p$plane_y[3] - p$plane_y[1]))
    net_piece(fc, list(A = A, b = c(p$plane_x[1], p$plane_y[1])), tile = fc$face)
  })
}

#' A net of whole faces unfolded along a spanning tree of the faces
#'
#' `parent[f]` is the face (from 0) face f - 1 unfolds from, NA at the root.
#' @noRd
layout_tree <- function(faces, parent) {
  iso <- vector("list", length(faces))
  root <- which(is.na(parent))
  iso[[root]] <- identity_iso()
  queue <- root
  while (length(queue) > 0L) {
    f <- queue[1]
    queue <- queue[-1]
    kids <- which(parent == f - 1L)
    for (k in kids) iso[[k]] <- unfold_iso(faces[[f]], iso[[f]], faces[[k]])
    queue <- c(queue, kids)
  }
  lapply(seq_along(faces), function(f) net_piece(faces[[f]], iso[[f]], tile = f - 1L))
}

#' The share of each great-circle arc from `from[i, ]` to `to[i, ]` (unit
#' vectors) that runs over land, from points every quarter degree along it
#' @noRd
arc_land_share <- function(from, to, land) {
  land <- surface_land(land)
  if (is.null(land)) stop("land must be TRUE or an sf/sfc object", call. = FALSE)
  probe <- lapply(seq_len(nrow(from)), function(j) {
    slerp(from[j, ], to[j, ], 0.25 * pi / 180)
  })
  ll <- vec_lonlat(do.call(rbind, probe))
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  pts <- sf::st_as_sf(data.frame(lon = ll[, 1], lat = ll[, 2]),
                      coords = c("lon", "lat"), crs = 4326)
  on_land <- lengths(sf::st_intersects(pts, sf::st_union(land))) > 0L
  as.vector(tapply(on_land, rep(seq_along(probe), vapply(probe, nrow, integer(1))), mean))
}

#' The share of each edge of the solid (rows of icosa_solid()$edges) that
#' runs over land
#' @noRd
edge_land_share <- function(icosa, land) {
  solid <- icosa_solid(icosa)
  V <- solid$vertices
  arc_land_share(V[solid$edges[, "v1"], , drop = FALSE],
                 V[solid$edges[, "v2"], , drop = FALSE], land)
}

#' The maximum spanning tree of a graph on nodes 1..n (Kruskal, ties to the
#' lower link), as the order in which a breadth-first walk from `root`
#' reaches the nodes: a data frame of node, the node it is reached from and
#' the link (row of from/to/weight) it is reached by, the root first with NA
#' @noRd
max_spanning_tree <- function(n, from, to, weight, root = 1L) {
  group <- seq_len(n)
  find <- function(i) {
    while (group[i] != i) i <- group[i]
    i
  }
  tree <- logical(length(weight))
  for (j in order(-weight, seq_along(weight))) {
    a <- find(from[j])
    b <- find(to[j])
    if (a == b) next
    group[a] <- b
    tree[j] <- TRUE
  }
  walk <- data.frame(node = root, parent = NA_integer_, link = NA_integer_)
  seen <- seq_len(n) == root
  i <- 1L
  while (i <= nrow(walk)) {
    f <- walk$node[i]
    for (j in which(tree & (from == f | to == f))) {
      k <- if (from[j] == f) to[j] else from[j]
      if (seen[k]) next
      seen[k] <- TRUE
      walk <- rbind(walk, data.frame(node = k, parent = f, link = j))
    }
    i <- i + 1L
  }
  walk
}

#' A net cut along the edges that cross the least land
#'
#' Each edge of the solid weighs the share of it that runs over land; the
#' faces unfold along the maximum spanning tree of those weights, so the cuts
#' fall at sea where they can. The net is turned to lie along its longest
#' extent.
#' @noRd
layout_land <- function(faces, icosa, land) {
  E <- icosa_solid(icosa)$edges
  walk <- max_spanning_tree(length(faces), E[, "f1"], E[, "f2"],
                            edge_land_share(icosa, land))
  parent <- rep(NA_real_, length(faces))
  parent[walk$node[-1]] <- walk$parent[-1] - 1L
  pieces <- layout_tree(faces, parent)

  # Turn the longest extent of the net (principal axis of its corners) along x.
  xy <- do.call(rbind, lapply(pieces, function(p) place_points(p$region, p)))
  ax <- eigen(stats::cov(xy), symmetric = TRUE)$vectors[, 1]
  turn <- -atan2(ax[2], ax[1])
  rot <- list(A = matrix(c(cos(turn), sin(turn), -sin(turn), cos(turn)), 2L),
              b = c(0, 0))
  lo <- apply(place_points(xy, rot), 2, min)
  move <- compose_iso(list(A = diag(2), b = -lo), rot)
  lapply(pieces, function(p) {
    m <- compose_iso(move, p)
    p$A <- m$A
    p$b <- m$b
    p
  })
}

#' The corners of each icosahedron face in DGGAL's 5 x 6 rhombic space, one
#' row per face (from 0) as x1, y1, x2, y2, x3, y3 for its vertices in the
#' order of the solid's face table, which is DGGAL's (vertices5x6 and
#' icoIndices in src/projections/ri5x6.ec, https://github.com/ecere/dggal).
#' Each unit square is a diamond of two faces, the ten squares a staircase.
#' @noRd
RHOMBIC_5X6 <- matrix(c(
  1, 0, 0, 0, 1, 1,   2, 1, 1, 1, 2, 2,   3, 2, 2, 2, 3, 3,   4, 3, 3, 3, 4, 4,
  5, 4, 4, 4, 5, 5,   0, 1, 1, 1, 0, 0,   1, 2, 2, 2, 1, 1,   2, 3, 3, 3, 2, 2,
  3, 4, 4, 4, 3, 3,   4, 5, 5, 5, 4, 4,   1, 1, 0, 1, 1, 2,   2, 2, 1, 2, 2, 3,
  3, 3, 2, 3, 3, 4,   4, 4, 3, 4, 4, 5,   5, 5, 4, 5, 5, 6,   0, 2, 1, 2, 0, 1,
  1, 3, 2, 3, 1, 2,   2, 4, 3, 4, 2, 3,   3, 5, 4, 5, 3, 4,   4, 6, 5, 6, 4, 5
), ncol = 6L, byrow = TRUE)

#' The affine map carrying the three points (rows) of `s` to those of `d`
#' @noRd
affine_of <- function(s, d) {
  A <- cbind(d[2, ] - d[1, ], d[3, ] - d[1, ]) %*%
    solve(cbind(s[2, ] - s[1, ], s[3, ] - s[1, ]))
  list(A = A, b = d[1, ] - drop(A %*% s[1, ]))
}

#' DGGAL's 5 x 6 rhombic layout of the icosahedron
#'
#' Each face is sheared onto half of a unit square, so the map keeps areas to
#' one constant and the ten diamonds become a staircase of unit squares.
#' DGGAL's y axis points down the page; the layout places a point at (x, -y)
#' so the map is not mirrored.
#' @noRd
layout_rhombic <- function(faces) {
  lapply(faces, function(fc) {
    d <- matrix(RHOMBIC_5X6[fc$face + 1L, ], 3L, 2L, byrow = TRUE)
    d[, 2] <- -d[, 2]
    net_piece(fc, affine_of(fc$tri, d), tile = fc$face)
  })
}

#' Van de Sande's Gosper World on the octahedron
#'
#' The faces of the octahedron fall into two sets of four, each face bordered
#' only by faces of the other set. Each face of one set is the middle of a
#' regular hexagon that takes a third of each of its three neighbours, so the
#' four hexagons hold the sphere once. Any two hexagons share one corner, a
#' vertex w of the octahedron where both have a 120-degree angle, and the two
#' sides there, each running from w to the centre of a face n of the other
#' set. A tile is joined to another across such a side by unfolding n, whose
#' thirds meet along it.
#'
#' The middle tile is the one holding `centre`. `arrange = "star"` joins the
#' other three across one side at each corner of the middle tile; "flower"
#' joins a copy across each of its six sides; "land" joins the four tiles
#' along the maximum spanning tree of the land along the twelve sides.
#' @noRd
layout_gosper <- function(faces, icosa, centre, mirror, arrange, land = TRUE) {
  if (is.null(centre)) centre <- c(0, 90)
  # Two-colour the faces; the hexagons sit on the set of face 0, or the other.
  colour <- rep(NA, length(faces))
  colour[1] <- TRUE
  queue <- 1L
  while (length(queue) > 0L) {
    f <- queue[1]
    queue <- queue[-1]
    for (n in faces[[f]]$nb + 1L) {
      if (is.na(colour[n])) {
        colour[n] <- !colour[f]
        queue <- c(queue, n)
      }
    }
  }
  middle <- colour != mirror

  # The tile holding the centre: its face, or the hexagon its third belongs to.
  P <- unit_vec(centre[1], centre[2])
  f0 <- point_faces(P, icosa) + 1L
  if (!middle[f0]) {
    fc <- faces[[f0]]
    t <- cpp_lonlat_to_face_solid(icosa, fc$face, centre[1], centre[2])[, c("tx", "ty"), drop = FALSE]
    k <- which(vapply(1:3, function(k) {
      in_triangle(t, rbind(colMeans(fc$tri), fc$tri[k, ], fc$tri[k %% 3L + 1L, ]))
    }, logical(1)))[1]
    f0 <- fc$nb[k] + 1L
  }

  tile <- function(f, iso) {
    fc <- faces[[f]]
    c(list(net_piece(fc, iso, tile = fc$face)),
      lapply(1:3, function(k) {
        n <- faces[[fc$nb[k] + 1L]]
        edge <- which(n$nb == fc$face)
        net_piece(n, unfold_iso(fc, iso, n), edge = edge, tile = fc$face)
      }))
  }
  # The other tile on the side from vertex w to the centre of face n (both
  # from 1), and its placement when tile f lies at `iso`.
  other_tile <- function(f, w, n) {
    nb <- faces[[n]]
    nb$nb[vapply(1:3, function(j) {
      w %in% nb$verts[c(j, j %% 3L + 1L)] && nb$nb[j] + 1L != f
    }, logical(1))] + 1L
  }
  across <- function(f, iso, w, n) {
    o <- other_tile(f, w, n)
    n_iso <- unfold_iso(faces[[f]], iso, faces[[n]])
    list(face = o, iso = unfold_iso(faces[[n]], n_iso, faces[[o]]))
  }
  # Edge k of a face runs from its corner k to corner k + 1, edge k - 1 into
  # corner k: the two sides of the hexagon at that corner.
  prev_edge <- function(k) (k + 1L) %% 3L + 1L

  home <- identity_iso()
  pieces <- tile(f0, home)
  fc <- faces[[f0]]
  if (arrange %in% c("star", "flower")) {
    for (k in 1:3) {
      edges <- if (arrange == "flower") c(k, prev_edge(k)) else k
      for (e in edges) {
        t <- across(f0, home, fc$verts[k], fc$nb[e] + 1L)
        pieces <- c(pieces, tile(t$face, t$iso))
      }
    }
  } else {
    # The twelve sides: (middle face, vertex, face whose centre ends it),
    # each listed once from both of its tiles.
    sides <- do.call(rbind, lapply(which(middle), function(f) {
      do.call(rbind, lapply(1:3, function(k) {
        n <- faces[[f]]$nb[c(k, prev_edge(k))] + 1L
        cbind(f = f, w = faces[[f]]$verts[k], n = n, o = vapply(n, function(m) {
          other_tile(f, faces[[f]]$verts[k], m)
        }, numeric(1)))
      }))
    }))
    sides <- sides[sides[, "f"] < sides[, "o"], , drop = FALSE]
    solid <- icosa_solid(icosa)
    weight <- arc_land_share(solid$vertices[sides[, "w"], , drop = FALSE],
                             solid$normals[sides[, "n"], , drop = FALSE], land)
    walk <- max_spanning_tree(length(faces), sides[, "f"], sides[, "o"], weight, f0)
    iso <- list()
    iso[[f0]] <- home
    for (i in seq_len(nrow(walk))[-1]) {
      s <- sides[walk$link[i], ]
      t <- across(walk$parent[i], iso[[walk$parent[i]]], s[["w"]], s[["n"]])
      iso[[t$face]] <- t$iso
      pieces <- c(pieces, tile(t$face, t$iso))
    }
  }
  orient_layout(pieces, icosa, centre[1], centre[2], snap = TRUE)
}

# =============================================================================
# DRAWING HELPERS
# =============================================================================

#' The edges of a layout's pieces: for each, the piece, its corner labels,
#' its two ends in the plane, whether it is an edge of the solid, and whether
#' the map is continuous across it (another piece has the same edge at the
#' same place)
#' @noRd
net_edges <- function(layout, tol = 1e-9) {
  rows <- lapply(seq_along(layout$pieces), function(i) {
    p <- layout$pieces[[i]]
    xy <- place_points(p$region, p)
    k2 <- c(2L, 3L, 1L)
    data.frame(piece = i, a = p$labels, b = p$labels[k2],
               x0 = xy[, 1], y0 = xy[, 2], x1 = xy[k2, 1], y1 = xy[k2, 2])
  })
  e <- do.call(rbind, rows)
  e$key <- ifelse(e$a < e$b, paste(e$a, e$b), paste(e$b, e$a))
  e$solid_edge <- startsWith(e$a, "v") & startsWith(e$b, "v")
  e$joined <- vapply(seq_len(nrow(e)), function(i) {
    o <- which(e$key == e$key[i] & e$piece != e$piece[i])
    if (length(o) == 0L) return(FALSE)
    same <- abs(e$x0[o] - e$x0[i]) < tol & abs(e$y0[o] - e$y0[i]) < tol &
      abs(e$x1[o] - e$x1[i]) < tol & abs(e$y1[o] - e$y1[i]) < tol
    swap <- abs(e$x0[o] - e$x1[i]) < tol & abs(e$y0[o] - e$y1[i]) < tol &
      abs(e$x1[o] - e$x0[i]) < tol & abs(e$y1[o] - e$y0[i]) < tol
    any(same | swap)
  }, logical(1))
  e
}

#' Glue tabs on the seams of a net, to fold a printed net into the solid
#'
#' The seams are the piece edges across which the map is not continuous.
#' Each edge of the sphere along a seam is printed on two (or more) pieces,
#' and one of them gets a tab: a trapezoid on the side away from its piece,
#' as deep as a fifth of the edge with sides at 45 degrees. The tab goes on
#' the first copy, in the order of the pieces, whose tab covers no piece,
#' or on the first copy when every tab would. Returns one closed ring (rows)
#' per tab.
#' @noRd
net_tabs <- function(layout) {
  e <- net_edges(layout)
  e <- e[!e$joined, , drop = FALSE]
  placed <- lapply(layout$pieces, function(p) place_points(p$region, p))
  tab_of <- function(i) {
    p0 <- c(e$x0[i], e$y0[i])
    p1 <- c(e$x1[i], e$y1[i])
    inside <- colMeans(placed[[e$piece[i]]])
    len <- sqrt(sum((p1 - p0)^2))
    along <- (p1 - p0) / len
    out <- c(-along[2], along[1])
    if (sum(out * (inside - p0)) > 0) out <- -out
    h <- 0.2 * len
    rbind(p0, p0 + h * (out + along), p1 + h * (out - along), p1, p0)
  }
  covers <- function(tab) {
    probe <- rbind(tab[2:3, ], colMeans(tab[1:4, ]),
                   (tab[2, ] + tab[3, ]) / 2)
    any(vapply(placed, function(r) any(in_triangle(probe, r, tol = -1e-9)),
               logical(1)))
  }
  lapply(split(seq_len(nrow(e)), factor(e$key, levels = unique(e$key))), function(rows) {
    tabs <- lapply(rows, tab_of)
    free <- which(!vapply(tabs, covers, logical(1)))
    tabs[[if (length(free)) free[1] else 1L]]
  })
}

#' Draw glue tabs
#' @noRd
draw_net_tabs <- function(tabs) {
  for (t in tabs) {
    graphics::polygon(t, col = "#E6E9EC", border = "#8C969F", lwd = 0.8)
  }
}

#' Clip segments to a triangle (Cyrus-Beck), with slack `tol`
#'
#' `p0` and `p1` are the segments' ends (rows). Returns the kept parts as a
#' matrix with columns x0, y0, x1, y1.
#' @noRd
clip_segments <- function(p0, p1, tri, tol = 1e-9) {
  n <- nrow(p0)
  t0 <- rep(0, n)
  t1 <- rep(1, n)
  s <- sign((tri[2, 1] - tri[1, 1]) * (tri[3, 2] - tri[1, 2]) -
              (tri[2, 2] - tri[1, 2]) * (tri[3, 1] - tri[1, 1]))
  for (k in 1:3) {
    a <- tri[k, ]
    e <- tri[k %% 3L + 1L, ] - a
    f0 <- s * (e[1] * (p0[, 2] - a[2]) - e[2] * (p0[, 1] - a[1])) + tol
    f1 <- s * (e[1] * (p1[, 2] - a[2]) - e[2] * (p1[, 1] - a[1])) + tol
    d <- f1 - f0
    cut <- -f0 / d
    enter <- d > 0
    leave <- d < 0
    t0[enter] <- pmax(t0[enter], cut[enter])
    t1[leave] <- pmin(t1[leave], cut[leave])
    t1[d == 0 & f0 < 0] <- -1
  }
  keep <- t1 >= t0
  d <- p1[keep, , drop = FALSE] - p0[keep, , drop = FALSE]
  cbind(p0[keep, , drop = FALSE] + t0[keep] * d, p0[keep, , drop = FALSE] + t1[keep] * d)
}

#' Polylines on the faces (a path matrix with columns cell, face, tx, ty) as
#' segments of a layout: each run of points on one face is clipped to each
#' piece of that face and placed. Returns a matrix with columns x0, y0, x1, y1.
#' @noRd
net_segments <- function(paths, layout) {
  n <- nrow(paths)
  if (n < 2L) return(matrix(numeric(0), 0L, 4L))
  s <- which(paths[-1L, "cell"] == paths[-n, "cell"] &
               paths[-1L, "face"] == paths[-n, "face"])
  out <- lapply(layout$pieces, function(p) {
    i <- s[paths[s, "face"] == p$face]
    if (length(i) == 0L) return(NULL)
    p0 <- paths[i, c("tx", "ty"), drop = FALSE]
    p1 <- paths[i + 1L, c("tx", "ty"), drop = FALSE]
    seg <- if (all(startsWith(p$labels, "v"))) cbind(p0, p1) else clip_segments(p0, p1, p$region)
    cbind(place_points(seg[, 1:2, drop = FALSE], p), place_points(seg[, 3:4, drop = FALSE], p))
  })
  out <- do.call(rbind, out)
  if (is.null(out)) matrix(numeric(0), 0L, 4L) else unname(out)
}

#' Check a c(lon, lat) argument
#' @noRd
check_lonlat <- function(x, what) {
  if (!is.numeric(x) || length(x) != 2L || !all(is.finite(x)) ||
      x[2] < -90 || x[2] > 90) {
    stop(what, " must be c(lon, lat) in degrees", call. = FALSE)
  }
  unname(x)
}
