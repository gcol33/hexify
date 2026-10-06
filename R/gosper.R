# =============================================================================
# Gosper-island and descendant cell shapes
# =============================================================================
#
# Two ways to draw a hexagonal cell by a finer boundary, both of which cover
# exactly the cell's area and tile the sphere.
#
# "gosper": every cell edge is replaced, `depth` times over, by three segments
# 1/sqrt(7) as long, the first and last turned by atan(sqrt(3) / 5) and the
# middle one by 60 degrees more. The three segments sum to the edge and are
# symmetric under a half turn about its midpoint, so the curve is the same
# whichever cell draws it and moves as much area out of the cell as into it.
# Repeating the step with the same turn gives a Gosper island (flowsnake).
# The curve is drawn on the face plane of the projection, where cell edges are
# straight and areas are those on the sphere, so the island keeps the cell's
# area exactly; a curve leaving its face is followed onto the next face by
# unfolding across their shared edge.
#
# "descendants": the outline of the cell's descendants `depth` resolutions
# down, which nests across resolutions. An aperture-7 cell refines into its
# central child and the six around it, but the child lattice of an ISEA7H grid
# turns one way at one resolution and back at the next, so the outline stays
# close to the hexagon rather than approaching a flowsnake.

#' Validate the depth of a Gosper or descendant cell shape
#'
#' @param g HexGridInfo object
#' @param shape "gosper" or "descendants"
#' @param depth Number of replacement steps, or of resolutions down
#' @return Integer depth
#' @noRd
cell_shape_depth <- function(g, shape, depth) {
  if (is_h3_grid(g)) {
    stop(sprintf("shape = \"%s\" needs an ISEA grid", shape))
  }
  if (!is.numeric(depth) || length(depth) != 1L || !is.finite(depth) ||
      depth < 1 || depth != round(depth)) {
    stop("depth must be a single positive whole number")
  }
  depth <- as.integer(depth)
  if (shape != "descendants") {
    return(depth)
  }

  # A mixed sequence's hierarchy is geometric, a parent being the coarser cell
  # holding a cell's centre, and its levels are scaled copies of one lattice
  # rather than rotated ones, so its cells do not nest past one level.
  if (is_mixed_aperture(g@aperture) || as.integer(g@aperture) != 7L) {
    stop(sprintf(
      "shape = \"descendants\" needs a grid of aperture 7, whose cells nest exactly; this grid's aperture is %s",
      g@aperture))
  }
  if (g@resolution + depth > MAX_RESOLUTION) {
    stop(sprintf("resolution %d plus depth %d exceeds the maximum resolution %d",
                 g@resolution, depth, MAX_RESOLUTION))
  }
  depth
}

#' Polygons of ISEA cells drawn as Gosper islands or descendant outlines
#'
#' @param cell_id Numeric vector of cell IDs of grid `g`
#' @param g HexGridInfo object
#' @param shape "gosper" or "descendants"
#' @param depth Validated depth
#' @param tolerance Edge tolerance, as for `isea_cell_rings()`
#' @return An sfc of POLYGON geometries, one per cell ID, in input order
#' @noRd
shaped_cells_to_sfc <- function(cell_id, g, shape, depth,
                                tolerance = CELL_EDGE_TOLERANCE) {
  rings <- switch(shape,
    gosper = gosper_cell_rings(cell_id, g, depth, tolerance),
    descendants = descendant_cell_rings(cell_id, g, depth, tolerance)
  )
  sf::st_sfc(lapply(rings, function(r) {
    sf::st_polygon(lapply(r, lonlat_ring_coords))
  }), crs = grid_crs(g))
}

# -----------------------------------------------------------------------------
# Corners and edges shared between cells
# -----------------------------------------------------------------------------

#' Merge points closer than a tolerance
#'
#' Two points within `eps / 4` of each other in every coordinate fall in the
#' same box of side `eps` in at least one of the eight lattices offset by half
#' a box along each axis, and so does every cluster of such points. Taking,
#' over the eight lattices, the smallest index sharing a point's box labels
#' every cluster by its smallest member. Points further apart than
#' `eps * sqrt(3)` are never merged.
#'
#' @param xyz Numeric matrix of points, one per row
#' @param eps Box side
#' @return Integer vector: for each point, the index of the first point it is
#'   merged with
#' @noRd
snap_points <- function(xyz, eps) {
  n <- nrow(xyz)
  label <- seq_len(n)
  q <- xyz / eps
  for (s in 0:7) {
    shift <- c(s %% 2L, (s %/% 2L) %% 2L, s %/% 4L) / 2
    k1 <- floor(q[, 1] + shift[1])
    k2 <- floor(q[, 2] + shift[2])
    k3 <- floor(q[, 3] + shift[3])
    o <- order(k1, k2, k3)
    start <- c(TRUE, diff(k1[o]) != 0 | diff(k2[o]) != 0 | diff(k3[o]) != 0)
    first <- o[start][cumsum(start)]
    smallest <- integer(n)
    smallest[o] <- first
    label <- pmin(label, smallest)
  }
  label
}

#' Unit vectors of lon/lat points
#' @noRd
lonlat_to_xyz <- function(lon, lat) {
  rad <- pi / 180
  cbind(cos(lat * rad) * cos(lon * rad),
        cos(lat * rad) * sin(lon * rad),
        sin(lat * rad))
}

#' Lon/lat of vectors
#' @noRd
xyz_to_lonlat <- function(xyz) {
  r <- sqrt(rowSums(xyz^2))
  cbind(lon = atan2(xyz[, 2], xyz[, 1]) * 180 / pi,
        lat = asin(pmax(-1, pmin(1, xyz[, 3] / r))) * 180 / pi)
}

#' The edges of a set of cells, with corners shared between cells identified
#'
#' @param rings Open corner rings, one lon/lat matrix per cell, every ring
#'   running the same way round
#' @param n_cells Cell count of the grid the rings belong to, which sets the
#'   distance below which two corners are one
#' @return A list: `corners` (all ring corners stacked), and per directed edge
#'   `from` and `to` (row of `corners` labelling each end), `cell` (ring index)
#'   and `k` (the corner of the ring it starts at)
#' @noRd
ring_edges <- function(rings, n_cells) {
  n_corner <- vapply(rings, nrow, integer(1))
  corners <- do.call(rbind, rings)

  # Corners of one cell lie about a cell spacing apart
  spacing <- sqrt(4 * pi / n_cells)
  label <- snap_points(lonlat_to_xyz(corners[, 1], corners[, 2]), spacing * 1e-3)

  cell <- rep(seq_along(rings), n_corner)
  first <- cumsum(c(0L, n_corner[-length(n_corner)]))
  k <- sequence(n_corner)
  following <- first[cell] + k %% n_corner[cell] + 1L
  list(corners = corners, from = label, to = label[following], cell = cell,
       k = k)
}

#' Open corner rings of ISEA cells
#' @noRd
open_corner_rings <- function(cell_id, g) {
  lapply(isea_cell_rings(cell_id, g@resolution, g@aperture, icosa_arg(g),
                         tolerance = 0),
         function(r) r[-nrow(r), , drop = FALSE])
}

# -----------------------------------------------------------------------------
# Gosper islands
# -----------------------------------------------------------------------------

#' Outline rings of Gosper islands
#'
#' @param cell_id Numeric vector of cell IDs of grid `g`
#' @param g HexGridInfo object of an ISEA grid
#' @param depth Number of replacement steps
#' @param tolerance Edge tolerance, as for `isea_cell_rings()`
#' @return A list with one element per cell ID, each a list holding one closed
#'   lon/lat ring
#' @noRd
gosper_cell_rings <- function(cell_id, g, depth,
                              tolerance = CELL_EDGE_TOLERANCE) {
  icosa <- icosa_arg(g)
  frames <- icosa_face_frames(icosa)

  rings <- open_corner_rings(cell_id, g)
  topo <- ring_edges(rings, grid_n_cells(g))

  # Each edge once, from its lower corner label to its higher
  span <- nrow(topo$corners) + 1
  lo <- pmin(topo$from, topo$to)
  hi <- pmax(topo$from, topo$to)
  key <- lo * span + hi
  unique_key <- unique(key)
  edge <- match(key, unique_key)
  first <- match(unique_key, key)
  a_ll <- topo$corners[lo[first], , drop = FALSE]
  b_ll <- topo$corners[hi[first], , drop = FALSE]

  paths <- gosper_edge_paths(a_ll, b_ll, depth, frames, icosa, tolerance)

  forward <- topo$from == lo
  n_corner <- vapply(rings, nrow, integer(1))
  by_cell <- split(seq_along(edge), rep(seq_along(rings), n_corner))
  lapply(by_cell, function(e) {
    pieces <- lapply(e, function(i) {
      p <- paths[[edge[i]]]
      if (!forward[i]) p <- p[rev(seq_len(nrow(p))), , drop = FALSE]
      p[-nrow(p), , drop = FALSE]
    })
    xy <- do.call(rbind, pieces)
    list(rbind(xy, xy[1, , drop = FALSE]))
  })
}

#' Gosper curves of cell edges, in lon/lat
#'
#' @param a_ll,b_ll Lon/lat matrices of the two ends of each edge
#' @param depth Number of replacement steps
#' @param frames Face frames from `icosa_face_frames()`
#' @param icosa Orientation argument of the C++ layer
#' @param tolerance Edge tolerance; 0 keeps the curve's own vertices alone
#' @return A list of lon/lat matrices, one per edge, from a to b
#' @noRd
gosper_edge_paths <- function(a_ll, b_ll, depth, frames, icosa, tolerance) {
  n_edge <- nrow(a_ll)
  mid <- xyz_to_lonlat(lonlat_to_xyz(a_ll[, 1], a_ll[, 2]) +
                       lonlat_to_xyz(b_ll[, 1], b_ll[, 2]))

  # Every edge is drawn on the face holding its midpoint, its ends carried
  # onto that face's plane across the face edge they may lie beyond
  face <- locate_face(frames, icosa, mid)
  z_mid <- project_to_face(icosa, face, mid)
  z_a <- unfold_onto_face(frames, icosa, face, a_ll)
  z_b <- unfold_onto_face(frames, icosa, face, b_ll)

  z <- gosper_curve(z_a, z_b, depth)
  n_point <- ncol(z)

  pts <- data.frame(
    edge = rep(seq_len(n_edge), n_point),
    s = rep((seq_len(n_point) - 1) / (n_point - 1), each = n_edge),
    z = as.vector(z)
  )
  pts$face <- rep(face, n_point)
  pts$start <- rep(z_mid, n_point)
  ll <- face_plane_to_lonlat(frames, icosa, pts$face, pts$start, pts$z)
  pts$lon <- ll[, 1]
  pts$lat <- ll[, 2]

  if (tolerance > 0) {
    pts <- densify_plane_paths(pts, frames, icosa, tolerance)
  }

  pts <- pts[order(pts$edge, pts$s), ]
  lapply(split(seq_len(nrow(pts)), pts$edge), function(i) {
    cbind(pts$lon[i], pts$lat[i])
  })
}

#' Gosper curve between two points of a plane
#'
#' @param z_a,z_b Complex vectors, the ends of each curve
#' @param depth Number of replacement steps
#' @return Complex matrix, one row per curve, 3^depth + 1 points from z_a to
#'   z_b
#' @noRd
gosper_curve <- function(z_a, z_b, depth) {
  turn <- exp(-1i * atan(sqrt(3) / 5)) / sqrt(7)
  bend <- exp(1i * pi / 3)
  z <- cbind(z_a, z_b)
  for (step in seq_len(depth)) {
    n <- ncol(z)
    p <- z[, -n, drop = FALSE]
    v <- z[, -1L, drop = FALSE] - p
    u1 <- v * turn
    u2 <- u1 * bend
    out <- matrix(0i, nrow(z), 3L * (n - 1L) + 1L)
    at <- 3L * (seq_len(n - 1L) - 1L) + 1L
    out[, at] <- p
    out[, at + 1L] <- p + u1
    out[, at + 2L] <- p + u1 + u2
    out[, 3L * (n - 1L) + 1L] <- z[, n]
    z <- out
  }
  z
}

#' Split plane paths until each lon/lat chord follows the path
#'
#' A piece is halved in the plane while its plane midpoint lies further from
#' the middle of its lon/lat chord than `tolerance` times the chord's length,
#' the test cell edges are densified by.
#'
#' @param pts Data frame of path points: edge, s (position along the edge),
#'   z, face, start, lon, lat
#' @return The data frame with the added points
#' @noRd
densify_plane_paths <- function(pts, frames, icosa, tolerance) {
  pts <- pts[order(pts$edge, pts$s), ]
  n <- nrow(pts)
  same <- pts$edge[-1L] == pts$edge[-n]
  left <- which(same)
  right <- left + 1L

  for (level in seq_len(12L)) {
    if (length(left) == 0L) break
    dlon <- (pts$lon[right] - pts$lon[left] + 180) %% 360 - 180
    over_pole <- abs(abs(pts$lon[right] - pts$lon[left]) - 180) < 1e-3
    keep <- !over_pole
    left <- left[keep]
    right <- right[keep]
    dlon <- dlon[keep]
    if (length(left) == 0L) break

    half <- data.frame(
      edge = pts$edge[left],
      s = (pts$s[left] + pts$s[right]) / 2,
      z = (pts$z[left] + pts$z[right]) / 2,
      face = pts$face[left],
      start = pts$start[left]
    )
    ll <- face_plane_to_lonlat(frames, icosa, half$face, half$start, half$z)
    half$lon <- ll[, 1]
    half$lat <- ll[, 2]

    clat <- (pts$lat[left] + pts$lat[right]) / 2
    coslat <- cos(clat * pi / 180)
    off_lon <- (half$lon - (pts$lon[left] + dlon / 2) + 180) %% 360 - 180
    off <- sqrt((off_lon * coslat)^2 + (half$lat - clat)^2)
    len <- sqrt((dlon * coslat)^2 + (pts$lat[right] - pts$lat[left])^2)
    split_here <- len > 0 & off > tolerance * len
    if (!any(split_here)) break

    added <- nrow(pts) + seq_len(sum(split_here))
    pts <- rbind(pts, half[split_here, ])
    left_new <- c(left[split_here], added)
    right <- c(added, right[split_here])
    left <- left_new
  }
  pts
}

# -----------------------------------------------------------------------------
# Face planes of the icosahedron
# -----------------------------------------------------------------------------

#' The faces of the active polyhedron in their own plane coordinates
#'
#' Each face's plane carries its corners, as complex numbers; across each
#' face edge a neighbour, and the rotation and shift `z -> a z + b` taking a
#' point of this face's plane, unfolded across that edge, to the neighbour's
#' plane.
#'
#' @param icosa Orientation argument of the C++ layer
#' @return List, one row or entry per face: `corner` (complex), `orient`
#'   (sign of the corner order), `neighbour`, `a`, `b` (column k for the edge
#'   from corner k to corner k + 1), `centre` (unit vector) and `reach` (cosine
#'   of the angle from the centre to the corners)
#' @noRd
icosa_face_frames <- function(icosa) {
  solid <- cpp_icosa_solid(icosa)
  vertex <- xyz_to_lonlat(solid$vertices)
  fv <- solid$faces

  n_face <- nrow(fv)
  corner <- matrix(0i, n_face, 3L)
  for (f in seq_len(n_face)) {
    t <- cpp_lonlat_to_face_solid(icosa, f - 1L, vertex[fv[f, ], 1],
                                  vertex[fv[f, ], 2])
    corner[f, ] <- complex(real = t[, "tx"], imaginary = t[, "ty"])
  }
  orient <- sign(Im(Conj(corner[, 2] - corner[, 1]) * (corner[, 3] - corner[, 1])))

  neighbour <- matrix(0L, n_face, 3L)
  a <- matrix(0i, n_face, 3L)
  b <- matrix(0i, n_face, 3L)
  for (f in seq_len(n_face)) for (k in 1:3) {
    k2 <- k %% 3L + 1L
    va <- fv[f, k]
    vb <- fv[f, k2]
    g <- setdiff(which(rowSums(fv == va) + rowSums(fv == vb) == 2L), f)
    qa <- corner[g, match(va, fv[g, ])]
    qb <- corner[g, match(vb, fv[g, ])]
    pa <- corner[f, k]
    pb <- corner[f, k2]
    neighbour[f, k] <- g
    a[f, k] <- (qb - qa) / (pb - pa)
    b[f, k] <- qa - a[f, k] * pa
    third <- a[f, k] * corner[f, 6L - k - k2] + b[f, k]
    if (Mod(third - corner[g, -match(c(va, vb), fv[g, ])]) < 1e-6) {
      stop("hexify internal error: two face planes have opposite orientation")
    }
  }

  centre <- t(vapply(seq_len(n_face), function(f) {
    m <- colSums(solid$vertices[fv[f, ], , drop = FALSE])
    m / sqrt(sum(m^2))
  }, numeric(3)))
  reach <- vapply(seq_len(n_face), function(f) {
    min(solid$vertices[fv[f, ], , drop = FALSE] %*% centre[f, ])
  }, numeric(1))

  list(corner = corner, orient = orient, neighbour = neighbour, a = a, b = b,
       centre = centre, reach = reach)
}

#' Which side of each edge of a face a plane point lies on
#'
#' @return Matrix, one column per face edge: positive inside the face
#' @noRd
face_edge_sides <- function(frames, face, z) {
  c1 <- frames$corner[cbind(face, 1L)]
  c2 <- frames$corner[cbind(face, 2L)]
  c3 <- frames$corner[cbind(face, 3L)]
  side <- function(p, q) Im(Conj(q - p) * (z - p))
  cbind(side(c1, c2), side(c2, c3), side(c3, c1)) * frames$orient[face]
}

#' Plane coordinates of lon/lat points on given faces
#' @noRd
project_to_face <- function(icosa, face, ll) {
  z <- complex(length(face))
  for (f in unique(face)) {
    at <- which(face == f)
    t <- cpp_lonlat_to_face_solid(icosa, f - 1L, ll[at, 1], ll[at, 2])
    z[at] <- complex(real = t[, "tx"], imaginary = t[, "ty"])
  }
  z
}

#' The face holding each lon/lat point
#'
#' @param prefer Optional face per point, taken when it holds the point, then
#'   any neighbour of it that does
#' @return Integer face index per point
#' @noRd
locate_face <- function(frames, icosa, ll, prefer = NULL) {
  n <- nrow(ll)
  xyz <- lonlat_to_xyz(ll[, 1], ll[, 2])
  slack <- 1e-9
  n_face <- nrow(frames$corner)
  holds <- matrix(FALSE, n, n_face)
  for (f in seq_len(n_face)) {
    near <- which(xyz %*% frames$centre[f, ] > frames$reach[f] - 0.1)
    if (length(near) == 0L) next
    z <- project_to_face(icosa, rep(f, length(near)), ll[near, , drop = FALSE])
    holds[near, f] <- apply(face_edge_sides(frames, rep(f, length(near)), z), 1,
                            min) >= -slack
  }
  closeness <- xyz %*% t(frames$centre)
  closeness[!holds] <- -Inf
  face <- max.col(closeness, ties.method = "first")

  if (!is.null(prefer)) {
    own <- holds[cbind(seq_len(n), prefer)]
    face[own] <- prefer[own]
    for (k in 1:3) {
      nb <- frames$neighbour[prefer, k]
      take <- !own & holds[cbind(seq_len(n), nb)]
      face[take] <- nb[take]
      own <- own | take
    }
  }
  if (any(!is.finite(closeness[cbind(seq_len(n), face)]))) {
    stop("hexify internal error: a point lies on no icosahedron face")
  }
  face
}

#' Plane coordinates of lon/lat points on a face, unfolded across its edges
#'
#' A point on the face itself is projected onto it; a point on a neighbouring
#' face is projected there and carried across their shared edge.
#'
#' @param face Face per point
#' @return Complex plane coordinates on `face`
#' @noRd
unfold_onto_face <- function(frames, icosa, face, ll) {
  home <- locate_face(frames, icosa, ll, prefer = face)
  z <- project_to_face(icosa, home, ll)
  away <- which(home != face)
  if (length(away) > 0L) {
    k <- vapply(away, function(i) match(home[i], frames$neighbour[face[i], ]),
                integer(1))
    if (anyNA(k)) {
      stop("hexify internal error: a cell edge spans faces that share no edge")
    }
    idx <- cbind(face[away], k)
    z[away] <- (z[away] - frames$b[idx]) / frames$a[idx]
  }
  z
}

#' Lon/lat of points of a face plane, the plane unfolded across face edges
#'
#' Each point is followed from a point of its starting face along the straight
#' line to it: crossing a face edge carries it onto the neighbour's plane.
#'
#' @param face Starting face per point
#' @param start Complex plane points inside the starting faces
#' @param z Complex plane points to locate
#' @return Lon/lat matrix
#' @noRd
face_plane_to_lonlat <- function(frames, icosa, face, start, z) {
  slack <- 1e-12
  for (step in seq_len(8L)) {
    sides <- face_edge_sides(frames, face, z)
    out <- which(apply(sides, 1, min) < -slack)
    if (length(out) == 0L) break

    # The first edge the line start -> z leaves the face by. The start lies
    # in the face, so the crossing of an edge z lies beyond is at t in (0, 1];
    # a start on an edge can read as just outside it, giving t just below 0.
    f <- face[out]
    m <- start[out]
    d <- z[out] - m
    cross <- function(u, v) Im(Conj(u) * v)
    t_exit <- vapply(1:3, function(k) {
      p <- frames$corner[cbind(f, k)]
      q <- frames$corner[cbind(f, k %% 3L + 1L)]
      t <- pmax(cross(p - m, q - p) / cross(d, q - p), 0)
      t[!(sides[out, k] < -slack) | !is.finite(t)] <- Inf
      t
    }, numeric(length(out)))
    t_exit <- matrix(t_exit, ncol = 3L)
    k <- max.col(-t_exit, ties.method = "first")
    t <- t_exit[cbind(seq_along(out), k)]

    idx <- cbind(f, k)
    a <- frames$a[idx]
    b <- frames$b[idx]
    start[out] <- a * (m + pmin(t, 1) * d) + b
    z[out] <- a * z[out] + b
    face[out] <- frames$neighbour[idx]
  }
  if (any(apply(face_edge_sides(frames, face, z), 1, min) < -slack)) {
    stop("hexify internal error: a curve point could not be placed on a face")
  }

  ll <- vapply(seq_along(z), function(i) {
    cpp_face_xy_to_ll(icosa, Re(z[i]), Im(z[i]), face[i] - 1L)
  }, numeric(2))
  cbind(lon = ll[1, ], lat = ll[2, ])
}

# -----------------------------------------------------------------------------
# Descendant outlines
# -----------------------------------------------------------------------------

#' Outline rings of the union of each cell's descendants
#'
#' The outline of a union of cells is the set of cell edges with a cell of the
#' union on one side only. Every cell ring runs the same way round, so an edge
#' shared by two cells of the union appears once in each direction, and the
#' edges left over chain into rings: three cells meet at every corner, so a
#' corner on the outline starts exactly one outline edge of that union.
#'
#' @param cell_id Numeric vector of cell IDs of grid `g`
#' @param g HexGridInfo object of aperture 7
#' @param depth Number of levels the outline is drawn down
#' @param tolerance Edge tolerance, as for `isea_cell_rings()`
#' @return A list with one element per cell ID, each a list of closed lon/lat
#'   rings: the outline first, then any holes
#' @noRd
descendant_cell_rings <- function(cell_id, g, depth,
                                  tolerance = CELL_EDGE_TOLERANCE) {
  fine_grid <- grid_at_resolution(g, g@resolution + depth)

  kids <- get_children(cell_id, g, levels = depth)
  fine <- unlist(kids, use.names = FALSE)
  owner <- rep(seq_along(cell_id), lengths(kids))

  rings <- open_corner_rings(fine, fine_grid)
  topo <- ring_edges(rings, grid_n_cells(fine_grid))
  from <- topo$from
  to <- topo$to
  island <- owner[topo$cell]

  span <- nrow(topo$corners) + 1
  reverse <- match(to * span + from, from * span + to)
  outline <- is.na(reverse) | island[reverse] != island

  e_from <- from[outline]
  e_to <- to[outline]
  e_island <- island[outline]
  e_row <- which(outline)

  start_key <- e_island * span + e_from
  if (anyDuplicated(start_key)) {
    stop("hexify internal error: an outline meets itself at a corner")
  }
  nxt <- match(e_island * span + e_to, start_key)
  if (anyNA(nxt)) {
    stop("hexify internal error: an outline does not close")
  }

  # Walk every ring at once, each from its smallest edge
  n_edge <- length(nxt)
  root <- seq_len(n_edge)
  jump <- nxt
  for (step in seq_len(ceiling(log2(n_edge)) + 1L)) {
    root <- pmin(root, root[jump])
    jump <- jump[jump]
  }
  starts <- which(root == seq_len(n_edge))
  order_in_ring <- integer(n_edge)
  ring_of <- integer(n_edge)
  current <- starts
  position <- 0L
  active <- rep(TRUE, length(starts))
  while (any(active)) {
    position <- position + 1L
    at <- current[active]
    order_in_ring[at] <- position
    ring_of[at] <- which(active)
    current[active] <- nxt[at]
    active[active] <- current[active] != starts[active]
  }

  # The path of each outline edge from its first corner up to its last
  pieces <- if (tolerance > 0) {
    dense_edge_pieces(fine, topo$cell[e_row], topo$k[e_row], rings, fine_grid,
                      tolerance)
  } else {
    lapply(e_row, function(r) topo$corners[r, , drop = FALSE])
  }

  by_ring <- split(seq_len(n_edge), ring_of)
  ring_island <- vapply(by_ring, function(e) e_island[e[1]], integer(1))
  ring_coords <- lapply(by_ring, function(e) {
    e <- e[order(order_in_ring[e])]
    xy <- do.call(rbind, pieces[e])
    rbind(xy, xy[1, , drop = FALSE])
  })

  lapply(seq_along(cell_id), function(i) {
    own <- ring_coords[ring_island == i]
    own[order(-vapply(own, nrow, integer(1)))]
  })
}

#' Densified paths of cell edges
#'
#' Each cell's densified ring holds its corners, so the path of edge k runs
#' from the ring point at corner k up to the one before corner k + 1.
#'
#' @param fine Cell IDs of the fine grid
#' @param cell,k For each edge, the index into `fine` and the corner it starts
#'   at
#' @param rings Open corner rings of `fine`
#' @param fine_grid HexGridInfo of `fine`
#' @param tolerance Edge tolerance
#' @return A list of lon/lat matrices, one per edge
#' @noRd
dense_edge_pieces <- function(fine, cell, k, rings, fine_grid, tolerance) {
  used <- unique(cell)
  dense <- isea_cell_rings(fine[used], fine_grid@resolution, fine_grid@aperture,
                           icosa_arg(fine_grid), tolerance = tolerance)
  at <- lapply(seq_along(used), function(u) {
    d <- dense[[u]]
    d <- d[-nrow(d), , drop = FALSE]
    dxyz <- lonlat_to_xyz(d[, 1], d[, 2])
    cxyz <- lonlat_to_xyz(rings[[used[u]]][, 1], rings[[used[u]]][, 2])
    idx <- vapply(seq_len(nrow(cxyz)), function(j) {
      which.min(colSums((t(dxyz) - cxyz[j, ])^2))
    }, integer(1))
    list(ring = d, idx = idx)
  })
  slot <- match(cell, used)

  lapply(seq_along(cell), function(e) {
    a <- at[[slot[e]]]
    n <- nrow(a$ring)
    i0 <- a$idx[k[e]]
    i1 <- a$idx[k[e] %% length(a$idx) + 1L]
    len <- (i1 - i0) %% n
    a$ring[(i0 - 1L + seq_len(len) - 1L) %% n + 1L, , drop = FALSE]
  })
}
