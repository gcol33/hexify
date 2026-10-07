# net_cells.R
# Cells on a net layout: the parts a cell is cut into by the pieces of the
# net, and which of them the map shows joined.

#' The parts of cells on a net layout
#'
#' Cuts each cell into its parts on the pieces of a net layout and places
#' them in the plane. A cell on one piece has one part. A cell that crosses
#' from one piece to another stays whole on the map where the two pieces are
#' joined, and is split where a seam of the net runs through it: its parts
#' then lie apart, and Carr et al. (1997) call it a cloned cell. A cell on a
#' piece the layout shows twice also has a part on each copy.
#'
#' Cell walls are straight on the faces of the solid and the layout moves
#' pieces rigidly (or, for \code{"rhombic"}, by one shear that keeps areas
#' to one constant), so a part's map area is exact. Snyder's projection is
#' equal-area, so on a grid built on it the map areas of a cell's parts add
#' up to the cell's share of the sphere: a face of the net has area
#' \eqn{\sqrt{3}/4}, and \eqn{4\pi/n} of the unit sphere for a solid of
#' \eqn{n} faces.
#'
#' @param layout A \code{hexify_net} object from \code{\link{net_layout}}.
#' @param grid The HexGridInfo object the layout was built for.
#' @param cells Cell IDs; \code{NULL} for every cell of the grid.
#'
#' @return A data frame with one row per part: \code{cell_id}, \code{piece}
#'   (position in \code{layout$pieces}), \code{group} (parts of a cell that
#'   the map joins share a group, numbered from 1 within the cell),
#'   \code{area} (map area in units of a face edge squared) and the centroid
#'   \code{x}, \code{y}. A cell whose parts fall into more than one group is
#'   cut by a seam.
#'
#' @references Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA
#'   discrete global grids. Statistical Computing & Graphics Newsletter
#'   8(2/3): 31-39.
#'
#' @seealso \code{\link{net_layout}}, \code{\link{net_project}}
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' net <- net_layout(g, "plane")
#' parts <- net_cells(net, g)
#' cut <- tapply(parts$group, parts$cell_id, max) > 1
#' sum(cut)
#'
#' # On Snyder's equal-area projection every hexagon has the same map area,
#' # however the net cuts it
#' area <- tapply(parts$area, parts$cell_id, sum)
#' hexagon <- !is_pentagon(as.numeric(names(area)), g)
#' range(area[hexagon])
net_cells <- function(layout, grid, cells = NULL) {
  if (!inherits(layout, "hexify_net")) {
    stop("layout must be a hexify_net object from net_layout()", call. = FALSE)
  }
  g <- extract_grid(grid)
  resolve_layout(layout, g)
  cells <- grid_cells(g, cells)
  paths <- grid_surface_paths(g, cells, step = 1)
  joins <- net_joins(layout)

  rows <- split(seq_len(nrow(paths)), paths[, "cell"])
  out <- lapply(names(rows), function(k) {
    p <- paths[rows[[k]], , drop = FALSE]
    parts <- cell_net_parts(p, layout)
    if (is.null(parts)) return(NULL)
    data.frame(cell_id = rep(cells[as.integer(k)], length(parts$piece)),
               piece = parts$piece,
               group = join_groups(parts, joins), area = parts$area,
               x = parts$x, y = parts$y)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' The face triangle in every face's triangle coordinates, anticlockwise
#' @noRd
FACE_TRIANGLE <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2))

#' A point's place along the face triangle's boundary, anticlockwise from
#' corner (0, 0): 0 to 1 along the first edge, 1 to 2 the second, 2 to 3 the
#' third
#' @noRd
boundary_position <- function(pt) {
  V <- FACE_TRIANGLE
  d <- vapply(1:3, function(e) {
    a <- V[e, ]
    b <- V[e %% 3L + 1L, ]
    t <- sum((pt - a) * (b - a))
    c(e - 1 + t, sqrt(sum((a + t * (b - a) - pt)^2)))
  }, numeric(2))
  d[1, which.min(d[2, ])]
}

#' The boundary of face part of a cell: the cell's run of points on the
#' face, from where it enters the face to where it leaves, closed along the
#' face's boundary anticlockwise from the exit back to the entry
#' @noRd
close_on_face <- function(run) {
  n <- nrow(run)
  if (n < 2L) return(matrix(numeric(0), 0L, 2L))
  from <- boundary_position(run[n, ])
  to <- boundary_position(run[1, ])
  if (to <= from) to <- to + 3
  lo <- ceiling(from + 1e-12)
  hi <- floor(to - 1e-12)
  corners <- if (lo <= hi) seq(lo, hi) else numeric(0)
  rbind(run, FACE_TRIANGLE[(corners %% 3) + 1L, , drop = FALSE])
}

#' Whether points (rows) all lie on one edge of the face triangle
#' @noRd
on_one_edge <- function(m, tol = 1e-9) {
  V <- FACE_TRIANGLE
  any(vapply(1:3, function(e) {
    a <- V[e, ]
    d <- V[e %% 3L + 1L, ] - a
    all(abs(d[1] * (m[, 2] - a[2]) - d[2] * (m[, 1] - a[1])) < tol)
  }, logical(1)))
}

#' A polygon (rows) clipped to a triangle (Sutherland-Hodgman)
#' @noRd
clip_polygon <- function(poly, tri) {
  s <- sign((tri[2, 1] - tri[1, 1]) * (tri[3, 2] - tri[1, 2]) -
              (tri[2, 2] - tri[1, 2]) * (tri[3, 1] - tri[1, 1]))
  for (k in 1:3) {
    if (nrow(poly) == 0L) break
    a <- tri[k, ]
    e <- tri[k %% 3L + 1L, ] - a
    f <- s * (e[1] * (poly[, 2] - a[2]) - e[2] * (poly[, 1] - a[1]))
    n <- nrow(poly)
    nxt <- c(seq_len(n)[-1L], 1L)
    keep <- list()
    for (i in seq_len(n)) {
      j <- nxt[i]
      if (f[i] >= 0) keep[[length(keep) + 1L]] <- poly[i, ]
      if ((f[i] >= 0) != (f[j] >= 0)) {
        t <- f[i] / (f[i] - f[j])
        keep[[length(keep) + 1L]] <- poly[i, ] + t * (poly[j, ] - poly[i, ])
      }
    }
    poly <- if (length(keep)) do.call(rbind, keep) else matrix(numeric(0), 0L, 2L)
  }
  poly
}

#' Signed area of a polygon (rows), positive anticlockwise
#' @noRd
polygon_area <- function(p) {
  n <- nrow(p)
  if (n < 3L) return(0)
  nxt <- c(seq_len(n)[-1L], 1L)
  sum(p[, 1] * p[nxt, 2] - p[nxt, 1] * p[, 2]) / 2
}

#' The parts of one cell on the pieces of a layout: its boundary path (rows
#' of cpp_cell_surface_paths()) cut into the run on each face, each run
#' closed on its face and clipped to the pieces of that face. Returns the
#' piece, map area, centroid and placed polygon of each part.
#' @noRd
cell_net_parts <- function(p, layout) {
  face <- p[, "face"]
  xy <- unname(p[, c("tx", "ty"), drop = FALSE])
  # A wall along a face edge can leave a run on the face across it with no
  # area; dropping it joins the runs on either side.
  if (!all(face == face[1])) {
    run <- cumsum(c(TRUE, face[-1] != face[-length(face)]))
    flat <- vapply(split(seq_along(face), run), function(i) on_one_edge(xy[i, , drop = FALSE]),
                   logical(1))
    if (any(flat) && !all(flat)) {
      keep <- !flat[run]
      face <- face[keep]
      xy <- xy[keep, , drop = FALSE]
      if (any(face != face[1])) {
        face <- c(face, face[1])
        xy <- rbind(xy, xy[1, ])
      }
    }
  }
  if (all(face == face[1])) {
    polys <- list(xy[-nrow(xy), , drop = FALSE])
    faces <- face[1]
  } else {
    # The path's last row repeats its first; drop it, and start the walk
    # where the path crosses onto another face.
    n <- nrow(xy) - 1L
    face <- face[seq_len(n)]
    xy <- xy[seq_len(n), , drop = FALSE]
    start <- which(face != face[c(n, seq_len(n - 1L))])[1]
    ord <- c(start:n, seq_len(start - 1L))
    face <- face[ord]
    xy <- xy[ord, , drop = FALSE]
    run <- cumsum(c(TRUE, face[-1] != face[-n]))
    polys <- lapply(split(seq_len(n), run), function(i) close_on_face(xy[i, , drop = FALSE]))
    faces <- face[!duplicated(run)]
  }

  parts <- list()
  for (k in seq_along(polys)) {
    for (j in seq_along(layout$pieces)) {
      piece <- layout$pieces[[j]]
      if (piece$face != faces[k]) next
      if (nrow(polys[[k]]) < 3L) next
      part <- clip_polygon(polys[[k]], piece$region)
      if (nrow(part) < 3L) next
      placed <- place_points(part, piece)
      area <- polygon_area(placed)
      if (abs(area) < 1e-14) next
      cx <- sum((placed[, 1] + placed[c(2:nrow(placed), 1L), 1]) *
                  (placed[, 1] * placed[c(2:nrow(placed), 1L), 2] -
                     placed[c(2:nrow(placed), 1L), 1] * placed[, 2])) / (6 * area)
      cy <- sum((placed[, 2] + placed[c(2:nrow(placed), 1L), 2]) *
                  (placed[, 1] * placed[c(2:nrow(placed), 1L), 2] -
                     placed[c(2:nrow(placed), 1L), 1] * placed[, 2])) / (6 * area)
      parts[[length(parts) + 1L]] <- list(piece = j, area = abs(area), x = cx, y = cy,
                                          polygon = placed)
    }
  }
  if (length(parts) == 0L) return(NULL)
  list(piece = vapply(parts, `[[`, numeric(1), "piece"),
       area = vapply(parts, `[[`, numeric(1), "area"),
       x = vapply(parts, `[[`, numeric(1), "x"),
       y = vapply(parts, `[[`, numeric(1), "y"),
       polygon = lapply(parts, `[[`, "polygon"))
}

#' The pairs of pieces the map joins: two pieces that place the same edge of
#' the sphere at the same place. Returns a two-column matrix of piece
#' positions, each pair once.
#' @noRd
net_joins <- function(layout, tol = 1e-9) {
  e <- net_edges(layout, tol)
  e <- e[e$joined, , drop = FALSE]
  pairs <- do.call(rbind, lapply(seq_len(nrow(e)), function(i) {
    o <- which(e$key == e$key[i] & e$piece > e$piece[i] &
                 ((abs(e$x0 - e$x0[i]) < tol & abs(e$y0 - e$y0[i]) < tol &
                     abs(e$x1 - e$x1[i]) < tol & abs(e$y1 - e$y1[i]) < tol) |
                    (abs(e$x0 - e$x1[i]) < tol & abs(e$y0 - e$y1[i]) < tol &
                       abs(e$x1 - e$x0[i]) < tol & abs(e$y1 - e$y0[i]) < tol)))
    if (length(o) == 0L) return(NULL)
    cbind(e$piece[i], e$piece[o])
  }))
  if (is.null(pairs)) matrix(integer(0), 0L, 2L) else unique(pairs)
}

#' Groups of a cell's parts the map shows as one: two parts on joined pieces
#' are one where they share a stretch of the joining edge, that is, two
#' placed vertices
#' @noRd
join_groups <- function(parts, joins, tol = 1e-9) {
  n <- length(parts$piece)
  group <- seq_len(n)
  find <- function(i) {
    while (group[i] != i) i <- group[i]
    i
  }
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        a <- parts$piece[i]
        b <- parts$piece[j]
        if (!any((joins[, 1] == a & joins[, 2] == b) | (joins[, 1] == b & joins[, 2] == a))) next
        P <- parts$polygon[[i]]
        Q <- parts$polygon[[j]]
        shared <- sum(vapply(seq_len(nrow(P)), function(r) {
          any(abs(Q[, 1] - P[r, 1]) < tol & abs(Q[, 2] - P[r, 2]) < tol)
        }, logical(1)))
        if (shared >= 2L) group[find(i)] <- find(j)
      }
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  match(roots, unique(roots))
}
