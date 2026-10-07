# hex_sample.R
# Random points inside cells

#' Random points inside cells
#'
#' Draws points uniformly over the area of each cell on the sphere, as
#' DGGRID's \code{randpts} output does. Each point is assigned back to its
#' cell by \code{\link{lonlat_to_cell}}, so the points of a cell are exactly
#' the points the grid places in it. Use \code{set.seed()} for a reproducible
#' draw.
#'
#' Points are drawn uniformly in a spherical cap around the cell centre that
#' holds the whole cell, and those outside the cell are drawn again.
#'
#' @param cell_id Cell IDs. integer64 for ISEA grids, character for H3 grids.
#' @param grid A HexGridInfo or HexData object.
#' @param n Number of points per cell: one number, or one per cell.
#'
#' @return A data frame with columns \code{cell_id}, \code{lon} and
#'   \code{lat}, the points of each cell together in the order of
#'   \code{cell_id}.
#'
#' @seealso \code{\link{cell_to_lonlat}} for cell centres
#'
#' @export
#' @examples
#' grid <- hex_grid(resolution = 3, aperture = 3)
#' cell <- lonlat_to_cell(10, 50, grid)
#' set.seed(1)
#' pts <- hex_sample(cell, grid, n = 200)
#' plot(sf::st_geometry(cell_to_sf(cell, grid)))
#' points(pts$lon, pts$lat, pch = 20, cex = 0.5)
hex_sample <- function(cell_id, grid, n = 1L) {
  g <- extract_grid(grid)
  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  if (!is.numeric(n) || !(length(n) %in% c(1L, length(cell_id))) ||
      anyNA(n) || any(n < 0) || any(n != round(n))) {
    stop("n must be a non-negative whole number, or one per cell")
  }
  if (anyNA(cell_id)) stop("cell_id must not contain NA")
  n <- rep_len(as.integer(n), length(cell_id))

  ctr <- cell_to_lonlat(cell_id, g)
  C <- unit_vec(ctr$lon_deg, ctr$lat_deg)
  cos_cap <- cos(cell_cap_radius(cell_id, g, C))

  got <- list()
  need <- n
  while (any(need > 0L)) {
    k <- which(need > 0L)
    # A hexagon fills about 0.83 of its circumscribed cap, so this draw
    # usually finishes in one round.
    draw <- rep(k, ceiling(need[k] * 1.3) + 2L)
    P <- sample_cap(C[draw, , drop = FALSE], cos_cap[draw])
    ll <- vec_lonlat(P)
    lon <- ll[, 1]
    lat <- ll[, 2]
    hit <- which(lonlat_to_cell(lon, lat, g) == cell_id[draw])
    # draw is sorted, so the hits of each cell are consecutive
    rank <- sequence(rle(draw[hit])$lengths)
    hit <- hit[rank <= need[draw[hit]]]
    got[[length(got) + 1L]] <- cbind(draw[hit], lon[hit], lat[hit])
    need <- need - tabulate(draw[hit], nbins = length(need))
  }

  pts <- do.call(rbind, got)
  if (is.null(pts)) pts <- matrix(numeric(0), 0, 3)
  pts <- pts[order(pts[, 1]), , drop = FALSE]
  data.frame(cell_id = cell_id[pts[, 1]], lon = pts[, 2], lat = pts[, 3])
}

#' Angular radius of a cap around each cell centre that holds the cell
#'
#' The largest angle from the centre to the cell's boundary, traced to within
#' 0.001 of each edge's length, with 1 percent added for that tolerance.
#' @noRd
cell_cap_radius <- function(cell_id, g, C) {
  rings <- if (is_h3_grid(g)) {
    cpp_densify_great_circle(cpp_h3_cellToBoundary(as.character(cell_id)),
                             CELL_EDGE_TOLERANCE)
  } else {
    isea_cell_rings(cell_id, g@resolution, g@aperture, icosa_arg(g))
  }
  vapply(seq_along(rings), function(i) {
    R <- unit_vec(rings[[i]][, 1], rings[[i]][, 2])
    1.01 * acos(max(-1, min(1, min(R %*% C[i, ]))))
  }, numeric(1))
}

#' Points drawn uniformly in spherical caps
#'
#' One point per row of `C`, in the cap around the unit vector `C[i, ]` whose
#' angular radius has cosine `cos_cap[i]`.
#' @noRd
sample_cap <- function(C, cos_cap) {
  m <- nrow(C)
  z <- 1 - stats::runif(m) * (1 - cos_cap)
  r <- sqrt(pmax(0, 1 - z^2))
  phi <- stats::runif(m, 0, 2 * pi)
  # Two unit vectors perpendicular to each centre, from the axis it lies
  # furthest from
  axis <- matrix(0, m, 3)
  axis[cbind(seq_len(m), max.col(-abs(C), ties.method = "first"))] <- 1
  e1 <- cross3(C, axis)
  e1 <- e1 / sqrt(rowSums(e1^2))
  e2 <- cross3(C, e1)
  z * C + (r * cos(phi)) * e1 + (r * sin(phi)) * e2
}
