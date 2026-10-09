# hex_glyphs.R
# Glyphs drawn inside cells: ray glyphs and the triangles towards each
# neighbour (Carr et al. 1997)

#' Ray glyphs for cell values
#'
#' Draws one ray per cell from its centre whose angle shows a value, as in
#' Carr et al. (1997): the ray points down for the smallest value, level
#' for the middle of the range and up for the largest. A confidence arc at
#' the ray's tip spans the angles of a lower and an upper bound. A second
#' variable gets a ray of its own pointing to the left from the same centre,
#' so one glyph shows two values.
#'
#' Rays run in the plane tangent to the sphere at the cell centre, with up
#' towards north, and are returned as lines in longitude and latitude, ready
#' for any map of the cells.
#'
#' @param cell_id Cell IDs, each listed once.
#' @param value Values shown by the right-hand rays.
#' @param grid A HexGridInfo or HexData object.
#' @param lower,upper Optional bounds of \code{value} (such as a confidence
#'   interval), drawn as an arc at the tip of the ray.
#' @param value2,lower2,upper2 Optional second variable and its bounds,
#'   drawn to the left.
#' @param limits Values drawn straight down and straight up: \code{c(min,
#'   max)}. \code{NULL} uses the range of the values and bounds; with a
#'   second variable both share it unless \code{limits2} is given.
#' @param limits2 Limits of the second variable.
#' @param length Length of a ray as a fraction of the distance from the
#'   cell centre to the middle of a wall.
#'
#' @return An sf object of LINESTRINGs with columns \code{cell_id},
#'   \code{side} (\code{"right"} or \code{"left"}) and \code{part}
#'   (\code{"ray"} or \code{"arc"}).
#'
#' @references Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA
#'   discrete global grids. Statistical Computing & Graphics Newsletter
#'   8(2/3): 31-39.
#'
#' @seealso \code{\link{hex_triangles}} for the change towards each
#'   neighbour
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 4, aperture = 3)
#' cells <- lonlat_to_cell(c(5, 10, 15), c(45, 47, 49), g)
#' trend <- c(-1, 0.2, 1.4)
#' rays <- hex_rays(cells, trend, g, lower = trend - 0.5, upper = trend + 0.5)
#' plot(sf::st_geometry(cell_to_sf(cells, g)))
#' plot(sf::st_geometry(rays), add = TRUE, col = "#D55E00")
hex_rays <- function(cell_id, value, grid, lower = NULL, upper = NULL,
                     value2 = NULL, lower2 = NULL, upper2 = NULL,
                     limits = NULL, limits2 = NULL, length = 0.85) {
  g <- extract_grid(grid)
  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  n <- length(cell_id)
  check_glyph_values(n, value = value, lower = lower, upper = upper,
                     value2 = value2, lower2 = lower2, upper2 = upper2)
  if ((is.null(lower) != is.null(upper)) || (is.null(lower2) != is.null(upper2))) {
    stop("give both bounds, lower and upper, or neither", call. = FALSE)
  }
  if (!is.null(lower2) && is.null(value2)) {
    stop("lower2 and upper2 bound value2", call. = FALSE)
  }
  if (!is.numeric(length) || base::length(length) != 1L || !is.finite(length) ||
      length <= 0) {
    stop("length must be a positive number", call. = FALSE)
  }
  limits <- glyph_limits(limits, c(value, lower, upper, if (is.null(limits2)) c(value2, lower2, upper2)))
  limits2 <- if (is.null(value2)) NULL else
    glyph_limits(if (is.null(limits2)) limits else limits2, c(value2, lower2, upper2))

  ctr <- cell_to_lonlat(cell_id, g)
  icosa <- icosa_arg(g)
  reach <- length * cell_inradius(cell_id, g)
  sides <- list(list(side = "right", sign = 1, value = value, lower = lower,
                     upper = upper, limits = limits))
  if (!is.null(value2)) {
    sides[[2]] <- list(side = "left", sign = -1, value = value2, lower = lower2,
                       upper = upper2, limits = limits2)
  }

  geoms <- list()
  meta <- list()
  for (s in sides) {
    for (i in seq_len(n)) {
      if (is.na(s$value[i])) next
      a <- glyph_angle(s$value[i], s$limits, s$sign)
      ray <- tangent_points(ctr$lon_deg[i], ctr$lat_deg[i], c(0, reach[i]), c(a, a),
                            icosa)
      geoms[[base::length(geoms) + 1L]] <- sf::st_linestring(ray)
      meta[[base::length(meta) + 1L]] <- c(i, s$side, "ray")
      if (!is.null(s$lower) && !is.na(s$lower[i]) && !is.na(s$upper[i])) {
        span <- glyph_angle(c(s$lower[i], s$upper[i]), s$limits, s$sign)
        ang <- seq(span[1], span[2], length.out = 25L)
        arc <- tangent_points(ctr$lon_deg[i], ctr$lat_deg[i], rep(reach[i], 25L), ang,
                              icosa)
        geoms[[base::length(geoms) + 1L]] <- sf::st_linestring(arc)
        meta[[base::length(meta) + 1L]] <- c(i, s$side, "arc")
      }
    }
  }
  meta <- do.call(rbind, meta)
  out <- sf::st_sf(cell_id = cell_id[as.integer(meta[, 1])], side = meta[, 2],
                   part = meta[, 3], geometry = sf::st_sfc(geoms, crs = grid_crs(g)))
  wrap_cells_at_dateline(out)
}

#' The triangles of cells towards each neighbour
#'
#' Cuts each cell into one triangle per wall, from the centre to that wall:
#' six for a hexagon, five for a pentagon. With values, each triangle carries
#' the change from its cell's value to the neighbour's across the wall, so
#' colouring the triangles shows the change from cell to cell in every
#' direction (Carr et al. 1997, after K. and R. Keister).
#'
#' @param cell_id Cell IDs, each listed once.
#' @param grid A HexGridInfo or HexData object.
#' @param value Optional values, one per cell. A neighbour not in
#'   \code{cell_id} or with an \code{NA} value gives an \code{NA} change.
#'
#' @return An sf object of POLYGONs with columns \code{cell_id},
#'   \code{neighbor_id} and, with values, \code{value}, \code{neighbor_value}
#'   and \code{change} (\code{neighbor_value - value}).
#'
#' @references Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA
#'   discrete global grids. Statistical Computing & Graphics Newsletter
#'   8(2/3): 31-39.
#'
#' @seealso \code{\link{hex_rays}}, \code{\link{get_neighbors}},
#'   \code{\link{wall_metrics}}
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' cells <- grid_rect(c(0, 40, 30, 60), g)$cell_id
#' ctr <- cell_to_lonlat(cells, g)
#' tri <- hex_triangles(cells, g, value = ctr$lat_deg)
#' plot(tri["change"], border = NA)
hex_triangles <- function(cell_id, grid, value = NULL) {
  g <- extract_grid(grid)
  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  if (anyDuplicated(cell_id)) stop("cell_id must list each cell once", call. = FALSE)
  if (!is.null(value)) {
    check_glyph_values(length(cell_id), value = value)
  }
  rings <- cell_rings_lonlat(cell_id, g)
  ctr <- cell_to_lonlat(cell_id, g)

  geoms <- list()
  from <- integer(0)
  nb <- list()
  for (i in seq_along(cell_id)) {
    walls <- ring_walls(rings[[i]], c(ctr$lon_deg[i], ctr$lat_deg[i]), g)
    for (w in walls) {
      geoms[[length(geoms) + 1L]] <- sf::st_polygon(list(
        triangle_ring(c(ctr$lon_deg[i], ctr$lat_deg[i]), w$points)))
      from <- c(from, i)
      nb[[length(nb) + 1L]] <- w$neighbor
    }
  }
  nb <- if (is_h3_grid(g)) as.character(unlist(nb)) else cell_id_unlist(nb)
  out <- data.frame(cell_id = cell_id[from], neighbor_id = nb)
  if (!is.null(value)) {
    out$value <- value[from]
    out$neighbor_value <- value[match(out$neighbor_id, cell_id)]
    out$change <- out$neighbor_value - out$value
  }
  out <- sf::st_sf(out, geometry = sf::st_sfc(geoms, crs = grid_crs(g)))
  wrap_cells_at_dateline(out)
}

# =============================================================================
# HELPERS
# =============================================================================

#' Check glyph value vectors: numeric, one per cell (NULL skipped)
#' @noRd
check_glyph_values <- function(n, ...) {
  args <- list(...)
  for (nm in names(args)) {
    x <- args[[nm]]
    if (is.null(x)) next
    if (!is.numeric(x) || length(x) != n) {
      stop(nm, " must be numeric with one value per cell", call. = FALSE)
    }
  }
  invisible()
}

#' Limits of a ray glyph: given, or the range of the values
#' @noRd
glyph_limits <- function(limits, values) {
  if (is.null(limits)) limits <- range(values, na.rm = TRUE)
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[1] >= limits[2]) {
    stop("limits must be c(min, max) with min < max", call. = FALSE)
  }
  limits
}

#' Ray angle in degrees anticlockwise from east: -90 (down) at the lower
#' limit, 90 (up) at the upper, mirrored to the left for `sign = -1`
#' @noRd
glyph_angle <- function(v, limits, sign) {
  t <- pmin(pmax((v - limits[1]) / (limits[2] - limits[1]), 0), 1)
  a <- -90 + 180 * t
  if (sign > 0) a else 180 - a
}

#' Points at distances `d` (radians) from (lon, lat) in directions `angle`
#' (degrees anticlockwise from east), along great circles of the grid's
#' sphere (latitudes geodetic on its ellipsoid, if `icosa` carries one);
#' longitudes kept continuous with the centre's
#' @noRd
tangent_points <- function(lon, lat, d, angle, icosa = numeric(0)) {
  bearing <- (90 - angle) * pi / 180
  la <- sphere_lat(lat, icosa) * pi / 180
  lat2 <- asin(sin(la) * cos(d) + cos(la) * sin(d) * cos(bearing))
  dlon <- atan2(sin(bearing) * sin(d) * cos(la), cos(d) - sin(la) * sin(lat2))
  cbind(lon + dlon * 180 / pi, geodetic_lat(lat2 * 180 / pi, icosa))
}

#' Angular distance (radians) from each cell's centre to the middle of a
#' wall: the inradius of the regular hexagon with the cell's area
#' @noRd
cell_inradius <- function(cell_id, g) {
  omega <- cell_area(cell_id, g) / body_surface_km2(grid_radius_km(g)) * 4 * pi
  hexagon <- !is_pentagon(cell_id, g)
  r <- sqrt(omega / (2 * sqrt(3)))
  # A pentagon has five sixths of a hexagon's area at the same inradius.
  r[!hexagon] <- sqrt(omega[!hexagon] * 6 / 5 / (2 * sqrt(3)))
  unname(r)
}

#' Cell boundaries as closed lon/lat rings, densified as cell_to_sf() draws
#' them
#' @noRd
cell_rings_lonlat <- function(cell_id, g) {
  if (is_h3_grid(g)) {
    return(cpp_densify_great_circle(cpp_h3_cellToBoundary(as.character(cell_id)),
                                    CELL_EDGE_TOLERANCE))
  }
  isea_cell_rings(cell_id, g@resolution, g@aperture, icosa_arg(g))
}

#' The walls of a cell: its ring cut into runs of points that border the same
#' neighbour. A ring segment borders the cell found a little outside its
#' midpoint, away from the centre. Returns one list per wall with the
#' neighbour's ID and the wall's points, in ring order.
#' @noRd
ring_walls <- function(ring, centre, g) {
  ring <- ring[c(TRUE, rowSums(abs(diff(ring[, 1:2, drop = FALSE]))) > 0), , drop = FALSE]
  n <- nrow(ring) - 1L
  icosa <- icosa_arg(g)
  P <- unit_vec(ring[seq_len(n + 1L), 1], ring[seq_len(n + 1L), 2], icosa)
  C <- drop(unit_vec(centre[1], centre[2], icosa))
  M <- P[seq_len(n), , drop = FALSE] + P[seq_len(n) + 1L, , drop = FALSE]
  M <- M / sqrt(rowSums(M^2))
  probe <- M + 0.02 * sweep(M, 2, C)
  probe <- vec_lonlat(probe / sqrt(rowSums(probe^2)), icosa)
  nb <- lonlat_to_cell(probe[, 1], probe[, 2], g)
  # Start at a change of neighbour so no wall wraps around the ring's start
  start <- which(nb != nb[c(n, seq_len(n - 1L))])[1]
  if (is.na(start)) start <- 1L
  ord <- c(start:n, seq_len(start - 1L))
  run <- cumsum(c(TRUE, nb[ord][-1] != nb[ord][-n]))
  lapply(split(ord, run), function(seg) {
    pts <- ring[c(seg, seg[length(seg)] + 1L), 1:2, drop = FALSE]
    list(neighbor = nb[seg[1]], points = pts)
  })
}

#' The closed lon/lat ring of the triangle from a cell centre to one wall. A
#' centre at a pole has no longitude, so it stands at each end of the wall's
#' longitudes.
#' @noRd
triangle_ring <- function(centre, wall) {
  ring <- if (abs(abs(centre[2]) - 90) < 1e-9) {
    rbind(c(wall[1, 1], centre[2]), wall, c(wall[nrow(wall), 1], centre[2]))
  } else {
    rbind(centre, wall)
  }
  ring <- unname(rbind(ring, ring[1, ]))
  lonlat_ring_coords(ring)
}
