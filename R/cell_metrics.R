# Shape and spacing metrics of grid cells

#' Shape metrics of cells: area, perimeter and compactness
#'
#' Measures each cell on the sphere: its area, the length of its boundary and
#' how close its shape is to a circle.
#'
#' @param cell_id Cell IDs. When \code{NULL}, the cells of a HexData
#'   \code{grid}, or every cell of a HexGridInfo \code{grid}.
#' @param grid A HexGridInfo or HexData object.
#'
#' @return Data frame with one row per \code{cell_id}: \code{cell_id},
#'   \code{area_km2} (as \code{\link{cell_area}} gives it),
#'   \code{normalized_area}, \code{perimeter_km}, \code{compactness} and
#'   \code{ipq}.
#'
#' @details
#' The perimeter follows each wall as the grid draws it. An ISEA wall is
#' straight on its face of the solid and curved on the sphere, so it is
#' followed on the sphere until its length is exact to about 1e-8. An H3
#' wall is a great-circle arc between corners, with an extra corner where it
#' crosses an edge of the icosahedron.
#'
#' Compactness (White et al. 1998) is the perimeter of a spherical cap with the
#' cell's area over the cell's perimeter, \eqn{\sqrt{4\pi a - a^2/r^2} / p} for
#' area \eqn{a}, perimeter \eqn{p} and radius \eqn{r}. A circle scores 1; a
#' small regular hexagon 0.952, pentagon 0.930, square 0.886 and triangle
#' 0.778.
#'
#' The normalized area is the cell's area over the mean area of the grid's
#' cells, the body's surface over the cell count; an equal-area grid of
#' \eqn{N} cells gives \eqn{N / (N - 2)} for its hexagons, 5/6 of that for
#' the pentagons of the icosahedron, 4/6 for the squares of the octahedron and
#' 3/6 for the triangles of the tetrahedron. The isoperimetric quotient
#' \code{ipq} is \eqn{4\pi a / p^2}, for a circle 1, a small regular hexagon
#' 0.907, square 0.785 and triangle 0.605. Kmoch et al. (2022) compare
#' grids by these two measures, reading area and perimeter in an equal-area
#' plane centred on each cell; here both are read on the sphere, which adds
#' \eqn{(a/(rp))^2} to \code{compactness}^2, a few parts in a million for
#' cells of a few thousand square kilometres.
#'
#' @references
#' White, D., Kimerling, A. J., Sahr, K. and Song, L. (1998). Comparing area
#' and shape distortion on polyhedral-based recursive partitions of the sphere.
#' International Journal of Geographical Information Science 12(8), 805-827.
#' \doi{10.1080/136588198241518}
#'
#' Kmoch, A., Vasilyev, I., Virro, H. and Uuemaa, E. (2022). Area and shape
#' distortions in open-source discrete global grid systems. Big Earth Data
#' 6(3), 256-275. \doi{10.1080/20964471.2022.2094926}
#'
#' @seealso \code{\link{wall_metrics}} for the spacing of neighbouring cells,
#'   \code{\link{cell_area}}
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' m <- cell_metrics(grid = g)
#' summary(m$compactness)
#'
#' # Fuller's projection is not equal-area: compare the spread of areas
#' gf <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
#' mf <- cell_metrics(grid = gf)
#' sd(mf$area_km2) / mean(mf$area_km2)
cell_metrics <- function(cell_id = NULL, grid) {
  resolved <- resolve_cells_grid(cell_id, grid, all = TRUE)
  g <- resolved$grid
  cell_id <- metric_ids(resolved$cell_id, g)
  ids <- unique(cell_id[!is.na(cell_id)])

  radius <- grid_radius_km(g)
  area <- unname(cell_area(ids, g))
  perimeter <- cell_wall_measures(ids, g, walls = FALSE)$perimeter
  # Solid angle from the area keeps compactness free of the body's shape:
  # an Earth grid's areas are read on the ellipsoid, its lengths on the sphere.
  omega <- area / body_surface_km2(radius) * 4 * pi

  at <- match(cell_id, ids)
  data.frame(cell_id = cell_id,
             area_km2 = area[at],
             normalized_area = area[at] / body_surface_km2(radius) * grid_n_cells(g),
             perimeter_km = perimeter[at] * radius,
             compactness = (sqrt(4 * pi * omega - omega^2) / perimeter)[at],
             ipq = (4 * pi * omega / perimeter^2)[at],
             stringsAsFactors = FALSE)
}

#' Spacing metrics of neighbouring cells: wall length, centre distance and
#' cell wall midpoint ratio
#'
#' Measures every wall between a cell and an adjacent cell: how long it is,
#' how far apart the two cell centres are, and how far the wall's midpoint lies
#' from the line between the centres.
#'
#' @inheritParams cell_metrics
#'
#' @return Data frame with one row per wall: \code{cell_id},
#'   \code{neighbor_id}, \code{wall_km}, \code{center_distance_km} and
#'   \code{midpoint_ratio}. A wall between two of the given cells appears once,
#'   on the row of the cell given first; a wall to a cell not given appears on
#'   the given cell's row.
#'
#' @details
#' The centre distance is the great-circle distance between the two cell
#' centres. Gregory et al. (2008) measure the spread of a grid's intercell
#' distances by their coefficient of variation.
#'
#' The cell wall midpoint ratio (Gregory et al. 2008, after Heikes and Randall
#' 1995) is the distance from the wall's midpoint to the midpoint of the
#' great-circle arc joining the two centres, over the wall's length. It is 0
#' where the arc between the centres crosses the wall at its middle, as it does
#' between two regular hexagons; a finite-volume scheme that evaluates a flux
#' at the wall's midpoint is most accurate there. The wall's midpoint is the
#' point halfway along the wall as the grid draws it (see
#' \code{\link{cell_metrics}}).
#'
#' @references
#' Gregory, M. J., Kimerling, A. J., White, D. and Sahr, K. (2008). A
#' comparison of intercell metrics on discrete global grid systems. Computers,
#' Environment and Urban Systems 32(3), 188-203.
#' \doi{10.1016/j.compenvurbsys.2007.11.003}
#'
#' Heikes, R. and Randall, D. A. (1995). Numerical integration of the
#' shallow-water equations on a twisted icosahedral grid. Part II. A detailed
#' description of the grid and an analysis of numerical accuracy. Monthly
#' Weather Review 123(6), 1881-1887.
#' \doi{10.1175/1520-0493(1995)123<1881:NIOTSW>2.0.CO;2}
#'
#' @seealso \code{\link{cell_metrics}} for cell shape,
#'   \code{\link{get_neighbors}}
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 3, aperture = 3)
#' w <- wall_metrics(grid = g)
#' sd(w$center_distance_km) / mean(w$center_distance_km)
#' summary(w$midpoint_ratio)
#'
#' # The walls around one cell
#' cell <- lonlat_to_cell(16.37, 48.21, g)
#' wall_metrics(cell, g)
wall_metrics <- function(cell_id = NULL, grid) {
  resolved <- resolve_cells_grid(cell_id, grid, all = TRUE)
  g <- resolved$grid
  cell_id <- metric_ids(resolved$cell_id, g)
  ids <- unique(cell_id[!is.na(cell_id)])

  w <- cell_wall_measures(ids, g, walls = TRUE)$walls
  from <- w$cell
  to <- match(w$neighbor_id, ids)
  keep <- is.na(to) | from < to
  w <- w[keep, , drop = FALSE]

  radius <- grid_radius_km(g)
  data.frame(cell_id = ids[w$cell],
             neighbor_id = w$neighbor_id,
             wall_km = w$wall * radius,
             center_distance_km = w$centre_distance * radius,
             midpoint_ratio = w$midpoint_offset / w$wall,
             stringsAsFactors = FALSE)
}

#' Cell IDs in the type a grid's backend takes
#' @noRd
metric_ids <- function(cell_id, g) {
  if (is_h3_grid(g)) as.character(cell_id) else as_cell_id(cell_id)
}

#' Perimeters and walls of cells on the unit sphere
#'
#' @param ids Unique cell IDs, without NA
#' @param g HexGridInfo object
#' @param walls Whether to measure each wall against its neighbour
#' @return List with `perimeter` (radians, one per id) and, with `walls`, a
#'   data frame `walls` of `cell` (position in `ids`), `neighbor_id`, `wall`,
#'   `centre_distance` and `midpoint_offset`, all in radians
#' @noRd
cell_wall_measures <- function(ids, g, walls) {
  if (is_h3_grid(g)) return(cpp_h3_cell_walls(ids, walls))
  lv <- isea_levels(g@aperture, g@resolution)
  cpp_cell_walls(icosa_arg(g), ids, lv$resolution, lv$aperture, lv$ap_seq,
                 CELL_WALL_TOLERANCE, walls)
}
