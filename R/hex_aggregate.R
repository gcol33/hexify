# R/hex_aggregate.R
# Aggregation of cell values into the cells of a coarser resolution

#' Aggregate cell values to a coarser resolution
#'
#' Adds up, or averages, values given on cells into the cells \code{levels}
#' resolutions coarser. The hexagons of successive resolutions do not nest, so
#' a cell can lie partly in two or three coarser cells; \code{rule} says how
#' its value is divided among them.
#'
#' @param cell_id Cell IDs, each listed once.
#' @param value Numeric vector with one value per cell, or a numeric matrix
#'   or data frame with one row per cell, each column aggregated on its own.
#' @param grid A HexGridInfo or HexData object at the resolution of
#'   \code{cell_id}.
#' @param levels Number of resolutions up.
#' @param rule \code{"area"} divides each cell among the coarser cells it
#'   overlaps by the share of its area in each; \code{"centre"} gives each
#'   cell wholly to the coarser cell holding its centre, the one
#'   \code{\link{get_parent}} returns.
#' @param measure \code{"total"} adds the values up, for counts and other
#'   amounts; \code{"mean"} takes their mean weighted by the area each cell
#'   contributes, for densities and other values per unit area.
#' @param na.rm If \code{TRUE}, a cell whose value is \code{NA} counts as
#'   absent from that column.
#'
#' @return A data frame with one row per coarser cell that receives a share
#'   of any listed cell, sorted by \code{cell_id}, and the aggregated values:
#'   a column \code{value} for a vector, the columns of \code{value}
#'   otherwise. With \code{measure = "mean"}, a coarser cell whose
#'   contributing cells are all \code{NA} gets \code{NA}.
#'
#' @details
#' \strong{The two rules.} Write \eqn{|c \cap p| / |c|} for the share of
#' cell \eqn{c}'s area that lies in the coarser cell \eqn{p}. The area rule
#' gives \eqn{p} that share of \eqn{c}'s value, which assumes the amount is
#' spread evenly over \eqn{c}; the centre rule gives \eqn{p} all of it when
#' \eqn{p} holds \eqn{c}'s centre and nothing otherwise. A cell's shares add
#' up to one under both, so both keep the total of the values. Only the area
#' rule keeps where the amount lies: a uniform density gives every
#' hexagonal coarser cell the same total under it, and each vertex cell
#' (the twelve pentagons of the icosahedron) its sides over six of that.
#' A mean is the total of value times area over the total area, both
#' aggregated by the same rule.
#'
#' \strong{Shares on the lattice.} Within a face, the cells of one resolution
#' are the Voronoi cells of a hexagonal lattice and those of the next coarser
#' resolution the Voronoi cells of a sublattice, of index 3, 4 or 7, the
#' aperture of the step. After an aperture-3 step the cells are centred on
#' the parent centres and corners: the three parent edges meeting at a
#' corner run through every second edge midpoint of the cell centred there
#' and cut it into thirds. After an aperture-4 step the cells are centred on
#' the parent centres and edge midpoints: a parent edge runs through the
#' centre and two opposite corners of the cell on its midpoint and cuts it in
#' half. After an aperture-7 step every centre lies inside a parent. Take the
#' finer lattice as the Eisenstein integers \eqn{Z[\omega]},
#' \eqn{\omega = e^{2\pi i/3}}, with unit spacing, so a cell has circumradius
#' \eqn{s = 1/\sqrt{3}}. The parent lattice is \eqn{\mu Z[\omega]} with
#' \eqn{\mu = 3 + \omega}, \eqn{|\mu|^2 = 7}, turned by
#' \eqn{\arg \mu = \arctan(\sqrt{3}/5)} (the mirror image
#' \eqn{3 + \bar\omega} at the other twist). The centre child, of
#' circumradius \eqn{s}, lies inside its parent, of inradius
#' \eqn{\sqrt{7}/2}. The parent edge facing the neighbour at \eqn{\mu} runs
#' from the parent corner \eqn{\mu(1 - \omega)/3 = 1 + (1 - \omega)/3}, which
#' is a corner \eqn{A} of the ring child centred at 1, through the parent's
#' edge midpoint \eqn{\mu/2 = (1 + (2 + \omega))/2}, which is the midpoint of
#' the edge the children at 1 and \eqn{2 + \omega} share. That edge is the
#' second from \eqn{A}, so the part of the child at 1 beyond its parent is
#' the triangle of \eqn{A}, the next corner \eqn{B} and \eqn{\mu/2}: sides
#' \eqn{s} and \eqn{s/2} at 120 degrees, area \eqn{\sqrt{3}s^2/8} against
#' \eqn{3\sqrt{3}s^2/2} for the cell, one twelfth. The half-turn about
#' \eqn{\mu/2} swaps the two parents and the two children, so the edge cuts
#' one twelfth off the child at \eqn{2 + \omega} the other way. Each of the
#' six ring children lies 11/12 in its parent and 1/12 in one neighbour, and
#' the parent's area checks: \eqn{1 + 6 \cdot 11/12 + 6 \cdot 1/12 = 7}.
#'
#' A Hex9 step (aperture 9 on the octahedron, Griffin 2026) divides each
#' triangle of the solid's lattice into nine. A cell is the six lattice
#' triangles around its centre, so a parent's edges are lines of the coarser
#' lattice and therefore of the finer one. Six of a parent's nine children
#' lie inside it; the other three are centred on its edges, and the edge
#' through such a child's centre cuts it into its two half-hexagons, three
#' triangles each: half in each of two parents. The parent's area checks:
#' \eqn{6 + 6 \cdot 1/2 = 9} children's worth.
#'
#' The vertices of the solid are cell centres at every resolution, and the
#' turn by the angle deficit about a vertex maps both lattices onto
#' themselves, so the cells around a vertex are cut as in the plane; the cell
#' at a vertex lies wholly inside the vertex cell above it. On Hex9 the
#' vertices are lattice points no cell is centred on, and the cells there are
#' cut along lattice lines as everywhere. Snyder's projection
#' and the vertex-oriented one (\code{"ivea"}) are equal-area, so these plane
#' shares are the shares on the sphere, for every cell, across quad and face
#' edges and at the vertices. \code{\link{get_parent}} with
#' \code{overlapping = TRUE} finds the coarser cells a cell overlaps, its own
#' parent first, and the share follows from their number and the aperture of
#' the step: one parent, 1; three after an aperture-3 step, 1/3 each; two
#' after aperture 4, 1/2 each; two after aperture 7, 11/12 to the parent and
#' 1/12 to the neighbour; two after a Hex9 step, 1/2 each, the parent holding
#' the child's mode-0 half first.
#'
#' \strong{Several levels.} \code{levels > 1} applies the one-level rule level
#' by level, each step with the shares of its own aperture: a grid built with
#' \code{aperture = c(3, 7, 4)} aggregates from resolution 3 with halves, then
#' twelfths, then thirds. A value at an intermediate resolution is again taken
#' as spread evenly over its cell, so several levels up the area rule is the
#' one-level rule repeated rather than the share of a fine cell's area in a
#' coarse cell.
#'
#' \strong{Fuller's projection and H3.} Neither is equal-area, so the plane
#' shares do not carry over to the sphere, and each share is measured. Both
#' cells are clipped against each other where they are straight, which makes
#' the pieces exact: an ISEA-family cell on each face plane it covers, an H3
#' cell on the gnomonic plane at the finer cell's centre, where its
#' great-circle edges are straight. An H3 piece is a spherical polygon whose
#' area is a sum of spherical triangles, exact to rounding. A face-plane piece
#' has curved edges on the sphere; its area is that of the polygon through
#' equally spaced points of its edges, extrapolated to infinitely many points
#' (Romberg), which leaves an error near rounding. The share is a piece's area
#' over that of all the cell's pieces. Measured this way on Snyder's
#' projection, the shares come out as the lattice's.
#'
#' @references Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA
#'   discrete global grids. Statistical Computing & Graphics Newsletter
#'   8(2/3): 31-39.
#'
#' @seealso \code{\link{get_parent}} for the coarser cells a cell overlaps,
#'   \code{\link{hex_summarize}} to aggregate points into cells
#'
#' @export
#' @examples
#' g <- hex_grid(resolution = 4, aperture = 3)
#' ids <- grid_global(g)$cell_id
#' set.seed(1)
#' counts <- rpois(length(ids), 5)
#'
#' area <- hex_aggregate(ids, counts, g)
#' centre <- hex_aggregate(ids, counts, g, rule = "centre")
#' c(sum(counts), sum(area$value), sum(centre$value))
#'
#' # A uniform density: every parent gets its own area, the hexagons one
#' # value and the twelve pentagons 5/6 of it
#' density <- cell_area(ids, g)
#' range(hex_aggregate(ids, density, g, levels = 2)$value)
#'
#' # Values per unit area average instead
#' hex_aggregate(ids, data.frame(a = counts, b = 2 * counts), g,
#'               measure = "mean")[1:3, ]
hex_aggregate <- function(cell_id, value, grid, levels = 1L,
                          rule = c("area", "centre"),
                          measure = c("total", "mean"), na.rm = FALSE) {
  g <- extract_grid(grid)
  rule <- match.arg(rule)
  measure <- match.arg(measure)
  cell_id <- if (is_h3_grid(g)) as.character(cell_id) else as_cell_id(cell_id)
  if (anyDuplicated(cell_id)) stop("cell_id must list each cell once", call. = FALSE)
  if (length(levels) != 1L || is.na(levels) || levels < 1 || levels != round(levels) ||
      levels > g@resolution) {
    stop("levels must be a whole number from 1 to the grid's resolution", call. = FALSE)
  }
  if (length(na.rm) != 1L || is.na(na.rm) || !is.logical(na.rm)) {
    stop("na.rm must be TRUE or FALSE", call. = FALSE)
  }

  amount <- value_matrix(value, length(cell_id))
  missing <- is.na(amount)
  if (na.rm) amount[missing] <- 0
  covered <- NULL
  if (measure == "mean") {
    covered <- matrix(unname(cell_area(cell_id, g)), nrow(amount), ncol(amount))
    if (na.rm) covered[missing] <- 0
    amount <- amount * covered
  }

  ids <- cell_id
  for (step in seq_len(levels)) {
    w <- parent_weights(ids, grid_at_resolution(g, g@resolution - step + 1L), rule)
    ids <- unique(w$parent)
    slot <- match(w$parent, ids)
    amount <- weighted_rowsum(amount, w, slot)
    if (!is.null(covered)) covered <- weighted_rowsum(covered, w, slot)
  }

  if (!is.null(covered)) {
    amount <- amount / covered
    amount[covered == 0] <- NA_real_
  }
  ord <- order(ids)
  out <- data.frame(cell_id = ids[ord], stringsAsFactors = FALSE)
  out <- cbind(out, as.data.frame(amount[ord, , drop = FALSE]))
  rownames(out) <- NULL
  out
}

#' Values as a numeric matrix with one row per cell and named columns
#' @noRd
value_matrix <- function(value, n) {
  if (is.data.frame(value)) {
    if (!all(vapply(value, is.numeric, logical(1)))) {
      stop("value columns must be numeric", call. = FALSE)
    }
    value <- as.matrix(value)
  }
  if (!is.numeric(value)) stop("value must be numeric", call. = FALSE)
  if (is.null(dim(value))) {
    value <- matrix(as.numeric(value), ncol = 1L, dimnames = list(NULL, "value"))
  }
  if (nrow(value) != n) {
    stop("value must have one entry (or row) per cell", call. = FALSE)
  }
  if (is.null(colnames(value))) colnames(value) <- paste0("value", seq_len(ncol(value)))
  storage.mode(value) <- "double"
  value
}

#' Rows of a matrix, weighted, summed into the coarser cells they go to
#' @param x Matrix, one row per finer cell
#' @param w Weights from parent_weights()
#' @param slot Row of the result each weight goes to
#' @noRd
weighted_rowsum <- function(x, w, slot) {
  out <- rowsum(x[w$child, , drop = FALSE] * w$weight, slot, reorder = TRUE)
  rownames(out) <- NULL
  out
}

#' Shares of a cell in the coarser cells it overlaps after one step, on the
#' plane lattice: by the step's aperture, then by the number of coarser cells
#' the cell overlaps, its own parent first (see hex_aggregate()). NULL marks a
#' count the lattice does not produce.
#' @noRd
LATTICE_PARENT_SHARES <- list(
  `3` = list(1, NULL, rep(1 / 3, 3)),
  `4` = list(1, c(1 / 2, 1 / 2)),
  `7` = list(1, c(11 / 12, 1 / 12)),
  `9` = list(1, c(1 / 2, 1 / 2))
)

#' How far the pieces of a cell may fall short of or pass its measured area
#' before the coarser cells found for it count as incomplete
#' @noRd
SHARE_COVER_SLACK <- 1e-6

#' The coarser cells each cell goes to one level up, with the share of the
#' cell each gets
#'
#' @param cell_id Cell IDs on `g`
#' @param g HexGridInfo object, the cells' grid
#' @param rule "area" or "centre"
#' @return List of `child` (position in `cell_id`), `parent` (coarser cell ID)
#'   and `weight`, one entry per pair, a cell's shares adding up to one
#' @noRd
parent_weights <- function(cell_id, g, rule) {
  if (rule == "centre") {
    return(list(child = seq_along(cell_id), parent = get_parent(cell_id, g),
                weight = rep(1, length(cell_id))))
  }
  ov <- get_parent(cell_id, g, overlapping = TRUE)
  parent <- if (is_h3_grid(g)) unlist(ov, use.names = FALSE) else cell_id_unlist(ov)
  weight <- if (!is_h3_grid(g) && is_equal_area_projection(grid_projection(g))) {
    lattice_shares(lengths(ov), step_aperture(g))
  } else {
    sphere_shares(cell_id, ov, g)
  }
  list(child = rep(seq_along(cell_id), lengths(ov)), parent = parent, weight = weight)
}

#' Plane-lattice shares, concatenated over cells, from the number of coarser
#' cells each overlaps
#' @noRd
lattice_shares <- function(n_parents, aperture) {
  shares <- LATTICE_PARENT_SHARES[[as.character(aperture)]]
  known <- n_parents >= 1L & n_parents <= length(shares)
  known[known] <- !vapply(shares[n_parents[known]], is.null, logical(1))
  if (!all(known)) {
    stop(sprintf(paste0("hexify internal error: a cell overlaps %d coarser cells ",
                        "after an aperture-%d step"),
                 n_parents[!known][1], aperture), call. = FALSE)
  }
  unlist(shares[n_parents], use.names = FALSE)
}

#' Aperture of the step from one resolution coarser to a grid's resolution
#'
#' A family spelling can take the step out of the middle of its sequence
#' (ISEA43H at resolution 3 is 4,3,3 and at 4 is 4,4,3,3), so the step is the
#' aperture the finer sequence has once more than the coarser one.
#' @noRd
step_aperture <- function(g) {
  if (!is_mixed_aperture(g@aperture)) return(as.integer(g@aperture))
  counts <- function(r) tabulate(isea_levels(g@aperture, r)$ap_seq[-1L], 7L)
  d <- counts(g@resolution) - counts(g@resolution - 1L)
  if (sum(d == 1L) != 1L || any(d != 0L & d != 1L)) {
    stop("hexify internal error: the coarser aperture sequence is not the ",
         "finer one with one step removed", call. = FALSE)
  }
  which(d == 1L)
}

#' Shares measured on the sphere: each cell clipped against every coarser
#' cell it overlaps where both are straight, and the pieces' areas measured
#'
#' ISEA-family cells are straight on each face (cpp_cell_overlap_solid_angles()),
#' H3 cells on the gnomonic plane, whose straight lines are great circles
#' (cpp_h3_overlap_solid_angles()).
#'
#' @param cell_id Cell IDs on `g`
#' @param ov List of the coarser cells each overlaps (get_parent(overlapping))
#' @param g HexGridInfo object
#' @return Numeric vector of shares, concatenated over cells as `ov` is
#' @noRd
sphere_shares <- function(cell_id, ov, g) {
  pg <- grid_at_resolution(g, g@resolution - 1L)
  flat <- if (is_h3_grid(g)) unlist(ov, use.names = FALSE) else cell_id_unlist(ov)
  parents <- unique(flat)
  pair_cell <- rep(seq_along(cell_id), lengths(ov))
  pair_parent <- match(flat, parents)
  omega <- if (is_h3_grid(g)) {
    cpp_h3_overlap_solid_angles(as.character(cell_id), as.character(parents),
                                pair_cell, pair_parent)
  } else {
    lv <- isea_levels(g@aperture, g@resolution)
    pl <- isea_levels(pg@aperture, pg@resolution)
    cpp_cell_overlap_solid_angles(icosa_arg(g), as_cell_id(cell_id), lv$resolution,
                                  lv$aperture, lv$ap_seq, parents, pl$resolution,
                                  pl$aperture, pl$ap_seq, pair_cell, pair_parent)
  }
  total <- as.numeric(rowsum(omega$piece, pair_cell, reorder = TRUE))
  cover <- total / omega$whole
  if (any(abs(cover - 1) > SHARE_COVER_SLACK)) {
    k <- which.max(abs(cover - 1))
    stop(sprintf(paste0("hexify internal error: the coarser cells found for cell %s ",
                        "cover %.12f of it"), as.character(cell_id[k]), cover[k]),
         call. = FALSE)
  }
  omega$piece / total[pair_cell]
}
