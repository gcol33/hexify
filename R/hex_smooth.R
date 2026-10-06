# R/hex_smooth.R
# Neighbour smoothing of cell values

#' Smooth cell values over neighbouring cells
#'
#' Replaces each cell's value by the weighted mean of its own value and the
#' values of its edge neighbours, and repeats this \code{steps} times: the
#' neighbour-matrix smoother of Carr et al. (1997). Each neighbour has weight
#' 1 and the cell itself \code{self_weight}. On an equal-area grid every cell
#' has the same area, so the mean is also the area-weighted mean.
#'
#' Only the cells in \code{cell_id} take part: a neighbour that is not listed,
#' or whose value is \code{NA}, does not contribute, so a cell at the edge of
#' the data averages over the neighbours it has. A pentagon averages over its
#' five neighbours. A cell whose value is \code{NA} stays \code{NA}.
#'
#' \code{group} restricts the smoothing to cells of the same group, such as
#' land and ocean: a neighbour in another group does not contribute, so values
#' do not mix across a coastline. A cell whose group is \code{NA} keeps its
#' value and contributes to no other cell.
#'
#' @param cell_id Cell IDs, each listed once.
#' @param value Numeric values, one per cell.
#' @param grid A HexGridInfo or HexData object.
#' @param steps Number of smoothing steps. Each step widens the reach by one
#'   ring of cells.
#' @param self_weight Weight of a cell's own value relative to each
#'   neighbour's. \code{0} replaces a value by the mean of its neighbours.
#' @param group Optional vector, one entry per cell; only neighbours in the
#'   same group contribute.
#'
#' @return Numeric vector of smoothed values, in the order of \code{cell_id}.
#'
#' @references Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA
#'   discrete global grids. Statistical Computing & Graphics Newsletter
#'   8(2/3): 31-39.
#'
#' @seealso [get_neighbors()] for the neighbours a step averages over
#'
#' @export
#' @examples
#' \donttest{
#' g <- hex_grid(resolution = 4, aperture = 3)
#' cells <- grid_global(g)
#' set.seed(1)
#' noise <- rnorm(nrow(cells))
#' smooth <- hex_smooth(cells$cell_id, noise, g, steps = 3)
#' c(sd(noise), sd(smooth))
#'
#' # Keep land and ocean apart
#' ctr <- sf::st_as_sf(cell_to_lonlat(cells$cell_id, g),
#'                     coords = c("lon_deg", "lat_deg"), crs = 4326)
#' land <- lengths(sf::st_intersects(ctr, hexify_world)) > 0
#' coastal <- hex_smooth(cells$cell_id, noise, g, steps = 3, group = land)
#' }
hex_smooth <- function(cell_id, value, grid, steps = 1L, self_weight = 1,
                       group = NULL) {
  g <- extract_grid(grid)
  n <- length(cell_id)
  validate_same_length(list(cell_id = cell_id, value = value), "hex_smooth")
  if (!is.numeric(value)) stop("value must be numeric", call. = FALSE)
  if (anyDuplicated(cell_id)) stop("cell_id must list each cell once", call. = FALSE)
  if (length(steps) != 1L || is.na(steps) || steps < 0 || steps != round(steps)) {
    stop("steps must be a non-negative whole number", call. = FALSE)
  }
  if (length(self_weight) != 1L || !is.finite(self_weight) || self_weight < 0) {
    stop("self_weight must be a non-negative number", call. = FALSE)
  }
  if (!is.null(group) && length(group) != n) {
    stop("group must have one entry per cell", call. = FALSE)
  }
  value <- as.numeric(value)
  if (n == 0L || steps == 0L) return(value)

  # Directed links i <- j from each cell to its listed neighbours
  nbrs <- get_neighbors(cell_id, g)
  i <- rep.int(seq_len(n), lengths(nbrs))
  j <- match(unlist(nbrs, use.names = FALSE), cell_id)
  keep <- !is.na(j)
  if (!is.null(group)) {
    keep <- keep & !is.na(group[i]) & !is.na(group[j]) & group[i] == group[j]
  }
  i <- i[keep]
  j <- j[keep]

  fixed <- is.na(value)
  if (!is.null(group)) fixed <- fixed | is.na(group)
  link_ok <- !fixed[j] & !fixed[i]
  i <- i[link_ok]
  j <- j[link_ok]

  for (step in seq_len(steps)) {
    total <- self_weight * value + link_sum(value[j], i, n)
    weight <- self_weight + link_sum(rep(1, length(j)), i, n)
    smoothed <- total / weight
    # A cell with no contributing neighbour and self weight 0 keeps its value
    smoothed[weight == 0] <- value[weight == 0]
    value[!fixed] <- smoothed[!fixed]
  }
  value
}

#' Sum of link values per target cell
#' @noRd
link_sum <- function(x, i, n) {
  out <- numeric(n)
  if (length(i)) {
    s <- rowsum(x, i, reorder = FALSE)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}
