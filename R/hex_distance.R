# R/hex_distance.R
# Grid distance between hex cells

#' Grid Distance Between Cells
#'
#' Computes the grid distance (minimum number of hops) between pairs
#' of hexagonal cells. This is a discrete distance measured in cell steps,
#' not a geodesic distance.
#'
#' @param cell_a,cell_b Cell IDs. Must be the same length (pairs) or one
#'   of them length 1 (broadcast). For H3 grids, character vectors.
#'   For ISEA grids, numeric vectors.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#'
#' @return An integer vector of grid distances. `NA` where the distance
#'   cannot be computed (e.g., pentagon path issues in H3).
#'
#' @details
#' **H3 backend:** Uses the vendored H3 `gridDistance` function.
#' Returns the exact shortest path length in cell hops.
#'
#' **ISEA backend:** For cells on the same quad, computes the cube-coordinate
#' distance: `max(|di|, |dj|, |di + dj|)`. Cross-quad distances use BFS
#' expansion and may return `NA` for very distant cells.
#'
#' For geodesic (geographic) distances between cell centers, convert to
#' lon/lat with [cell_to_lonlat()] and use `sf::st_distance()`.
#'
#' @seealso [get_neighbors()] for finding cells within a given distance,
#'   [cell_to_lonlat()] for geographic coordinates
#'
#' @export
#' @examples
#' \donttest{
#' # H3 grid distance
#' g <- hex_grid(resolution = 5, type = "h3")
#' a <- lonlat_to_cell(10, 50, g)
#' b <- lonlat_to_cell(10.1, 50.1, g)
#' hex_distance(a, b, g)
#' }
hex_distance <- function(cell_a, cell_b, grid) {
  g <- extract_grid(grid)

  # Recycle to equal length
  n <- max(length(cell_a), length(cell_b))
  cell_a <- rep_len(cell_a, n)
  cell_b <- rep_len(cell_b, n)

  if (is_h3_grid(g)) {
    return(cpp_h3_gridDistance(as.character(cell_a), as.character(cell_b)))
  }

  # ISEA backend: same-quad cube distance
  ap <- as.integer(g@aperture)
  res <- g@resolution

  # Get quad IJ for both cells
  info_a <- hexify_cell_to_quad_ij(cell_a, res, ap)
  info_b <- hexify_cell_to_quad_ij(cell_b, res, ap)

  result <- integer(n)

  for (idx in seq_len(n)) {
    qa <- info_a$quad[idx]
    qb <- info_b$quad[idx]

    if (is.na(qa) || is.na(qb)) {
      result[idx] <- NA_integer_
      next
    }

    if (qa == qb) {
      # Same quad: cube distance
      di <- info_a$i[idx] - info_b$i[idx]
      dj <- info_a$j[idx] - info_b$j[idx]
      result[idx] <- as.integer(max(abs(di), abs(dj), abs(di + dj)))
    } else {
      # Cross-quad: BFS from cell_a, searching for cell_b
      # Limited search to avoid excessive computation
      max_search <- 100L
      found <- FALSE
      visited <- cell_a[idx]
      current <- cell_a[idx]
      target <- cell_b[idx]

      for (dist in seq_len(max_search)) {
        nbrs <- cpp_get_neighbors_isea(current, res, ap)
        new_cells <- unique(unlist(nbrs))
        new_cells <- setdiff(new_cells, visited)

        if (target %in% new_cells) {
          result[idx] <- dist
          found <- TRUE
          break
        }

        if (length(new_cells) == 0) break
        visited <- c(visited, new_cells)
        current <- new_cells
      }

      if (!found) {
        result[idx] <- NA_integer_
      }
    }
  }

  result
}
