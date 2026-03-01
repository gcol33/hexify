# R/hexify_neighbors.R
# Neighbor finding for ISEA and H3 grids
#
# get_neighbors() returns the k-ring (disk) of neighboring cells
# around a given cell ID.

#' Get Neighboring Cells
#'
#' Returns the k-ring (disk) of cells neighboring the input cells.
#' For `k = 1`, returns the immediate 6 neighbors (5 for pentagons).
#' For `k > 1`, returns all cells within `k` grid hops.
#'
#' @param cell_id Cell IDs to find neighbors for. Numeric vector for ISEA
#'   grids, character vector for H3 grids.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#' @param k Integer. Ring distance (default 1). `k = 1` returns immediate
#'   neighbors, `k = 2` includes neighbors-of-neighbors, etc.
#' @param include_self Logical. If `TRUE`, include the input cell in the
#'   result (default `FALSE`).
#' @param distances Logical. If `TRUE`, return a data.frame with cell IDs
#'   and their ring distance from the origin (default `FALSE`).
#'
#' @return If `distances = FALSE` (default): a list of cell ID vectors, one
#'   per input cell. If `distances = TRUE`: a list of data.frames with columns
#'   `cell_id` and `ring_distance`.
#'
#' @details
#' For **ISEA grids**, neighbors are computed using axial coordinate offsets
#' in the quad IJ space. Cells at quad boundaries are handled by reprojection
#' through lon/lat coordinates.
#'
#' For **H3 grids**, neighbors use the vendored H3 `gridDisk` /
#' `gridDiskDistances` functions.
#'
#' Pentagon cells (the 12 icosahedron vertices) have only 5 neighbors instead
#' of the usual 6.
#'
#' @seealso [hexify()] for creating HexData objects,
#'   [hex_distance()] for grid distances between cells
#'
#' @export
#' @examples
#' \donttest{
#' # ISEA grid neighbors
#' g <- hex_grid(area_km2 = 1000)
#' cell <- lonlat_to_cell(10, 50, g)
#' nbrs <- get_neighbors(cell, g)
#' nbrs[[1]]
#'
#' # H3 grid neighbors
#' g_h3 <- hex_grid(resolution = 5, type = "h3")
#' cell_h3 <- lonlat_to_cell(10, 50, g_h3)
#' get_neighbors(cell_h3, g_h3, k = 2)
#'
#' # With distances
#' get_neighbors(cell, g, k = 2, distances = TRUE)
#' }
get_neighbors <- function(cell_id, grid, k = 1L, include_self = FALSE,
                           distances = FALSE) {
  g <- extract_grid(grid)
  k <- as.integer(k)

  if (k < 0L) stop("k must be a non-negative integer")
  if (k == 0L) {
    if (include_self) {
      if (distances) {
        return(lapply(cell_id, function(cid) {
          data.frame(cell_id = cid, ring_distance = 0L)
        }))
      }
      return(as.list(cell_id))
    }
    n <- length(cell_id)
    if (distances) {
      empty <- data.frame(cell_id = cell_id[0], ring_distance = integer(0))
      return(rep(list(empty), n))
    }
    return(rep(list(cell_id[0]), n))
  }

  if (is_h3_grid(g)) {
    .get_neighbors_h3(cell_id, g, k, include_self, distances)
  } else {
    .get_neighbors_isea(cell_id, g, k, include_self, distances)
  }
}

# --- H3 backend ---
.get_neighbors_h3 <- function(cell_id, grid, k, include_self, distances) {
  if (distances) {
    result <- cpp_h3_gridDiskDistances(cell_id, k)
    if (!include_self) {
      result <- lapply(seq_along(result), function(i) {
        df <- result[[i]]
        if (nrow(df) > 0) {
          df[df$ring_distance > 0L, , drop = FALSE]
        } else {
          df
        }
      })
    }
    return(result)
  }

  result <- cpp_h3_gridDisk(cell_id, k)
  if (!include_self) {
    result <- mapply(function(nbrs, origin) {
      nbrs[nbrs != origin]
    }, result, cell_id, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  }
  result
}

# --- ISEA backend ---
.get_neighbors_isea <- function(cell_id, grid, k, include_self, distances) {
  ap <- as.integer(grid@aperture)
  res <- grid@resolution

  if (k == 1L && !distances) {
    # Fast path: direct C++ neighbor lookup
    result <- cpp_get_neighbors_isea(cell_id, res, ap)
    if (include_self) {
      result <- mapply(function(nbrs, origin) c(origin, nbrs),
                        result, cell_id, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    return(result)
  }

  # General k-ring: iterative BFS expansion
  lapply(cell_id, function(origin) {
    visited <- origin
    ring_dist <- 0L
    current_ring <- origin

    for (ring in seq_len(k)) {
      # Get neighbors of all cells in current ring
      all_nbrs <- cpp_get_neighbors_isea(current_ring, res, ap)
      new_cells <- unique(unlist(all_nbrs))
      new_cells <- setdiff(new_cells, visited)

      if (length(new_cells) == 0) break
      visited <- c(visited, new_cells)
      ring_dist <- c(ring_dist, rep(ring, length(new_cells)))
      current_ring <- new_cells
    }

    if (!include_self) {
      keep <- ring_dist > 0L
      visited <- visited[keep]
      ring_dist <- ring_dist[keep]
    }

    if (distances) {
      data.frame(cell_id = visited, ring_distance = ring_dist)
    } else {
      visited
    }
  })
}
