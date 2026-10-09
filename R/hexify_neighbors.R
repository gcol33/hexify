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
#' @param cell_id Cell IDs to find neighbors for: integer64 (or whole
#'   numbers) for ISEA grids, character for H3 grids.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#' @param k Integer. Ring distance (default 1). `k = 1` returns immediate
#'   neighbors, `k = 2` includes neighbors-of-neighbors, etc.
#' @param include_self Logical. If `TRUE`, include the input cell in the
#'   result (default `FALSE`).
#' @param distances Logical. If `TRUE`, return a data.frame with cell IDs
#'   and their ring distance from the origin (default `FALSE`).
#' @param as_sf Logical. If `TRUE`, return the links from each cell to its
#'   neighbours as great-circle lines between cell centres.
#' @param ring Logical. If `TRUE`, return only the cells exactly `k` hops
#'   away, the hollow ring, rather than the whole disk (default `FALSE`).
#'
#' @return If `distances = FALSE` (default): a list of cell ID vectors, one
#'   per input cell. If `distances = TRUE`: a list of data.frames with columns
#'   `cell_id` and `ring_distance`. If `as_sf = TRUE`: an sf object of
#'   lines (MULTILINESTRINGs when one is split at the antimeridian), one
#'   row per cell and neighbour, with columns `cell_id`,
#'   `neighbor_id` and `ring_distance`.
#'
#' @details
#' For **ISEA grids**, neighbors are computed using axial coordinate offsets
#' in the quad IJ space. On aperture-7 grids of the icosahedron, a cell at a
#' quad boundary takes its neighbors from digit arithmetic on its IGEO7 Z7
#' index, which is exact and needs no projection; on other ISEA grids such a
#' cell's neighbors are found by reprojection through lon/lat coordinates.
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
#'
#' # Only the cells two hops away
#' get_neighbors(cell, g, k = 2, ring = TRUE)
#'
#' # Links to the neighbours, drawn over the cells
#' links <- get_neighbors(cell, g, k = 2, as_sf = TRUE)
#' plot(sf::st_geometry(cell_to_sf(c(cell, links$neighbor_id), g)))
#' plot(sf::st_geometry(links), col = "red", add = TRUE)
#' }
get_neighbors <- function(cell_id, grid, k = 1L, include_self = FALSE,
                           distances = FALSE, as_sf = FALSE, ring = FALSE) {
  g <- extract_grid(grid)
  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  k <- as.integer(k)
  if (ring) {
    if (as_sf) {
      links <- neighbor_links(cell_id, g, k, include_self = k == 0L)
      return(links[links$ring_distance == k, ])
    }
    disk <- get_neighbors(cell_id, g, k, include_self = k == 0L, distances = TRUE)
    disk <- lapply(disk, function(d) d[d$ring_distance == k, , drop = FALSE])
    if (distances) return(disk)
    return(lapply(disk, `[[`, "cell_id"))
  }
  if (as_sf) return(neighbor_links(cell_id, g, k, include_self))

  if (k < 0L) stop("k must be a non-negative integer")
  if (k == 0L) {
    if (include_self) {
      if (distances) {
        return(lapply(seq_along(cell_id), function(i) {
          data.frame(cell_id = cell_id[i], ring_distance = 0L)
        }))
      }
      return(lapply(seq_along(cell_id), function(i) cell_id[i]))
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
  # k = 1 fast path. C++ walks every aperture and every mixed sequence.
  get_k1 <- function(cids) grid_neighbors_isea(cids, grid)

  if (k == 1L && !distances) {
    result <- get_k1(cell_id)
    if (include_self) {
      result <- lapply(seq_along(result), function(i) c(cell_id[i], result[[i]]))
    }
    return(result)
  }

  lapply(seq_along(cell_id), function(i) {
    rings <- isea_rings(cell_id[i], grid, k)
    keep <- include_self | rings$ring_distance > 0L
    if (distances) {
      data.frame(cell_id = rings$cell_id[keep],
                 ring_distance = rings$ring_distance[keep])
    } else {
      rings$cell_id[keep]
    }
  })
}

#' Cells within k rings of a cell, by breadth-first search
#'
#' Each ring holds the neighbours of the previous ring not met before, in the
#' order the neighbour lookup returns them. The search stops early once every
#' cell of `targets` is reached, or when a ring comes back empty.
#'
#' @param source One cell ID
#' @param g HexGridInfo object of an ISEA grid
#' @param k Number of rings
#' @param targets Optional cell IDs whose rings are wanted
#' @return List with `cell_id` (source first, then ring by ring) and the
#'   integer `ring_distance` of each
#' @noRd
isea_rings <- function(source, g, k, targets = NULL) {
  visited <- source
  ring_distance <- 0L
  frontier <- source
  remaining <- if (!is.null(targets)) setdiff(targets, source)

  for (ring in seq_len(k)) {
    if (!is.null(targets) && length(remaining) == 0L) break
    frontier <- setdiff(unique(cell_id_unlist(grid_neighbors_isea(frontier, g))),
                        visited)
    if (length(frontier) == 0L) break
    visited <- c(visited, frontier)
    ring_distance <- c(ring_distance, rep(ring, length(frontier)))
    if (!is.null(targets)) remaining <- remaining[!remaining %in% frontier]
  }

  list(cell_id = visited, ring_distance = ring_distance)
}

#' Links from cells to their neighbours as sf lines
#'
#' Each link is the great-circle arc between the two cell centres, with a
#' point at least every degree, its longitudes continuous and split at the
#' antimeridian. When a link is split, every link becomes a MULTILINESTRING.
#' @noRd
neighbor_links <- function(cell_id, g, k, include_self) {
  rings <- get_neighbors(cell_id, g, k, include_self, distances = TRUE)
  from <- rep(cell_id, vapply(rings, nrow, integer(1)))
  to <- lapply(rings, `[[`, "cell_id")
  to <- if (is_h3_grid(g)) unlist(to, use.names = FALSE) else cell_id_unlist(to)
  ring_distance <- unlist(lapply(rings, `[[`, "ring_distance"), use.names = FALSE)

  ids <- unique(c(from, to))
  ctr <- cell_to_lonlat(ids, g)
  V <- unit_vec(ctr$lon_deg, ctr$lat_deg)
  a <- V[match(from, ids), , drop = FALSE]
  b <- V[match(to, ids), , drop = FALSE]
  lines <- lapply(seq_along(from), function(i) {
    P <- rbind(slerp(a[i, ], b[i, ], pi / 180), b[i, ])
    ll <- vec_lonlat(P)
    ll[, 1] <- ll[1, 1] + cumsum(c(0, (diff(ll[, 1]) + 180) %% 360 - 180))
    sf::st_linestring(ll)
  })
  links <- sf::st_sf(cell_id = from, neighbor_id = to,
                     ring_distance = as.integer(ring_distance),
                     geometry = sf::st_sfc(lines, crs = grid_crs(g)))
  links <- wrap_cells_at_dateline(links)
  if (inherits(sf::st_geometry(links), "sfc_GEOMETRY")) {
    links <- sf::st_cast(links, "MULTILINESTRING")
  }
  links
}
