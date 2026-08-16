# =============================================================================
# Hierarchical navigation for mixed aperture grids
# =============================================================================
#
# Mixed grids have no DGGRID-standard hierarchical index (DGGRID rejects
# hierarchical indexing for any non-"PURE" aperture), so there is no Z7/Z3-style
# string to port. hexify instead defines the hierarchy geometrically, which is
# the mathematically correct relationship for how these grids are built:
#
# quad_xy_to_ij_mixed() quantizes the projected point once, by a single scalar
# scale, the product of sqrt(aperture) over the levels. A cell at resolution r
# and its parent at r-1 are therefore both direct quantizations of the SAME
# continuous quad-space at different scales -- so a cell's parent is simply the
# coarser cell whose lattice point contains the cell's centre
# (center-containment), with no accumulated per-step rotation to go wrong.
# Parent/children/index all follow from re-projecting a cell centre through the
# validated forward pipeline.
#
# Every helper takes the aperture spelling alongside the resolution, since the
# coarser grid a parent lives on is the same spelling read at that resolution
# (see aperture_at_resolution()).

mixed_ap_seq <- function(aperture, resolution) {
  parse_aperture_seq(aperture_at_resolution(aperture, resolution), resolution)
}

mixed_cell_center <- function(cell_id, resolution, aperture) {
  cpp_cell_to_lonlat_seq(as.numeric(cell_id), mixed_ap_seq(aperture, resolution))
}

mixed_point_to_cell <- function(lon, lat, resolution, aperture) {
  cpp_lonlat_to_cell_seq(as.numeric(lon), as.numeric(lat),
                         mixed_ap_seq(aperture, resolution))
}

mixed_cell_qij <- function(cell_id, resolution, aperture) {
  cpp_cell_to_quad_ij_seq(as.numeric(cell_id), mixed_ap_seq(aperture, resolution))
}

mixed_qij_cell <- function(quad, i, j, resolution, aperture) {
  cpp_quad_ij_to_cell_seq(as.integer(quad), as.numeric(i), as.numeric(j),
                          mixed_ap_seq(aperture, resolution))
}

#' Maximum (i,j) substrate coordinate for a mixed grid at a resolution.
#' One less than the quad edge coordinate the C++ layer quantizes against.
#' @noRd
mixed_max_dim <- function(resolution, aperture) {
  if (resolution == 0) return(0L)
  as.integer(cpp_ap_seq_edge_dim(mixed_ap_seq(aperture, resolution))) - 1L
}

#' All valid cells within a band of the (i,j) boundary of the ten body quads,
#' plus the poles. Non-nested ISEA seams -- especially near the twelve
#' icosahedron vertices (pentagon points), where a parent's children spread into
#' a polar cap rather than onto a single edge line -- put children a few cells
#' inside the quad edge, so the band has width `w` (not just the edge line). The
#' set is O(w * sqrt(n_cells)), still cheap against the full grid.
#' @noRd
mixed_boundary_cells <- function(resolution, aperture, n_cells, w = 4L) {
  m <- mixed_max_dim(resolution, aperture)
  w <- min(w, as.integer(m) + 1L)
  lo <- 0:(w - 1L)
  hi <- (m - w + 1L):m
  band <- unique(c(lo, hi))            # (i,j) indices near either extreme
  band <- band[band >= 0 & band <= m]
  full <- 0:m
  # Two vertical strips (i in band, any j) and two horizontal strips (j in band).
  gi <- c(rep(band, times = length(full)), rep(full, times = length(band)))
  gj <- c(rep(full, each = length(band)), rep(band, each = length(full)))
  ne <- length(gi)
  q <- rep(1:10, each = ne)
  ii <- rep(gi, times = 10)
  jj <- rep(gj, times = 10)
  cells <- mixed_qij_cell(q, ii, jj, resolution, aperture)
  cells <- c(cells,
             mixed_point_to_cell(0, 90, resolution, aperture),
             mixed_point_to_cell(0, -90, resolution, aperture))
  cells <- unique(cells)
  cells <- cells[is.finite(cells) & cells >= 1 & cells <= n_cells]
  # Keep only (i,j) that round-trip to a real cell.
  if (length(cells) == 0) return(cells)
  q2 <- mixed_cell_qij(cells, resolution, aperture)
  rt <- mixed_qij_cell(q2$quad, q2$i, q2$j, resolution, aperture)
  cells[rt == cells]
}

#' Geometric parent of mixed cells (center-containment)
#' @noRd
mixed_get_parent <- function(cell_id, resolution, aperture, levels = 1L) {
  parent_res <- resolution - as.integer(levels)
  ll <- mixed_cell_center(cell_id, resolution, aperture)
  mixed_point_to_cell(ll$lon_deg, ll$lat_deg, parent_res, aperture)
}

#' Candidate neighbour cells of a set of mixed cells.
#'
#' Two sources unioned for robustness: exact axial (i,j) neighbours in each
#' cell's own quad (correct for quad interiors), plus directional lon/lat probes
#' at the cell spacing (which cross icosahedron seams into adjacent quads). Only
#' used to grow a candidate superset; correctness comes from the parent filter.
#' @noRd
mixed_neighbors <- function(cells, child_res, aperture, spacing_deg, n_cells_child) {
  qij <- mixed_cell_qij(cells, child_res, aperture)
  # 6 axial hex-neighbour offsets
  off <- rbind(c(1, 0), c(0, 1), c(-1, 1), c(-1, 0), c(0, -1), c(1, -1))
  nb <- integer(0)
  for (r in seq_len(nrow(off))) {
    nb <- c(nb, mixed_qij_cell(qij$quad, qij$i + off[r, 1], qij$j + off[r, 2],
                               child_res, aperture))
  }
  # Directional probes at the cell spacing (seam-crossing)
  cll <- mixed_cell_center(cells, child_res, aperture)
  if (is.finite(spacing_deg) && spacing_deg > 0) {
    ang <- seq(0, 2 * pi, length.out = 13)[-13]
    coslat <- cos(cll$lat_deg * pi / 180)
    for (rad in c(0.8, 1.15) * spacing_deg) {
      for (a in ang) {
        plat <- pmax(pmin(cll$lat_deg + rad * sin(a), 89.9999), -89.9999)
        dlon <- ifelse(coslat > 1e-6, rad * cos(a) / coslat, rad * cos(a))
        plon <- ((cll$lon_deg + dlon + 180) %% 360) - 180
        nb <- c(nb, mixed_point_to_cell(plon, plat, child_res, aperture))
      }
    }
  }
  nb <- unique(nb)
  nb <- nb[is.finite(nb) & nb >= 1 & nb <= n_cells_child]
  # Keep only cells whose (i,j) round-trips to a real cell
  if (length(nb) == 0) return(nb)
  q2 <- mixed_cell_qij(nb, child_res, aperture)
  rt <- mixed_qij_cell(q2$quad, q2$i, q2$j, child_res, aperture)
  nb[rt == nb]
}

#' Geometric children of a single mixed cell.
#'
#' Children are cells whose geometric parent (centre re-quantised at the coarser
#' resolution) is this cell -- an exact test. The work is producing a complete
#' candidate superset to run that test on. At coarse resolutions (few, huge
#' cells) every child cell is tested (cheap and exhaustive). Otherwise children
#' form a connected patch, so a breadth-first walk over neighbours, seeded from
#' an (i,j) box and grown only through cells that pass the parent test, visits
#' them all without relying on sampling density. Returns child IDs sorted.
#' @noRd
mixed_get_children_one <- function(cell_id, resolution, child_res, aperture,
                                   n_cells_child) {
  # Coarse levels: exhaustive filter (no locality assumptions, cheap).
  if (n_cells_child <= 2000) {
    all_child <- seq_len(n_cells_child)
    par <- mixed_get_parent(all_child, child_res, aperture, child_res - resolution)
    return(sort(all_child[par == cell_id]))
  }

  spacing_deg <- mixed_spacing_deg(child_res, aperture)
  ll <- mixed_cell_center(cell_id, resolution, aperture)

  # Seed set: (i,j) box around the central child, plus the pole cells.
  central <- mixed_point_to_cell(ll$lon_deg, ll$lat_deg, child_res, aperture)
  cq <- mixed_cell_qij(central, child_res, aperture)
  box <- 3L
  di <- rep(-box:box, times = 2 * box + 1)
  dj <- rep(-box:box, each = 2 * box + 1)
  seed <- mixed_qij_cell(cq$quad, cq$i + di, cq$j + dj, child_res, aperture)
  seed <- c(seed, central,
            mixed_point_to_cell(0, 90, child_res, aperture),
            mixed_point_to_cell(0, -90, child_res, aperture))
  seed <- unique(seed)
  seed <- seed[is.finite(seed) & seed >= 1 & seed <= n_cells_child]

  parent_of <- function(x) {
    mixed_get_parent(x, child_res, aperture, child_res - resolution)
  }

  found <- seed[parent_of(seed) == cell_id]
  if (length(found) == 0) {
    # Fall back to a broader probe if the box seed missed entirely.
    found <- central[parent_of(central) == cell_id]
  }
  visited <- unique(c(seed, found))
  frontier <- found

  # BFS over neighbours, keeping only genuine children.
  for (iter in seq_len(64)) {
    if (length(frontier) == 0) break
    nb <- mixed_neighbors(frontier, child_res, aperture, spacing_deg, n_cells_child)
    nb <- setdiff(nb, visited)
    visited <- c(visited, nb)
    if (length(nb) == 0) break
    kids <- nb[parent_of(nb) == cell_id]
    found <- c(found, kids)
    frontier <- kids
  }

  # Seam children: non-nested ISEA seams (and the twelve icosahedron vertices,
  # which are quad corners) put a few children on a distant quad boundary that
  # the interior walk cannot reach. This can only happen when the PARENT's own
  # footprint touches a seam, i.e. its (i,j) is near a quad edge or corner. For
  # such parents -- an O(sqrt(n)) minority -- sweep the O(sqrt(n)) boundary band
  # and add any cell whose exact parent is this cell. Interior parents keep all
  # children in-quad, so the fast walk above is already complete for them.
  pq <- mixed_cell_qij(cell_id, resolution, aperture)
  m_parent <- mixed_max_dim(resolution, aperture)
  gate <- 3L
  near_boundary <- min(pq$i, pq$j, m_parent - pq$i, m_parent - pq$j) < gate
  if (near_boundary) {
    bnd <- mixed_boundary_cells(child_res, aperture, n_cells_child)
    if (length(bnd) > 0) {
      found <- c(found, bnd[parent_of(bnd) == cell_id])
    }
  }

  sort(unique(found))
}

#' Approximate cell centre spacing in degrees at a resolution, for sizing the
#' child-footprint sampling disk.
#' @noRd
mixed_spacing_deg <- function(resolution, aperture) {
  # Mean cell area (km^2) from the grid's cell count over the body's surface,
  # then hexagon centre spacing s = sqrt(2 A / sqrt(3)) in degrees of arc. Area
  # goes with r^2 and both the spacing and the km per degree with r, so the
  # angle is the same on every body and Earth's figures give it.
  n <- aperture_n_cells(aperture, resolution)
  area_km2 <- EARTH_SURFACE_KM2 / n
  spacing_km <- sqrt(2 * area_km2 / sqrt(3))
  spacing_km / KM_PER_DEGREE
}

# -----------------------------------------------------------------------------
# Hierarchical index string
# -----------------------------------------------------------------------------
# The index is prefix-hierarchical: a 2-digit base cell (the resolution-0 cell,
# 1..12) followed by a 2-digit ordinal per resolution level, giving which child
# of its parent the cell is (children sorted ascending by cell ID). Two digits
# per level because near the twelve icosahedron vertices centre-containment can
# gather more than nine children into one coarse cell. Dropping the last two
# digits yields the parent's index, so the string sorts by spatial ancestry.
mixed_index_digit_width <- 2L

#' Encode a single mixed cell to its hierarchical index string.
#' @noRd
mixed_cell_to_index_one <- function(cell_id, resolution, aperture) {
  if (resolution == 0) {
    return(sprintf("%02d", as.integer(cell_id)))
  }
  # Ancestor chain a[0..resolution]; a[resolution] = cell_id.
  anc <- numeric(resolution + 1)
  anc[resolution + 1] <- cell_id
  for (k in resolution:1) {
    anc[k] <- mixed_get_parent(anc[k + 1], k, aperture, 1L)
  }
  w <- mixed_index_digit_width
  digits <- integer(resolution)
  for (k in 1:resolution) {
    kids <- mixed_get_children_one(anc[k], k - 1L, k, aperture,
                                   aperture_n_cells(aperture, k))
    pos <- match(anc[k + 1], kids)
    if (is.na(pos)) {
      stop("hexify internal error: mixed child enumeration missed a descendant ",
           "(cell ", format(cell_id, scientific = FALSE), ", level ", k, ")")
    }
    digits[k] <- pos - 1L
  }
  if (any(digits >= 10^w)) {
    stop("hexify internal error: mixed child ordinal exceeds index digit width")
  }
  paste0(sprintf("%02d", as.integer(anc[1])),
         paste(sprintf(paste0("%0", w, "d"), digits), collapse = ""))
}

#' Decode a single mixed hierarchical index string to a cell ID.
#' @noRd
mixed_index_to_cell_one <- function(index, resolution, aperture) {
  base0 <- as.integer(substr(index, 1, 2))
  cell <- base0
  if (resolution == 0) return(cell)
  w <- mixed_index_digit_width
  for (k in 1:resolution) {
    start <- 2L + (k - 1L) * w + 1L
    d <- as.integer(substr(index, start, start + w - 1L))
    kids <- mixed_get_children_one(cell, k - 1L, k, aperture,
                                   aperture_n_cells(aperture, k))
    if (d + 1L > length(kids)) {
      stop("hexify internal error: mixed index digit out of range on decode")
    }
    cell <- kids[d + 1L]
  }
  cell
}
