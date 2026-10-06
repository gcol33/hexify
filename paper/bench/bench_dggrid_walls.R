# The cell wall midpoint ratio (Gregory et al. 2008) of hexify's and DGGRID's
# aperture-7 cells at odd resolutions, where DGGRID places two corners of each
# cell straddling a face edge apart from where its neighbours place them.
#
# Both programs' cells are read the same way: corners only, walls as
# great-circle arcs between them, each program's own cell centres, and the
# neighbours hexify gives (cell numbers agree between the two). A wall's
# ratio is the distance from its midpoint to the midpoint of the arc joining
# the two centres, over its length. Every cell of the grid is measured.
#
# A cell has a corner gap where some DGGRID corner lies more than 1 m from
# every hexify corner of the same cell (CORNER_GAP_KM). The script counts the
# DGGRID walls whose ratio exceeds the largest hexify ratio, and how many of
# them belong to cells with a corner gap.
#
# Usage: Rscript paper/bench/bench_dggrid_walls.R

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "bench_common.R"))
source(file.path(here, "bench_dggrid_common.R"))

cfg <- CONFIGS[[3]]
stopifnot(identical(cfg[[1]], "7"))

wall_ratios <- function(rings, centres, nb) {
  cm <- as.matrix(centres)
  w <- hexify:::cpp_ring_walls(rings, cm, lapply(nb, function(x) cm[x, , drop = FALSE]))
  # Rows run cell by cell, each cell's neighbours in the order given
  stopifnot(nrow(w) == length(unlist(nb)))
  w$neighbor_id <- unlist(nb)
  w$ratio <- w$midpoint_offset / w$wall
  w
}

one_resolution <- function(res) {
  g <- hex_grid(resolution = res, aperture = 7)
  ids <- seq_len(hexify:::grid_n_cells(g))
  lv <- hexify:::isea_levels(g@aperture, g@resolution)
  nb <- hexify:::grid_neighbors_isea(ids, g)

  hr <- hexify:::cpp_cell_to_corners(hexify:::icosa_arg(g), ids, lv$resolution,
                                     lv$aperture, lv$ap_seq, 0)
  hc <- cell_to_lonlat(ids, g)
  dg <- dggrid_generate(cfg, res, ids, cells = TRUE)
  dr <- dg$cells$corners[match(ids, dg$cells$seqnum)]
  dc <- dg$centres[match(ids, dg$centres$seqnum), c("lon", "lat")]
  stopifnot(!anyNA(match(ids, dg$cells$seqnum)))

  wh <- wall_ratios(hr, hc, nb)
  wd <- wall_ratios(dr, dc, nb)
  stopifnot(identical(wh$cell, wd$cell), identical(wh$neighbor_id, wd$neighbor_id))

  gap <- vapply(ids, function(i) corner_gap_km(dr[[i]], hr[[i]]), numeric(1)) > CORNER_GAP_KM
  hmax <- max(wh$ratio)
  over <- wd$ratio > hmax
  over_gap <- over & (gap[wd$cell] | gap[wd$neighbor_id])

  q <- function(x, p) unname(quantile(x, p))
  rbind(
    data.frame(resolution = res, program = "hexify", n_cells = length(ids),
               n_walls = nrow(wh), n_gap_cells = sum(gap),
               ratio_median = median(wh$ratio), ratio_q999 = q(wh$ratio, 0.999),
               ratio_max = hmax, n_walls_over_hexify_max = 0L,
               n_over_at_gap_cells = 0L),
    data.frame(resolution = res, program = "DGGRID", n_cells = length(ids),
               n_walls = nrow(wd), n_gap_cells = sum(gap),
               ratio_median = median(wd$ratio), ratio_q999 = q(wd$ratio, 0.999),
               ratio_max = max(wd$ratio), n_walls_over_hexify_max = sum(over),
               n_over_at_gap_cells = sum(over_gap)))
}

out <- do.call(rbind, lapply(c(3L, 5L), one_resolution))
print(out, digits = 4)
write_result(out, "dggrid_wall_ratio")
