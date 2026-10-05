# Agreement of hexify's ISEA backend with DGGRID, the reference implementation
# of the ISEA discrete global grid systems (Sahr et al. 2003).
#
# DGGRID is run as its command-line program, built from source, so that every
# aperture hexify supports can be compared: 3, 4, 7, the ISEA43H family and a
# mixed per-level sequence (dggridR exposes apertures 3 and 4 only). For each
# grid the same random points are assigned to cells by both programs. The
# script records how many cell IDs (DGGRID SEQNUMs) agree, the largest and
# median distance between the two programs' centres of the same cells, and, for each point
# assigned differently, its distance to the boundary of the hexify cell. Cell
# boundaries are compared on a random subset of cells as the largest distance
# from a DGGRID corner to the nearest vertex of the hexify polygon.
#
# Usage: Rscript paper/bench/bench_dggrid_agreement.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

N_POINTS <- 200000L
N_POLY <- 500L
SEED <- 20261005L

corner_gap_km <- function(cd, hexify_geom) {
  ch <- sf::st_coordinates(hexify_geom)[, 1:2]
  max(vapply(seq_len(nrow(cd)), function(i)
    min(gc_km(cd[i, 1], cd[i, 2], ch[, 1], ch[, 2])), numeric(1)))
}

pts <- sphere_points(N_POINTS, seed = SEED)
rows <- list()
differing <- list()
for (cfg in CONFIGS) for (res in cfg[[4]]) {
  ap_arg <- if (length(cfg[[2]]) > 1) head(cfg[[2]], res) else cfg[[2]]
  g <- hex_grid(resolution = res, aperture = ap_arg)
  h_id <- lonlat_to_cell(pts$lon, pts$lat, g)
  d_id <- dggrid_seqnum(cfg, res, pts$lon, pts$lat)
  same <- h_id == d_id

  ok_cells <- unique(h_id[same])
  set.seed(res)
  centre_cells <- if (length(ok_cells) > 20000) sample(ok_cells, 20000) else ok_cells
  hc <- cell_to_lonlat(centre_cells, g)
  dc <- dggrid_centres(cfg, res, centre_cells)
  centre_gap_m <- gc_km(hc[[1]], hc[[2]], dc$lon, dc$lat) * 1000

  diff_idx <- which(!same)
  edge_m <- NA_real_
  if (length(diff_idx)) {
    differing[[length(differing) + 1]] <- data.frame(
      aperture = cfg[[1]], resolution = res, lon = pts$lon[diff_idx], lat = pts$lat[diff_idx],
      hexify_cell = h_id[diff_idx], dggrid_cell = d_id[diff_idx])
    dpts <- sf::st_as_sf(pts[diff_idx, ], coords = c("lon", "lat"), crs = 4326)
    hpoly <- cell_to_sf(h_id[diff_idx], g)
    bnd <- sf::st_boundary(hpoly[match(h_id[diff_idx], hpoly$cell_id), ])
    edge_m <- max(as.numeric(sf::st_distance(dpts, bnd, by_element = TRUE)))
  }

  poly_cells <- sample(ok_cells, min(N_POLY, length(ok_cells)))
  dp <- dggrid_polygons(cfg, res, poly_cells)
  hp <- cell_to_sf(dp$seqnum, g)
  hp <- hp[match(dp$seqnum, hp$cell_id), ]
  cgap <- vapply(seq_along(dp$seqnum), function(i)
    corner_gap_km(dp$corners[[i]], sf::st_geometry(hp)[i]), numeric(1))

  message(sprintf("ap %-12s res %2d: %d/%d agree, centre gap max %.3g m, corner gap %.3g m, edge distance %.3g m",
                  cfg[[1]], res, sum(same), N_POINTS, max(centre_gap_m), max(cgap) * 1000, edge_m))
  rows[[length(rows) + 1]] <- data.frame(
    aperture = cfg[[1]], resolution = res, cell_area_km2 = g@area_km2,
    n_points = N_POINTS, n_agree = sum(same), n_differ = sum(!same),
    n_centre_cells = length(centre_cells),
    max_centre_gap_m = max(centre_gap_m), median_centre_gap_m = median(centre_gap_m),
    max_edge_distance_m = edge_m, n_polygons = length(dp$seqnum),
    max_corner_gap_m = max(cgap) * 1000)
}

dg_version <- tryCatch(system2("git", c("-C", shQuote(dirname(dirname(dirname(dirname(dirname(DGGRID_EXE)))))),
                                        "log", "-1", "--format=%h"), stdout = TRUE),
                       error = function(e) NA)
write_result(do.call(rbind, rows), "dggrid_agreement",
             extra = c(paste("n_points:", N_POINTS), paste("seed:", SEED),
                       paste("DGGRID executable:", DGGRID_EXE),
                       paste("DGGRID commit:", dg_version)))
write_result(do.call(rbind, differing), "dggrid_disagreements",
             extra = c(paste("seed:", SEED), "points the two programs assign to different cells"))
