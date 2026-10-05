# Agreement of hexify's ISEA backend with DGGRID under orientations other than
# the standard one.
#
# DGGRID places the icosahedron either at a given vertex 0 and azimuth
# (dggs_orient_specify_type SPECIFIED) or about a region centre
# (REGION_CENTER). For each placement and each grid the same random points are
# assigned to cells by both programs, as in bench_dggrid_agreement.R: the
# script records how many cell IDs (DGGRID SEQNUMs) agree, the largest and
# median distance between the two programs' centres of the same cells, and the
# largest distance from a DGGRID cell corner to the nearest vertex of the
# hexify polygon on a random subset of cells. A REGION_CENTER grid is built in
# hexify with hex_grid(orientation = "region"), so its rows also test that
# hexify places the icosahedron where DGGRID does.
#
# Usage: Rscript paper/bench/bench_dggrid_orientation.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

N_POINTS <- 50000L
N_POLY <- 200L
SEED <- 20261005L

# One resolution per aperture configuration of CONFIGS.
GRIDS <- list(list(CONFIGS[[1]], 9), list(CONFIGS[[2]], 7), list(CONFIGS[[3]], 5),
              list(CONFIGS[[4]], 8), list(CONFIGS[[5]], 6))

# label, orientation or NULL, region centre or NULL
PLACEMENTS <- list(
  list("vert0 (-40, 20), azimuth 33", c(-40, 20, 33), NULL),
  list("vert0 (0, 90), azimuth 0", c(0, 90, 0), NULL),
  list("vert0 (150, -60), azimuth 200", c(150, -60, 200), NULL),
  list("region (10, 46.5)", NULL, c(10, 46.5)),
  list("region (-75, -10)", NULL, c(-75, -10)),
  list("region (0, 90)", NULL, c(0, 90))
)

corner_gap_km <- function(cd, hexify_geom) {
  ch <- sf::st_coordinates(hexify_geom)[, 1:2]
  max(vapply(seq_len(nrow(cd)), function(i)
    min(gc_km(cd[i, 1], cd[i, 2], ch[, 1], ch[, 2])), numeric(1)))
}

pts <- sphere_points(N_POINTS, seed = SEED)
rows <- list()
for (pl in PLACEMENTS) for (gr in GRIDS) {
  cfg <- gr[[1]]
  res <- gr[[2]]
  orient <- orient_lines(pl[[2]], pl[[3]])
  ap_arg <- if (length(cfg[[2]]) > 1) head(cfg[[2]], res) else cfg[[2]]
  g <- if (is.null(pl[[3]])) {
    hex_grid(resolution = res, aperture = ap_arg, orientation = pl[[2]])
  } else {
    hex_grid(resolution = res, aperture = ap_arg, orientation = "region", region = pl[[3]])
  }

  h_id <- lonlat_to_cell(pts$lon, pts$lat, g)
  d_id <- dggrid_seqnum(cfg, res, pts$lon, pts$lat, orient)
  same <- h_id == d_id

  ok_cells <- unique(h_id[same])
  set.seed(res)
  centre_cells <- if (length(ok_cells) > 20000) sample(ok_cells, 20000) else ok_cells
  centre_gap_m <- NA_real_
  cgap <- NA_real_
  if (length(centre_cells)) {
    hc <- cell_to_lonlat(centre_cells, g)
    dc <- dggrid_centres(cfg, res, centre_cells, orient)
    centre_gap_m <- gc_km(hc[[1]], hc[[2]], dc$lon, dc$lat) * 1000

    poly_cells <- sample(ok_cells, min(N_POLY, length(ok_cells)))
    dp <- dggrid_polygons(cfg, res, poly_cells, orient)
    hp <- cell_to_sf(dp$seqnum, g)
    hp <- hp[match(dp$seqnum, hp$cell_id), ]
    cgap <- vapply(seq_along(dp$seqnum), function(i)
      corner_gap_km(dp$corners[[i]], sf::st_geometry(hp)[i]), numeric(1)) * 1000
  }

  o <- g@orientation
  message(sprintf("%-30s ap %-12s res %2d: %d/%d agree, centre gap max %.3g m, corner gap %.3g m",
                  pl[[1]], cfg[[1]], res, sum(same), N_POINTS, max(centre_gap_m), max(cgap)))
  rows[[length(rows) + 1]] <- data.frame(
    placement = pl[[1]], vert0_lon = o[["vert0_lon"]], vert0_lat = o[["vert0_lat"]],
    azimuth = o[["azimuth"]], aperture = cfg[[1]], resolution = res,
    n_points = N_POINTS, n_agree = sum(same), n_differ = sum(!same),
    n_centre_cells = length(centre_cells),
    max_centre_gap_m = max(centre_gap_m), median_centre_gap_m = median(centre_gap_m),
    n_polygons = length(cgap), max_corner_gap_m = max(cgap))
}

dg_version <- tryCatch(system2("git", c("-C", shQuote(dirname(dirname(dirname(dirname(dirname(DGGRID_EXE)))))),
                                        "log", "-1", "--format=%h"), stdout = TRUE),
                       error = function(e) NA)
write_result(do.call(rbind, rows), "dggrid_orientation",
             extra = c(paste("n_points:", N_POINTS), paste("seed:", SEED),
                       paste("DGGRID executable:", DGGRID_EXE),
                       paste("DGGRID commit:", dg_version)))
