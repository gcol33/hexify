# Run time of point assignment and polygon generation, single R process.
#
# Point assignment: one million points uniform on the sphere, assigned to an
# ISEA aperture-3 resolution-10 grid by hexify and by dggridR, and to an H3
# resolution-7 grid by hexify's H3 backend, h3r and h3o. Polygon generation:
# 10,000 ISEA cells by hexify (cell_to_sf) and dggridR (dgcellstogrid).
# Each task runs REPS times; the median elapsed time is reported.
#
# Usage: Rscript paper/bench/bench_speed.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))
suppressPackageStartupMessages(library(dggridR))

N_POINTS <- 1e6
N_POLY <- 10000L
REPS <- 7L
pts <- sphere_points(N_POINTS, seed = 20261005)

g3 <- hex_grid(resolution = 10, aperture = 3)
d3 <- dgconstruct(res = 10, aperture = 3, show_info = FALSE)
gh <- hex_grid(resolution = 7, type = "h3")
invisible(suppressMessages(lonlat_to_cell(0, 0, gh)))

set.seed(1)
cells <- unique(lonlat_to_cell(pts$lon[1:(3 * N_POLY)], pts$lat[1:(3 * N_POLY)], g3))[1:N_POLY]
xy <- cbind(pts$lon, pts$lat)

tasks <- list(
  list("assign", "ISEA ap3 res10", "hexify", function() lonlat_to_cell(pts$lon, pts$lat, g3)),
  list("assign", "ISEA ap3 res10", "dggridR", function() dgGEO_to_SEQNUM(d3, pts$lon, pts$lat)$seqnum),
  list("assign", "H3 res7", "hexify", function() lonlat_to_cell(pts$lon, pts$lat, gh)),
  list("assign", "H3 res7", "h3r", function() h3r::latLngToCell(pts$lat, pts$lon, 7L)),
  list("assign", "H3 res7", "h3o", function() h3o::h3_from_xy(pts$lon, pts$lat, 7L)),
  list("polygons", "ISEA ap3 res10", "hexify", function() cell_to_sf(cells, g3)),
  list("polygons", "ISEA ap3 res10", "dggridR", function() dgcellstogrid(d3, cells))
)

rows <- lapply(tasks, function(t) {
  invisible(t[[4]]())
  el <- vapply(seq_len(REPS), function(i) system.time(t[[4]](), gcFirst = TRUE)[["elapsed"]],
               numeric(1))
  n <- if (t[[1]] == "assign") N_POINTS else N_POLY
  message(sprintf("%-8s %-15s %-8s median %.3f s", t[[1]], t[[2]], t[[3]], median(el)))
  data.frame(task = t[[1]], grid = t[[2]], package = t[[3]], n = n,
             median_s = median(el), min_s = min(el), max_s = max(el),
             per_second = n / median(el))
})

write_result(do.call(rbind, rows), "speed",
             extra = c(paste("reps:", REPS), paste("n_points:", N_POINTS),
                       paste("n_polygons:", N_POLY)))
