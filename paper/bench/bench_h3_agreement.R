# Agreement of hexify's H3 backend with the H3 reference library (through h3r,
# which links the h3lib build of Uber's C library). The same random points are
# indexed by both at several resolutions; the script records the number of
# identical cell IDs and the largest distance between the two cell centres.
#
# Usage: Rscript paper/bench/bench_h3_agreement.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

N_POINTS <- 200000L
RESOLUTIONS <- c(0, 3, 6, 9, 12, 15)
pts <- sphere_points(N_POINTS, seed = 20261005)

rows <- lapply(RESOLUTIONS, function(res) {
  g <- hex_grid(resolution = res, type = "h3")
  h_id <- suppressMessages(lonlat_to_cell(pts$lon, pts$lat, g))
  r_id <- h3r::latLngToCell(pts$lat, pts$lon, as.integer(res))
  same <- tolower(h_id) == tolower(r_id)
  cells <- unique(h_id)
  hc <- cell_to_lonlat(cells, g)
  rc <- h3r::cellToLatLng(cells)
  gap_m <- max(gc_km(hc[[1]], hc[[2]], rc$lng, rc$lat)) * 1000
  message(sprintf("H3 res %2d: %d/%d agree, centre gap %.3g m", res, sum(same), N_POINTS, gap_m))
  data.frame(resolution = res, n_points = N_POINTS, n_agree = sum(same),
             n_cells = length(cells), max_centre_gap_m = gap_m)
})

write_result(do.call(rbind, rows), "h3_agreement",
             extra = c(paste("n_points:", N_POINTS), "seed: 20261005"))
