# Shape and spacing metrics of whole grids: compactness (White et al. 1998),
# intercell distance and the cell wall midpoint ratio (Gregory et al. 2008),
# for apertures 3, 4 and 7 on Snyder's and Fuller's projections and for H3,
# every cell and every wall of each grid. Area spread is given beside them,
# since Fuller's projection and H3 are not equal-area.
#
# Resolutions are chosen so that each grid has some 5,000 to 25,000 cells.
#
# Usage: Rscript paper/bench/bench_grid_metrics.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

GRIDS <- list(
  list("ISEA3H", 3, 6, "isea"), list("FULLER3H", 3, 6, "fuller"),
  list("ISEA4H", 4, 5, "isea"), list("FULLER4H", 4, 5, "fuller"),
  list("ISEA7H", 7, 4, "isea"), list("FULLER7H", 7, 4, "fuller"),
  list("H3", NA, 2, NA)
)

summarise_grid <- function(spec) {
  g <- if (spec[[1]] == "H3") {
    hex_grid(resolution = spec[[3]], type = "h3")
  } else {
    hex_grid(resolution = spec[[3]], aperture = spec[[2]], projection = spec[[4]])
  }
  t0 <- proc.time()[["elapsed"]]
  m <- cell_metrics(grid = g)
  w <- wall_metrics(grid = g)
  secs <- proc.time()[["elapsed"]] - t0
  cv <- function(x) sd(x) / mean(x)
  data.frame(
    grid = spec[[1]], resolution = spec[[3]], n_cells = nrow(m), n_walls = nrow(w),
    area_cv = cv(m$area_km2),
    area_max_over_min = max(m$area_km2) / min(m$area_km2),
    compactness_mean = mean(m$compactness),
    compactness_min = min(m$compactness),
    compactness_max = max(m$compactness),
    distance_cv = cv(w$center_distance_km),
    distance_range_over_mean = diff(range(w$center_distance_km)) / mean(w$center_distance_km),
    midpoint_ratio_mean = mean(w$midpoint_ratio),
    midpoint_ratio_max = max(w$midpoint_ratio),
    seconds = secs)
}

out <- do.call(rbind, lapply(GRIDS, summarise_grid))
print(out, digits = 4)
write_result(out, "grid_metrics")
