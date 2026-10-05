# Cell area against latitude for three grids of about 12,400 km2 per cell:
# ISEA aperture 4 resolution 6 (hexify), H3 resolution 3 (hexify's H3 backend)
# and a 1 x 1 degree longitude-latitude grid.
#
# Hexagon areas are measured independently of hexify's own area functions:
# each cell polygon returned by cell_to_sf() is measured on the sphere with s2.
# Longitude-latitude cell areas use the exact spherical formula
# R^2 * dlon * (sin(lat2) - sin(lat1)). Pentagon cells are excluded.
#
# Usage: Rscript paper/bench/bench_cell_area.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

R_KM <- 6371.007180918475
sf::sf_use_s2(TRUE)
set.seed(20261005)
lats <- seq(-89.5, 89.5, by = 0.5)
pts <- data.frame(lat = rep(lats, each = 20), lon = runif(length(lats) * 20, -180, 180))

hex_areas <- function(g, label) {
  id <- suppressMessages(lonlat_to_cell(pts$lon, pts$lat, g))
  keep <- !duplicated(id) & !is_pentagon(id, g)
  cells <- id[keep]
  poly <- cell_to_sf(cells, g)
  poly <- poly[match(cells, poly$cell_id), ]
  ctr <- cell_to_lonlat(cells, g)
  data.frame(grid = label, cell = as.character(cells), lat = ctr[[2]],
             area_km2 = as.numeric(sf::st_area(poly)) / 1e6)
}

isea <- hex_areas(hex_grid(resolution = 6, aperture = 4), "ISEA aperture 4, res 6")
h3 <- hex_areas(hex_grid(resolution = 3, type = "h3"), "H3 res 3")

band <- seq(-90, 89)
ll <- data.frame(grid = "1 x 1 degree", cell = as.character(band), lat = band + 0.5,
                 area_km2 = R_KM^2 * (pi / 180) * (sin((band + 1) * pi / 180) - sin(band * pi / 180)))

out <- rbind(isea, h3, ll)
summ <- do.call(rbind, lapply(split(out, out$grid), function(d)
  data.frame(grid = d$grid[1], n = nrow(d), mean_km2 = mean(d$area_km2),
             min_km2 = min(d$area_km2), max_km2 = max(d$area_km2),
             max_over_min = max(d$area_km2) / min(d$area_km2),
             cv_pct = 100 * sd(d$area_km2) / mean(d$area_km2))))
print(summ)
write_result(out, "cell_area_by_latitude")
write_result(summ, "cell_area_summary")
