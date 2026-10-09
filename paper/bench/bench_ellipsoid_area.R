# Cell areas on the WGS84 ellipsoid, measured geodesically.
#
# Every cell of ISEA3H, ISEA4H and ISEA7H grids is traced with hexify's own
# boundary densification (relative chord tolerance 1e-4, and no piece longer
# than 1e-3 radians, which puts the measured areas within about 1e-6 of their
# converged values), and the ring's area is measured on WGS84 by GeographicLib's geodesic polygon area (Karney 2013,
# through geosphere::areaPolygon()). Each area is divided by the area hexify
# reports for the cell, cell_area(), which on Earth shares the ellipsoid's
# surface among the cells. The default grid's rings are also measured on the
# sphere of the grid's radius, where that ratio is 1 by construction, so its
# spread is the measurement's own error.
#
# The default grid reads geodetic latitude as latitude on the sphere; a grid
# built with ellipsoid = "WGS84" converts it to authalic latitude, which
# carries the ellipsoid onto its authalic sphere with areas kept.
#
# Usage: Rscript paper/bench/bench_ellipsoid_area.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

A_WGS84 <- 6378137
F_WGS84 <- 1 / 298.257223563
TOLERANCE <- 1e-4
MAX_ARC <- 1e-3
CHUNK <- 1000L
GRIDS <- list(
  list(name = "ISEA3H", ap = 3, res = 5),
  list(name = "ISEA4H", ap = 4, res = 4),
  list(name = "ISEA7H", ap = 7, res = 3)
)
has_ellipsoid <- "ellipsoid" %in% names(formals(hex_grid))
MODELS <- if (has_ellipsoid) list(sphere = NULL, WGS84 = "WGS84") else list(sphere = NULL)

# Areas of the cells' rings in m^2 on the ellipsoid (a, f), traced a chunk
# of cells at a time
ring_area_m2 <- function(ids, g, a, f) {
  chunks <- split(seq_along(ids), (seq_along(ids) - 1L) %/% CHUNK)
  unlist(lapply(chunks, function(k) {
    rings <- hexify:::isea_cell_rings(ids[k], g@resolution, g@aperture,
                                      hexify:::icosa_arg(g), tolerance = TOLERANCE,
                                      max_arc = MAX_ARC)
    vapply(rings, function(r) {
      abs(geosphere::areaPolygon(unique(r[, 1:2, drop = FALSE]), a = a, f = f))
    }, numeric(1))
  }), use.names = FALSE)
}

rows <- list()
bands <- list()
for (spec in GRIDS) for (model in names(MODELS)) {
  g <- if (is.null(MODELS[[model]])) {
    hex_grid(resolution = spec$res, aperture = spec$ap)
  } else {
    hex_grid(resolution = spec$res, aperture = spec$ap, ellipsoid = MODELS[[model]])
  }
  ids <- hexify:::grid_cells(g)
  nominal_m2 <- unname(cell_area(ids, g)) * 1e6
  wgs84 <- ring_area_m2(ids, g, A_WGS84, F_WGS84) / nominal_m2
  sphere <- if (model == "sphere") {
    sphere_m2 <- 4 * pi * (hexify:::grid_radius_km(g) * 1000)^2
    ring_area_m2(ids, g, hexify:::grid_radius_km(g) * 1000, 0) /
      (nominal_m2 * sphere_m2 / (hexify:::body_surface_km2(hexify:::grid_radius_km(g)) * 1e6))
  } else {
    NA_real_
  }
  lat <- cell_to_lonlat(ids, g)$lat_deg

  rows[[length(rows) + 1]] <- data.frame(
    grid = spec$name, resolution = spec$res, earth_model = model, cells = length(ids),
    min_ratio = min(wgs84), max_ratio = max(wgs84),
    max_abs_deviation_pct = 100 * max(abs(wgs84 - 1)),
    lat_of_min = lat[which.min(wgs84)], lat_of_max = lat[which.max(wgs84)],
    sphere_check_max_abs_deviation = max(abs(sphere - 1)))
  print(rows[[length(rows)]], digits = 6)

  band <- cut(abs(lat), breaks = seq(0, 90, by = 10), include.lowest = TRUE)
  bands[[length(bands) + 1]] <- data.frame(
    grid = spec$name, resolution = spec$res, earth_model = model,
    abs_lat_band = levels(band),
    mean_ratio = as.numeric(tapply(wgs84, band, mean)),
    min_ratio = as.numeric(tapply(wgs84, band, min)),
    max_ratio = as.numeric(tapply(wgs84, band, max)))
}

extra <- c(sprintf("boundary tolerance: %g, max arc %g rad", TOLERANCE, MAX_ARC),
           paste("geosphere:", as.character(packageVersion("geosphere"))))
write_result(do.call(rbind, rows), "ellipsoid_area", extra = extra)
write_result(do.call(rbind, bands), "ellipsoid_area_by_latitude", extra = extra)
