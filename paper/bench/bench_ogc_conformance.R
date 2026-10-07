# Two measurements behind the OGC Topic 21 conformance table.
#
# 1. Centroid (Req 27): the EAERS asks that a zone's position be its centroid,
#    "computed as the geodesic centre of surface area". cell_to_lonlat() gives
#    the cell centre, the inverse projection of the centre of the planar
#    hexagon. The area centroid of a region R of the unit sphere points along
#    the integral of x over R, which by Stokes' theorem is half the integral of
#    x cross dx around the boundary; along a great-circle piece from a to b
#    that is the arc length times the unit normal a x b / |a x b|. Each wall is
#    followed with hexify's own boundary densification, its pieces read as
#    great-circle arcs. The offset is the angle between the two directions,
#    in km and as a fraction of the distance between neighbouring centres.
#
# 2. Area on WGS84 (Req 28/29): hexify reads geodetic latitude as latitude on
#    a sphere of the ellipsoid's area. A small region at geodetic latitude phi
#    then has ellipsoidal area M N / R_q^2 times its spherical area (M, N the
#    meridional and prime-vertical radii of curvature, R_q the authalic
#    radius), so a cell's area on WGS84 departs from the sphere's by that
#    factor, averaged over the cell. Reported at each cell's area centroid,
#    for every cell of a grid.
#
# Usage: Rscript paper/bench/bench_ogc_conformance.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

unit <- function(lon, lat) {
  r <- pi / 180
  cbind(cos(lat * r) * cos(lon * r), cos(lat * r) * sin(lon * r), sin(lat * r))
}
cross <- function(a, b) {
  cbind(a[, 2] * b[, 3] - a[, 3] * b[, 2], a[, 3] * b[, 1] - a[, 1] * b[, 3],
        a[, 1] * b[, 2] - a[, 2] * b[, 1])
}

# Area centroid direction of each cell from its densified boundary
area_centroids <- function(ids, g) {
  rings <- hexify:::isea_cell_rings(ids, g@resolution, g@aperture, hexify:::icosa_arg(g),
                                    tolerance = 1e-5)
  ctr <- cell_to_lonlat(ids, g)
  mid <- unit(ctr$lon_deg, ctr$lat_deg)
  cen <- t(vapply(rings, function(r) {
    a <- unit(r[, 1], r[, 2])
    b <- a[c(2:nrow(a), 1), , drop = FALSE]
    n <- cross(a, b)
    s <- sqrt(rowSums(n^2))
    theta <- atan2(s, rowSums(a * b))
    use <- s > 0
    v <- colSums(n[use, , drop = FALSE] / s[use] * theta[use]) / 2
    v / sqrt(sum(v^2))
  }, numeric(3)))
  # A ring running clockwise gives the antipode
  cen * sign(rowSums(cen * mid))
}

angle <- function(a, b) atan2(sqrt(rowSums(cross(a, b)^2)), rowSums(a * b))

GRIDS <- list(
  list(name = "ISEA3H", ap = 3, res = c(3, 6, 9, 12)),
  list(name = "ISEA4H", ap = 4, res = c(2, 4, 6, 8)),
  list(name = "ISEA7H", ap = 7, res = c(2, 3, 4, 5, 7)),
  list(name = "FULLER7H", ap = 7, res = c(3, 5), projection = "fuller")
)
N_SAMPLE <- 20000

rows <- list()
for (spec in GRIDS) {
  for (res in spec$res) {
    g <- hex_grid(resolution = res, aperture = spec$ap,
                  projection = if (is.null(spec$projection)) "isea" else spec$projection)
    n <- hexify:::grid_n_cells(g)
    ids <- if (n <= N_SAMPLE) hexify:::grid_cells(g) else {
      p <- sphere_points(N_SAMPLE, seed = res)
      unique(lonlat_to_cell(p$lon, p$lat, g))
    }
    ctr <- cell_to_lonlat(ids, g)
    cen <- area_centroids(ids, g)
    off <- angle(cen, unit(ctr$lon_deg, ctr$lat_deg))
    radius <- hexify:::grid_radius_km(g)
    spacing <- mean(wall_metrics(ids[seq_len(min(length(ids), 2000))], g)$center_distance_km)
    rows[[length(rows) + 1]] <- data.frame(
      grid = spec$name, resolution = res, cells = length(ids),
      all_cells = length(ids) == n,
      offset_km_median = median(off) * radius, offset_km_max = max(off) * radius,
      offset_over_spacing_max = max(off) * radius / spacing,
      spacing_km = spacing)
    print(rows[[length(rows)]], digits = 4)
  }
}
write_result(do.call(rbind, rows), "centroid_offset",
             extra = sprintf("cells sampled from %d uniform points where a grid has more", N_SAMPLE))

# Area on WGS84 relative to the sphere of equal area
a <- 6378137; f <- 1 / 298.257223563; e2 <- f * (2 - f); e <- sqrt(e2)
q <- function(phi) (1 - e2) * (sin(phi) / (1 - e2 * sin(phi)^2) -
                                 log((1 - e * sin(phi)) / (1 + e * sin(phi))) / (2 * e))
rq2 <- a^2 * q(pi / 2) / 2
scale_at <- function(lat) {
  s2 <- sin(lat * pi / 180)^2
  (a * (1 - e2) / (1 - e2 * s2)^1.5) * (a / sqrt(1 - e2 * s2)) / rq2
}
area_rows <- list()
for (spec in list(list("ISEA3H", 3, 8), list("ISEA4H", 4, 6), list("ISEA7H", 7, 5))) {
  g <- hex_grid(resolution = spec[[3]], aperture = spec[[2]])
  ids <- hexify:::grid_cells(g)
  cen <- area_centroids(ids, g)
  lat <- asin(pmax(-1, pmin(1, cen[, 3]))) * 180 / pi
  s <- scale_at(lat)
  area_rows[[length(area_rows) + 1]] <- data.frame(
    grid = spec[[1]], resolution = spec[[3]], cells = length(ids),
    min_rel_area = min(s), max_rel_area = max(s),
    max_abs_deviation_pct = 100 * max(abs(s - 1)),
    equator = scale_at(0), lat45 = scale_at(45), pole = scale_at(90))
  print(area_rows[[length(area_rows)]], digits = 6)
}
write_result(do.call(rbind, area_rows), "wgs84_area_scale")
