# Grids on an ellipsoid of revolution: geodetic latitude carried to authalic
# latitude before the projection and back after it. The geodesic areas of
# the cells on WGS84 are measured in paper/bench/bench_ellipsoid_area.R.

wgs84_f <- 1 / 298.257223563
wgs84_icosa <- c(standard_icosa("icosahedron"), wgs84_f)

# The ellipsoid's area (a = 1) between the equator and latitude phi, per
# radian of longitude: F(phi) = q(phi) / 2, written out here independently of
# the C++ code.
ellipsoid_F <- function(phi, f) {
  e2 <- f * (2 - f)
  e <- sqrt(e2)
  s <- sin(phi)
  (1 - e2) / 2 * (s / (1 - e2 * s^2) + atanh(e * s) / e)
}

# Area of a lon/lat ring on the ellipsoid (a = 1) by Green's theorem in
# (lon, lat): the integral of F(lat) d(lon) around it. A ring around a pole
# takes the area beyond it towards that pole; a ring over a pole, of a cell
# with the pole on an edge, which it draws as two meridians, is not measured
# (NA).
ring_area_ellipsoid <- function(ring, f, centre_lat) {
  lon <- ring[, 1] * pi / 180
  lat <- ring[, 2] * pi / 180
  dlon <- diff(lon)
  dlon <- (dlon + pi) %% (2 * pi) - pi
  if (any(abs(dlon) > 0.99 * pi)) return(NA_real_)
  Fm <- (ellipsoid_F(lat[-1], f) + ellipsoid_F(lat[-length(lat)], f)) / 2
  pole <- ellipsoid_F(pi / 2, f) * sign(centre_lat)
  abs(sum((Fm - pole) * dlon))
}

test_that("latitude conversions are inverse to rounding on every ellipsoid", {
  lat <- seq(-90, 90, by = 0.01)
  for (e in list(ELLIPSOIDS$wgs84, ELLIPSOIDS$mars, ELLIPSOIDS$saturn, c(a_km = 1, f = 0.25))) {
    icosa <- c(standard_icosa("icosahedron"), e[["f"]])
    beta <- sphere_lat(lat, icosa)
    expect_lt(max(abs(geodetic_lat(beta, icosa) - lat)), 1e-13)
    # the series against the closed form it is fitted to
    exact <- cpp_sphere_latitude(icosa, lat, exact = TRUE)
    expect_lt(max(abs(beta - exact)), 1e-13)
    expect_lt(max(abs(cpp_sphere_latitude(icosa, exact, inverse = TRUE, exact = TRUE) - lat)),
              1e-13)
  }
})

test_that("the authalic latitude matches Karney's series and OGC's figures on WGS84", {
  # Karney (2023), On auxiliary latitudes, coefficients C_xi_phi (A19) to
  # sixth order in the third flattening n, as DGGAL's authalic.ec lists them
  n <- wgs84_f / (2 - wgs84_f)
  C <- list(c(-4/3, -4/45, 88/315, 538/4725, 20824/467775, -44732/2837835),
            c(34/45, 8/105, -2482/14175, -37192/467775, -12467764/212837625),
            c(-1532/2835, -898/14175, 54968/467775, 100320856/1915538625),
            c(6007/14175, 24496/467775, -5884124/70945875),
            c(-23356/66825, -839792/19348875),
            570284222/1915538625)
  cp <- vapply(1:6, function(k) sum(C[[k]] * n^(seq_along(C[[k]]) - 1)) * n^k, numeric(1))
  phi <- seq(-90, 90, by = 0.05) * pi / 180
  karney <- phi + vapply(phi, function(p) sum(cp * sin(2 * (1:6) * p)), numeric(1))
  expect_lt(max(abs(sphere_lat(phi * 180 / pi, wgs84_icosa) - karney * 180 / pi)), 1e-13)

  # OGC's ISEA3H and ISEA7H definitions put vertex 0 at authalic latitude
  # arctan(golden ratio), "~58.397145907431 N geodetic latitude"
  vertex <- geodetic_lat(atan((1 + sqrt(5)) / 2) * 180 / pi, wgs84_icosa)
  expect_equal(vertex, 58.397145907431, tolerance = 1e-12 / 58)
  # DGGRID's WGS84 authalic radius, WGS84_AUTHALIC_RADIUS_KM
  expect_equal(authalic_radius_km(ELLIPSOIDS$wgs84), 6371.007180918475, tolerance = 1e-14)
})

test_that("hex_grid() takes an ellipsoid by name or by its axes", {
  g <- hex_grid(resolution = 4, aperture = 3, ellipsoid = "WGS84")
  expect_identical(grid_ellipsoid(g), ELLIPSOIDS$wgs84)
  expect_identical(icosa_arg(g)[6], wgs84_f)
  expect_equal(grid_radius_km(g), authalic_radius_km(ELLIPSOIDS$wgs84))
  expect_identical(g@crs, 4326L)
  expect_true(is_earth_grid(g))
  expect_identical(summary(g)$ellipsoid, ELLIPSOIDS$wgs84)
  expect_output(print(g), "Ellipsoid:   WGS84, authalic radius 6371.0072 km")

  expect_identical(grid_ellipsoid(hex_grid(resolution = 4, ellipsoid = "earth")),
                   ELLIPSOIDS$wgs84)
  axes <- hex_grid(resolution = 4, ellipsoid = c(6378.137, 1 / 298.257222101))
  expect_identical(grid_ellipsoid(axes), ELLIPSOIDS$grs80)
  expect_output(print(axes), "GRS80")
  expect_identical(axes@crs, "+proj=longlat +ellps=GRS80 +no_defs")

  mars <- hex_grid(resolution = 3, aperture = 7, ellipsoid = "mars")
  expect_false(is_earth_grid(mars))
  expect_equal(as.numeric(sf::st_crs(mars)$SemiMajor), 3396190)
  expect_equal(sf::st_crs(mars)$InvFlattening, 3396.19 / (3396.19 - 3376.20))
  expect_equal(grid_radius_km(mars), authalic_radius_km(ELLIPSOIDS$mars))

  # no ellipsoid is today's grid
  g0 <- hex_grid(resolution = 4, aperture = 3)
  expect_identical(hex_grid(resolution = 4, aperture = 3, ellipsoid = NULL), g0)
  expect_length(icosa_arg(g0), 5L)
  expect_length(grid_ellipsoid(g0), 0L)
  # f = 0 is the sphere of radius a
  sphere <- hex_grid(resolution = 4, ellipsoid = c(1000, 0))
  expect_length(icosa_arg(sphere), 5L)
  expect_equal(grid_radius_km(sphere), 1000)
})

test_that("hex_grid() refuses an ellipsoid it cannot use", {
  expect_error(hex_grid(resolution = 3, type = "h3", ellipsoid = "WGS84"), "H3 reads latitude")
  expect_error(hex_grid(resolution = 3, ellipsoid = "pluto"), "Unknown ellipsoid")
  expect_error(hex_grid(resolution = 3, ellipsoid = c(6378, 1)), "0 <= f < 1")
  expect_error(hex_grid(resolution = 3, ellipsoid = c(-1, 0.1)), "a > 0")
  expect_error(hex_grid(resolution = 3, ellipsoid = "WGS84", radius_km = "mars"),
               "ellipsoid sets the radius")
  expect_s4_class(hex_grid(resolution = 3, ellipsoid = "WGS84",
                           radius_km = authalic_radius_km(ELLIPSOIDS$wgs84)), "HexGridInfo")
  g <- hex_grid(resolution = 3, ellipsoid = "WGS84")
  expect_error(as_dggrid(g), "DGGRID reads latitude on the sphere")
})

test_that("an ellipsoid moves points by the authalic latitude and keeps cell IDs", {
  g0 <- hex_grid(resolution = 5, aperture = 3)
  g <- hex_grid(resolution = 5, aperture = 3, ellipsoid = "WGS84")
  ids <- grid_cells(g)
  c0 <- cell_to_lonlat(ids, g0)
  c1 <- cell_to_lonlat(ids, g)
  expect_identical(c1$lon_deg, c0$lon_deg)
  expect_equal(c1$lat_deg, geodetic_lat(c0$lat_deg, wgs84_icosa), tolerance = 1e-13)
  expect_identical(lonlat_to_cell(c1$lon_deg, c1$lat_deg, g), ids)
  expect_lt(max(abs(c1$lat_deg - c0$lat_deg)), 0.13)
  expect_gt(max(abs(c1$lat_deg - c0$lat_deg)), 0.12)
  # a point is read where its authalic latitude puts it on the default grid
  set.seed(4)
  lon <- runif(500, -180, 180)
  lat <- asin(runif(500, -1, 1)) * 180 / pi
  expect_identical(lonlat_to_cell(lon, lat, g),
                   lonlat_to_cell(lon, sphere_lat(lat, wgs84_icosa), g0))
  # the hierarchy and neighbours do not depend on the earth model
  expect_identical(get_parent(ids[1:50], g), get_parent(ids[1:50], g0))
  expect_identical(get_neighbors(ids[1:20], g), get_neighbors(ids[1:20], g0))
})

test_that("cells are equal-area on the ellipsoid and add up to its area", {
  for (spec in list(list(ap = 3, res = 3), list(ap = 7, res = 2))) {
    for (proj in c("isea", "ivea")) {
      g <- hex_grid(resolution = spec$res, aperture = spec$ap, projection = proj,
                    ellipsoid = "WGS84")
      ids <- grid_cells(g)
      area <- cell_area(ids, g)
      # the ellipsoid's surface, 2 pi a^2 (1 + (1 - e^2) atanh(e) / e)
      e2 <- wgs84_f * (2 - wgs84_f)
      surface <- 2 * pi * 6378.137^2 * (1 + (1 - e2) * atanh(sqrt(e2)) / sqrt(e2))
      expect_equal(sum(area), surface, tolerance = 1e-13)
      expect_equal(sum(area), EARTH_SURFACE_KM2, tolerance = 1e-12)

      rings <- isea_cell_rings(ids, g@resolution, g@aperture, icosa_arg(g),
                               tolerance = 1e-5, max_arc = 2e-3)
      lat <- cell_to_lonlat(ids, g)$lat_deg
      on_ellipsoid <- vapply(seq_along(rings), function(i) {
        ring_area_ellipsoid(rings[[i]], wgs84_f, lat[i])
      }, numeric(1)) * 6378.137^2
      expect_lte(sum(is.na(on_ellipsoid)), 4L)
      expect_lt(max(abs(on_ellipsoid / area - 1), na.rm = TRUE), 1e-7)
    }
  }
  # without the ellipsoid the same measure departs by up to 0.9 percent
  g0 <- hex_grid(resolution = 3, aperture = 3)
  ids <- grid_cells(g0)
  rings <- isea_cell_rings(ids, 3L, "3", icosa_arg(g0), tolerance = 1e-5, max_arc = 2e-3)
  lat <- cell_to_lonlat(ids, g0)$lat_deg
  ratio <- vapply(seq_along(rings), function(i) {
    ring_area_ellipsoid(rings[[i]], wgs84_f, lat[i])
  }, numeric(1)) * 6378.137^2 / cell_area(ids, g0)
  expect_gt(max(ratio, na.rm = TRUE), 1.008)
  expect_lt(min(ratio, na.rm = TRUE), 0.996)
})

test_that("Fuller cells on an ellipsoid report their area on it", {
  g <- hex_grid(resolution = 3, aperture = 3, projection = "fuller", ellipsoid = "WGS84")
  ids <- grid_cells(g)
  area <- cell_area(ids, g)
  expect_equal(sum(area), EARTH_SURFACE_KM2, tolerance = 1e-9)
  rings <- isea_cell_rings(ids, 3L, "3", icosa_arg(g), tolerance = 1e-5, max_arc = 2e-3)
  lat <- cell_to_lonlat(ids, g)$lat_deg
  on_ellipsoid <- vapply(seq_along(rings), function(i) {
    ring_area_ellipsoid(rings[[i]], wgs84_f, lat[i])
  }, numeric(1)) * 6378.137^2
  expect_lt(max(abs(on_ellipsoid / area - 1), na.rm = TRUE), 1e-6)
})

test_that("the indicatrix on an ellipsoid keeps area under the equal-area projections", {
  set.seed(9)
  lon <- runif(300, -180, 180)
  lat <- asin(runif(300, -1, 1)) * 180 / pi
  for (proj in c("isea", "ivea")) {
    g <- hex_grid(resolution = 2, projection = proj, ellipsoid = "WGS84")
    d <- projection_distortion(g, lon, lat)
    expect_lt(max(abs(d$areal - 1)), 1e-9)
    d0 <- projection_distortion(hex_grid(resolution = 2, projection = proj), lon,
                                sphere_lat(lat, wgs84_icosa))
    # the authalic map adds at most its own distortion, about e^2 / 2
    expect_lt(max(abs(d$angular - d0$angular)), 0.4)
    expect_gt(max(abs(d$angular - d0$angular)), 0.01)
  }
})

test_that("a region orientation is placed from the region's geodetic centre", {
  g <- hex_grid(resolution = 4, orientation = "region", region = c(10, 46.5),
                ellipsoid = "WGS84")
  expect_identical(g@orientation,
                   region_orientation(10, sphere_lat(46.5, wgs84_icosa)))
  o <- hex_grid(resolution = 4, orientation = "ogc", ellipsoid = "WGS84")
  expect_identical(unname(o@orientation), c(11.20, 58.282525588538995, 0))
  expect_error(hex_grid(resolution = 3, orientation = "ogc", polyhedron = "octahedron",
                        aperture = 4), "places an icosahedron")
})

test_that("samples, links and polygons of an ellipsoid grid stay in their cells", {
  g <- hex_grid(resolution = 3, aperture = 3, ellipsoid = "WGS84")
  ids <- lonlat_to_cell(c(0, 120, -60), c(60, -45, 80), g)
  set.seed(2)
  pts <- hex_sample(ids, g, n = 50)
  expect_identical(lonlat_to_cell(pts$lon, pts$lat, g), pts$cell_id)
  polys <- cell_to_sf(ids, g)
  ctr <- cell_to_lonlat(ids, g)
  inside <- sf::st_within(sf::st_as_sf(ctr, coords = c("lon_deg", "lat_deg"), crs = 4326),
                          sf::st_geometry(polys), sparse = FALSE)
  expect_true(all(diag(inside)))
  expect_identical(grid_ellipsoid(grid_at_resolution(g, 2L)), ELLIPSOIDS$wgs84)
})

test_that("a grid saved before grids carried an ellipsoid reads the sphere", {
  g <- hex_grid(resolution = 5)
  attr(g, "ellipsoid") <- NULL
  expect_false(.hasSlot(g, "ellipsoid"))
  expect_length(grid_ellipsoid(g), 0L)
  expect_identical(extract_grid(g)@ellipsoid, numeric(0))
  expect_identical(lonlat_to_cell(10, 45, g), lonlat_to_cell(10, 45, hex_grid(resolution = 5)))
})
