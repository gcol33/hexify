
# tests/testthat/test-radius.R
# Tests for grids on bodies other than Earth

MARS_RADIUS_KM <- 3389.50
MARS_SURFACE_KM2 <- 4 * pi * MARS_RADIUS_KM^2

# =============================================================================
# RADIUS RESOLUTION
# =============================================================================

test_that("radius_km takes a number or a body name", {
  expect_equal(hexify:::resolve_radius_km(3389.5), 3389.5)
  expect_equal(hexify:::resolve_radius_km("mars"), MARS_RADIUS_KM)
  expect_equal(hexify:::resolve_radius_km("Mars"), MARS_RADIUS_KM)
  expect_equal(hexify:::resolve_radius_km(" mars "), MARS_RADIUS_KM)
  expect_equal(hexify:::resolve_radius_km("earth"), hexify:::EARTH_RADIUS_KM)
})

test_that("radius_km rejects invalid input", {
  expect_error(hexify:::resolve_radius_km(0), "positive")
  expect_error(hexify:::resolve_radius_km(-1), "positive")
  expect_error(hexify:::resolve_radius_km(NA), "positive")
  expect_error(hexify:::resolve_radius_km(Inf), "positive")
  expect_error(hexify:::resolve_radius_km(c(1000, 2000)), "positive")
  expect_error(hexify:::resolve_radius_km("pandora"), "Unknown body")
})

test_that("body surface area is the sphere area, Earth the ellipsoid area", {
  expect_equal(hexify:::body_surface_km2(MARS_RADIUS_KM), MARS_SURFACE_KM2)
  expect_equal(hexify:::body_surface_km2(hexify:::EARTH_RADIUS_KM),
               hexify:::EARTH_SURFACE_KM2)
})

test_that("km per degree scales with radius", {
  expect_equal(hexify:::km_per_degree(hexify:::EARTH_RADIUS_KM),
               hexify:::KM_PER_DEGREE)
  expect_equal(
    hexify:::km_per_degree(MARS_RADIUS_KM) / hexify:::km_per_degree(hexify:::EARTH_RADIUS_KM),
    MARS_RADIUS_KM / hexify:::EARTH_RADIUS_KM
  )
})

# =============================================================================
# GRID CONSTRUCTION
# =============================================================================

test_that("Earth stays the default", {
  g <- hex_grid(resolution = 6)
  expect_equal(hexify:::grid_radius_km(g), hexify:::EARTH_RADIUS_KM)
  expect_true(hexify:::is_earth_grid(g))

  named <- hex_grid(resolution = 6, radius_km = "earth")
  expect_equal(named@area_km2, g@area_km2)
  expect_equal(named@diagonal_km, g@diagonal_km)
})

test_that("a grid on another body scales its cell areas", {
  earth <- hex_grid(resolution = 6)
  mars <- hex_grid(resolution = 6, radius_km = "mars")

  expect_equal(hexify:::grid_radius_km(mars), MARS_RADIUS_KM)
  expect_equal(mars@area_km2, MARS_SURFACE_KM2 / (10 * 3^6 + 2))
  expect_equal(mars@area_km2 / earth@area_km2,
               MARS_SURFACE_KM2 / hexify:::EARTH_SURFACE_KM2)
  expect_equal(mars@diagonal_km, sqrt(mars@area_km2 * 2 / sqrt(3)))
})

test_that("a body name and its radius give the same grid", {
  by_name <- hex_grid(resolution = 4, radius_km = "mars")
  by_number <- hex_grid(resolution = 4, radius_km = MARS_RADIUS_KM)
  expect_equal(by_name@area_km2, by_number@area_km2)
  expect_equal(by_name@resolution, by_number@resolution)
})

test_that("resolution for a target area follows the body", {
  mars <- hex_grid(area_km2 = 1000, radius_km = "mars")
  n_cells <- 10 * 3^mars@resolution + 2

  # The chosen resolution is the one whose cells are closest to 1000 km^2 on
  # Mars, so a smaller body reaches that cell size at a coarser resolution.
  expect_equal(mars@area_km2, MARS_SURFACE_KM2 / n_cells)
  expect_lt(mars@resolution, hex_grid(area_km2 = 1000)@resolution)
  expect_lt(abs(log(mars@area_km2 / 1000)), log(3) / 2 + 1e-8)
})

test_that("resolution for a target area follows the body for mixed apertures", {
  mars <- hex_grid(area_km2 = 5000, aperture = "4/3", radius_km = "mars")
  n_cells <- hexify:::aperture_n_cells(mars@aperture, mars@resolution)
  expect_equal(mars@area_km2, MARS_SURFACE_KM2 / n_cells)
  expect_lt(abs(log(mars@area_km2 / 5000)), log(4) / 2 + 1e-8)
})

test_that("radius_km must be positive in a grid object", {
  expect_error(hex_grid(resolution = 4, radius_km = -100), "positive")
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = 4L, area_km2 = 100,
        diagonal_km = 10, crs = 4326L, grid_type = "isea", radius_km = -1),
    "radius_km must be positive"
  )
})

# =============================================================================
# CELL GEOMETRY IS ANGULAR, SO RADIUS-INDEPENDENT
# =============================================================================

test_that("cells and centres are the same on every body", {
  earth <- hex_grid(resolution = 5)
  mars <- hex_grid(resolution = 5, radius_km = "mars")

  lon <- c(0, 16.37, -70.5, 120)
  lat <- c(0, 48.21, -33.4, 65)

  expect_equal(lonlat_to_cell(lon, lat, mars), lonlat_to_cell(lon, lat, earth))
  expect_equal(cell_to_lonlat(lonlat_to_cell(lon, lat, mars), mars),
               cell_to_lonlat(lonlat_to_cell(lon, lat, earth), earth))
  expect_equal(get_neighbors(1000, mars), get_neighbors(1000, earth))
})

test_that("cell_area reports the body's cell area", {
  mars <- hex_grid(resolution = 5, radius_km = "mars")
  cells <- lonlat_to_cell(c(0, 30), c(0, 40), mars)
  expect_equal(unname(cell_area(cells, mars)), rep(mars@area_km2, 2))
})

# =============================================================================
# STATISTICS
# =============================================================================

test_that("dgearthstat reports the body the grid covers", {
  mars <- hex_grid(resolution = 6, radius_km = "mars")
  stats <- dgearthstat(mars)

  expect_equal(stats$area_km, MARS_SURFACE_KM2)
  expect_equal(stats$n_cells, 10 * 3^6 + 2)
  expect_equal(stats$cell_area_km2, MARS_SURFACE_KM2 / (10 * 3^6 + 2))
  expect_equal(stats$cell_spacing_km, sqrt(2 * stats$cell_area_km2 / sqrt(3)))
  expect_equal(stats$cls_km, 2 * sqrt(stats$cell_area_km2 / pi))
})

test_that("hexify_compare_resolutions takes a radius", {
  cmp <- hexify_compare_resolutions(aperture = 3, res_range = 0:5,
                                    radius_km = "mars")
  expect_equal(cmp$cell_area_km2, MARS_SURFACE_KM2 / cmp$n_cells)

  earth <- hexify_compare_resolutions(aperture = 3, res_range = 0:5)
  expect_equal(cmp$n_cells, earth$n_cells)
  expect_true(all(cmp$cell_area_km2 < earth$cell_area_km2))
})

# =============================================================================
# LEGACY GRID OBJECTS
# =============================================================================

test_that("legacy hexify_grid carries a radius", {
  g <- hexify_grid(area = 1000, aperture = 3, radius_km = "mars")
  expect_equal(g$radius_km, MARS_RADIUS_KM)

  stats <- dgearthstat(g)
  expect_equal(stats$area_km, MARS_SURFACE_KM2)
  expect_equal(stats$cell_area_km2, MARS_SURFACE_KM2 / stats$n_cells)
})

test_that("the radius survives conversion between grid representations", {
  s4 <- hex_grid(resolution = 5, radius_km = "mars")
  legacy <- hexify:::HexGridInfo_to_hexify_grid(s4)
  expect_equal(legacy$radius_km, MARS_RADIUS_KM)

  back <- hexify:::hexify_grid_to_HexGridInfo(legacy)
  expect_equal(hexify:::grid_radius_km(back), MARS_RADIUS_KM)
  expect_equal(back@area_km2, s4@area_km2)
})

test_that("a grid with no radius reads as Earth", {
  expect_equal(hexify:::grid_radius_km(list(aperture = 3, resolution = 5)),
               hexify:::EARTH_RADIUS_KM)
  expect_equal(hexify:::grid_radius_km(list(radius_km = NA_real_)),
               hexify:::EARTH_RADIUS_KM)
})

# =============================================================================
# HEXIFY AND SF EXPORT
# =============================================================================

test_that("hexify() passes a radius through to the grid", {
  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  hd <- hexify(df, lon = "lon", lat = "lat", resolution = 5, radius_km = "mars")

  expect_equal(hexify:::grid_radius_km(hd@grid), MARS_RADIUS_KM)
  expect_equal(as.data.frame(hd)$cell_area_km2,
               rep(hd@grid@area_km2, 3))

  earth <- hexify(df, lon = "lon", lat = "lat", resolution = 5)
  expect_equal(hd@cell_id, earth@cell_id)
})

test_that("polygon export samples the same cells on every body", {
  skip_if_not_installed("sf")

  earth <- hex_grid(resolution = 3)
  mars <- hex_grid(resolution = 3, radius_km = "mars")

  bbox <- c(-10, 35, 30, 60)
  expect_equal(sort(grid_rect(bbox, mars)$cell_id),
               sort(grid_rect(bbox, earth)$cell_id))
})

# =============================================================================
# H3 IS EARTH ONLY
# =============================================================================

test_that("H3 grids reject a non-Earth radius", {
  expect_error(hex_grid(resolution = 5, type = "h3", radius_km = "mars"),
               "Earth")
  expect_error(hexify_compare_resolutions(type = "h3", radius_km = "mars"),
               "Earth")
  expect_error(
    new("HexGridInfo", aperture = "7", resolution = 5L, area_km2 = 100,
        diagonal_km = 10, crs = 4326L, grid_type = "h3",
        radius_km = MARS_RADIUS_KM),
    "Earth"
  )
})

test_that("h3_crosswalk rejects a non-Earth ISEA grid", {
  mars <- hex_grid(resolution = 5, radius_km = "mars")
  cells <- lonlat_to_cell(c(0, 10), c(45, 50), mars)
  expect_error(h3_crosswalk(cells, mars, direction = "isea_to_h3"), "Earth")
})
