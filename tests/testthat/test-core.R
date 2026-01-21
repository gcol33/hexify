# test-core.R
# Fast tests using precomputed cached data
# All tests must complete in < 5 seconds total

# Load cached test data
cache <- readRDS(test_path("fixtures/test_cache.rds"))

# =============================================================================
# HEXIFY FUNCTION
# =============================================================================

test_that("hexify returns correct cell IDs", {
  df <- data.frame(lon = c(0, 10, 16.37, 2.35, -3.70), lat = c(45, 50, 48.21, 48.86, 40.42))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_s4_class(result, "HexData")
  expect_equal(result@cell_id, cache$hexify_basic$cell_id)
  expect_equal(result@grid@resolution, cache$hexify_basic$resolution)
})

# =============================================================================
# HEX_GRID
# =============================================================================

test_that("hex_grid returns correct parameters", {
  grid <- hex_grid(area_km2 = 1000)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@resolution, cache$grid_1000$resolution)
  expect_equal(grid@aperture, cache$grid_1000$aperture)
})

# =============================================================================
# CELL CONVERSIONS - APERTURE 3
# =============================================================================

test_that("cell conversions work for aperture 3", {
  cells <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 5, aperture = 3)
  expect_equal(cells, cache$cells_ap3$cells)

  coords <- hexify_cell_to_lonlat(cells, resolution = 5, aperture = 3)
  expect_equal(coords$lon_deg, cache$cells_ap3$lon, tolerance = 1e-6)
  expect_equal(coords$lat_deg, cache$cells_ap3$lat, tolerance = 1e-6)
})

# =============================================================================
# CELL CONVERSIONS - APERTURE 4
# =============================================================================

test_that("cell conversions work for aperture 4", {
  cells <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 5, aperture = 4)
  expect_equal(cells, cache$cells_ap4$cells)

  coords <- hexify_cell_to_lonlat(cells, resolution = 5, aperture = 4)
  expect_equal(coords$lon_deg, cache$cells_ap4$lon, tolerance = 1e-6)
  expect_equal(coords$lat_deg, cache$cells_ap4$lat, tolerance = 1e-6)
})

# =============================================================================
# CELL CONVERSIONS - APERTURE 7
# =============================================================================

test_that("cell conversions work for aperture 7", {
  cells <- hexify_lonlat_to_cell(c(0, 10, -30), c(45, 50, -15), resolution = 3, aperture = 7)
  expect_equal(cells, cache$cells_ap7$cells)

  coords <- hexify_cell_to_lonlat(cells, resolution = 3, aperture = 7)
  expect_equal(coords$lon_deg, cache$cells_ap7$lon, tolerance = 1e-6)
  expect_equal(coords$lat_deg, cache$cells_ap7$lat, tolerance = 1e-6)
})

# =============================================================================
# PROJECTION
# =============================================================================

test_that("snyder forward projection works", {
  proj <- cpp_snyder_forward(10, 45)

  expect_equal(proj[["face"]], cache$proj_forward$face)
  expect_equal(proj[["icosa_triangle_x"]], cache$proj_forward$icosa_triangle_x, tolerance = 1e-9)
  expect_equal(proj[["icosa_triangle_y"]], cache$proj_forward$icosa_triangle_y, tolerance = 1e-9)
})

# =============================================================================
# QUAD IJ
# =============================================================================

test_that("quad ij conversion works", {
  qij <- hexify_cell_to_quad_ij(c(100, 1000, 5000), resolution = 10, aperture = 3)

  expect_equal(qij$quad, cache$quad_ij$quad)
  expect_equal(qij$i, cache$quad_ij$i)
  expect_equal(qij$j, cache$quad_ij$j)
})

# =============================================================================
# INDEX CONVERSION
# =============================================================================

test_that("z3 index encoding works", {
  idx <- hexify_cell_to_index(0, 5, 4, 3, 3, "z3")
  expect_equal(idx, cache$index_z3_sample)
})

# =============================================================================
# FACE DETECTION
# =============================================================================

test_that("face centers are valid", {
  centers <- hexify_face_centers()

  expect_equal(nrow(centers), 20)
  expect_true(all(centers$lon >= -180 & centers$lon <= 180))
  expect_true(all(centers$lat >= -90 & centers$lat <= 90))
})
