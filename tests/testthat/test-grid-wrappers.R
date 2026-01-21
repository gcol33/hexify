
# tests/testthat/test-grid-wrappers.R
# Tests for grid-based coordinate conversion wrappers

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# GRID-BASED CELL CONVERSIONS
# =============================================================================

test_that("hexify_grid_to_cell works with grid object", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))

  expect_true(all(cell_ids > 0))
  expect_equal(length(cell_ids), 2)
})

test_that("hexify_grid_to_cell validates grid object", {
  expect_error(hexify_grid_to_cell("not a grid", lon = 0, lat = 0))
  expect_error(hexify_grid_to_cell(list(x = 1), lon = 0, lat = 0))
})

test_that("hexify_grid_cell_to_lonlat works with grid object", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))
  coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)

  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  expect_equal(length(coords$lon_deg), 2)
})

test_that("hexify_grid_cell_to_lonlat validates grid object", {
  expect_error(hexify_grid_cell_to_lonlat("not a grid", cell_id = 1))
})

test_that("grid-based round-trip is consistent", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  lon <- c(0, 45, -120)
  lat <- c(0, 30, -45)

  cell_ids <- hexify_grid_to_cell(grid, lon = lon, lat = lat)
  coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)
  cell_ids2 <- hexify_grid_to_cell(grid, lon = coords$lon_deg, lat = coords$lat_deg)

  expect_equal(cell_ids, cell_ids2)
})

# =============================================================================
# HIERARCHICAL INDEX CONVERSIONS
# =============================================================================

test_that("hexify_lonlat_to_h_index returns correct structure", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_lonlat_to_h_index(grid, lon = 0, lat = 45)

  expect_s3_class(result, "data.frame")
  expect_true("h_index" %in% names(result))
  expect_true("face" %in% names(result))
  expect_equal(nrow(result), 1)
})

test_that("hexify_lonlat_to_h_index validates grid object", {
  expect_error(
    hexify_lonlat_to_h_index("not a grid", lon = 0, lat = 0),
    "hexify_grid"
  )
})

test_that("hexify_lonlat_to_h_index validates input lengths", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  expect_error(
    hexify_lonlat_to_h_index(grid, lon = c(0, 1), lat = c(0)),
    "same length"
  )
})

test_that("hexify_lonlat_to_h_index validates numeric input", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  expect_error(
    hexify_lonlat_to_h_index(grid, lon = "0", lat = "45"),
    "numeric"
  )
})

test_that("hexify_lonlat_to_h_index warns on invalid coordinates", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  expect_warning(
    hexify_lonlat_to_h_index(grid, lon = 200, lat = 45),
    "longitude"
  )
  expect_warning(
    hexify_lonlat_to_h_index(grid, lon = 0, lat = 100),
    "latitude"
  )
})

test_that("hexify_lonlat_to_h_index handles NA values", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_lonlat_to_h_index(grid, lon = c(0, NA), lat = c(45, NA))

  expect_equal(nrow(result), 2)
  expect_true(is.na(result$h_index[2]))
  expect_true(is.na(result$face[2]))
})

test_that("hexify_h_index_to_lonlat returns correct structure", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  h_result <- hexify_lonlat_to_h_index(grid, lon = 0, lat = 45)
  coords <- hexify_h_index_to_lonlat(grid, h_result$h_index)

  expect_s3_class(coords, "data.frame")
  expect_true("lon" %in% names(coords))
  expect_true("lat" %in% names(coords))
})

test_that("hexify_h_index_to_lonlat validates grid object", {
  expect_error(
    hexify_h_index_to_lonlat("not a grid", h_index = "01"),
    "hexify_grid"
  )
})

test_that("hexify_h_index_to_lonlat validates character input", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  expect_error(
    hexify_h_index_to_lonlat(grid, h_index = 123),
    "character"
  )
})

test_that("hexify_h_index_to_lonlat handles NA values", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  h_result <- hexify_lonlat_to_h_index(grid, lon = 0, lat = 45)
  coords <- hexify_h_index_to_lonlat(grid, c(h_result$h_index, NA))

  expect_equal(nrow(coords), 2)
  expect_true(is.na(coords$lon[2]))
  expect_true(is.na(coords$lat[2]))
})

test_that("h_index round-trip is consistent", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  lon <- 10
  lat <- 50

  h_result <- hexify_lonlat_to_h_index(grid, lon = lon, lat = lat)
  coords <- hexify_h_index_to_lonlat(grid, h_result$h_index)
  h_result2 <- hexify_lonlat_to_h_index(grid, lon = coords$lon, lat = coords$lat)

  expect_equal(h_result$h_index, h_result2$h_index)
})

# =============================================================================
# ROUND-TRIP TEST FUNCTION
# =============================================================================

test_that("hexify_roundtrip_test returns correct structure", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_roundtrip_test(grid, lon = 10, lat = 50)

  expect_type(result, "list")
  expect_true("original" %in% names(result))
  expect_true("h_index" %in% names(result))
  expect_true("reconstructed" %in% names(result))
  expect_true("error" %in% names(result))
  expect_true("units" %in% names(result))
})

test_that("hexify_roundtrip_test default units are km", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_roundtrip_test(grid, lon = 10, lat = 50)

  expect_equal(result$units, "km")
})

test_that("hexify_roundtrip_test supports degrees units", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_roundtrip_test(grid, lon = 10, lat = 50, units = "degrees")

  expect_equal(result$units, "degrees")
})

test_that("hexify_roundtrip_test error is reasonable", {
  setup_icosa()

  grid <- hexify_grid(area = 1000, aperture = 3)
  result <- hexify_roundtrip_test(grid, lon = 10, lat = 50)

  # Error should be less than cell spacing (roughly sqrt(area))
  expect_lt(result$error, sqrt(1000))
})

# =============================================================================
# DIFFERENT APERTURES
# =============================================================================

test_that("grid wrappers work for all apertures", {
  skip_on_cran()  # Loop test across apertures
  setup_icosa()

  for (ap in c(3, 4, 7)) {
    grid <- hexify_grid(area = 10000, aperture = ap)

    cell_ids <- hexify_grid_to_cell(grid, lon = 0, lat = 45)
    coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)

    expect_true(cell_ids > 0)
    expect_true(is.finite(coords$lon_deg))
    expect_true(is.finite(coords$lat_deg))
  }
})

test_that("h_index works for different apertures", {
  skip_on_cran()  # Loop test across apertures
  setup_icosa()

  for (ap in c(3, 4, 7)) {
    grid <- hexify_grid(area = 10000, aperture = ap)

    result <- hexify_lonlat_to_h_index(grid, lon = 0, lat = 45)
    coords <- hexify_h_index_to_lonlat(grid, result$h_index)

    expect_true(nchar(result$h_index) > 0)
    expect_true(is.finite(coords$lon))
    expect_true(is.finite(coords$lat))
  }
})
