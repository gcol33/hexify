
# tests/testthat/test-constructors.R
# Tests for constructor functions

# =============================================================================
# hex_grid Tests
# =============================================================================

test_that("hex_grid creates HexGridInfo with area_km2", {
  grid <- hex_grid(area_km2 = 1000)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "3")
  expect_true(grid@resolution > 0)
})

test_that("hex_grid creates HexGridInfo with resolution", {
  grid <- hex_grid(resolution = 8, aperture = 3)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@resolution, 8L)
})

test_that("hex_grid handles aperture 4", {
  grid <- hex_grid(area_km2 = 1000, aperture = 4)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "4")
})

test_that("hex_grid handles aperture 7", {
  grid <- hex_grid(area_km2 = 10000, aperture = 7)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "7")
})

test_that("hex_grid handles mixed aperture 4/3", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "4/3")
})

test_that("hex_grid errors without area_km2 or resolution", {
  expect_error(hex_grid(), "Exactly one")
})

test_that("hex_grid errors with both area_km2 and resolution", {
  expect_error(hex_grid(area_km2 = 1000, resolution = 5), "not both")
})

test_that("hex_grid errors with invalid aperture", {
  expect_error(hex_grid(area_km2 = 1000, aperture = 5), "Aperture must be")
})

test_that("hex_grid errors with invalid area_km2", {
  expect_error(hex_grid(area_km2 = -100), "positive number")
  expect_error(hex_grid(area_km2 = "text"), "positive number")
})

test_that("hex_grid errors with invalid resolution", {
  expect_error(hex_grid(resolution = -1), "must be between")
  expect_error(hex_grid(resolution = 100), "must be between")
})

test_that("hex_grid works with resround = 'up'", {
  grid <- hex_grid(area_km2 = 1000, resround = "up")

  expect_s4_class(grid, "HexGridInfo")
})

test_that("hex_grid works with resround = 'down'", {
  grid <- hex_grid(area_km2 = 1000, resround = "down")

  expect_s4_class(grid, "HexGridInfo")
})

test_that("hex_grid errors with invalid resround", {
  expect_error(hex_grid(area_km2 = 1000, resround = "invalid"), "resround must be")
})

test_that("hex_grid sets custom CRS", {
  grid <- hex_grid(area_km2 = 1000, crs = 3857)

  expect_equal(grid@crs, 3857L)
})

test_that("hex_grid calculates area and diagonal", {
  grid <- hex_grid(area_km2 = 1000)

  expect_true(grid@area_km2 > 0)
  expect_true(grid@diagonal_km > 0)
})

# =============================================================================
# new_hex_data Tests (internal constructor)
# =============================================================================

test_that("new_hex_data creates HexData object", {
  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  cell_id <- c(12345, 67890)
  cell_center <- matrix(c(0, 10, 45, 50), ncol = 2, dimnames = list(NULL, c("lon", "lat")))

  result <- hexify:::new_hex_data(df, grid, cell_id, cell_center)

  expect_s4_class(result, "HexData")
  expect_equal(result@cell_id, c(12345, 67890))
})

test_that("new_hex_data errors on invalid data", {
  grid <- hex_grid(area_km2 = 1000)
  cell_id <- c(12345, 67890)
  cell_center <- matrix(c(0, 10, 45, 50), ncol = 2)

  expect_error(hexify:::new_hex_data(list(a = 1), grid, cell_id, cell_center),
               "must be a data.frame")
})

test_that("new_hex_data errors on invalid grid", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  cell_id <- c(12345, 67890)
  cell_center <- matrix(c(0, 10, 45, 50), ncol = 2)

  expect_error(hexify:::new_hex_data(df, list(a = 1), cell_id, cell_center),
               "must be a HexGridInfo")
})

test_that("new_hex_data adds column names to cell_center matrix", {
  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  cell_id <- c(12345, 67890)
  cell_center <- matrix(c(0, 10, 45, 50), ncol = 2)  # No column names

  result <- hexify:::new_hex_data(df, grid, cell_id, cell_center)

  expect_equal(colnames(result@cell_center), c("lon", "lat"))
})

test_that("new_hex_data works with sf data", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  cell_id <- c(12345, 67890)
  cell_center <- matrix(c(0, 10, 45, 50), ncol = 2, dimnames = list(NULL, c("lon", "lat")))

  result <- hexify:::new_hex_data(sf_data, grid, cell_id, cell_center)

  expect_s4_class(result, "HexData")
})
