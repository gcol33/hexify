# tests/testthat/test-constants.R
# Tests for constants.R validation functions and helpers

# =============================================================================
# VALIDATE_LON
# =============================================================================

test_that("validate_lon returns TRUE for valid values", {
  result <- hexify:::validate_lon(c(0, 45, -120, 180, -180), warn = FALSE)
  expect_true(all(result))
})

test_that("validate_lon returns FALSE for out-of-range values", {
  result <- hexify:::validate_lon(c(0, 200, -200), warn = FALSE)
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("validate_lon handles NA values", {
  result <- hexify:::validate_lon(c(0, NA, 45), warn = FALSE)
  expect_equal(result, c(TRUE, TRUE, TRUE))
})

test_that("validate_lon warns on invalid values", {
  expect_warning(
    hexify:::validate_lon(c(200)),
    "outside valid range"
  )
})

test_that("validate_lon errors on non-numeric input", {
  expect_error(
    hexify:::validate_lon("not numeric"),
    "must be numeric"
  )
})

# =============================================================================
# VALIDATE_LAT
# =============================================================================

test_that("validate_lat returns TRUE for valid values", {
  result <- hexify:::validate_lat(c(0, 45, -45, 90, -90), warn = FALSE)
  expect_true(all(result))
})

test_that("validate_lat returns FALSE for out-of-range values", {
  result <- hexify:::validate_lat(c(0, 100, -100), warn = FALSE)
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("validate_lat handles NA values", {
  result <- hexify:::validate_lat(c(0, NA, 45), warn = FALSE)
  expect_equal(result, c(TRUE, TRUE, TRUE))
})

test_that("validate_lat warns on invalid values", {
  expect_warning(
    hexify:::validate_lat(c(100)),
    "outside valid range"
  )
})

test_that("validate_lat errors on non-numeric input", {
  expect_error(
    hexify:::validate_lat("not numeric"),
    "must be numeric"
  )
})

# =============================================================================
# VALIDATE_RESOLUTION
# =============================================================================

test_that("validate_resolution accepts valid resolutions", {
  for (res in c(0, 1, 10, 20, 30)) {
    expect_true(hexify:::validate_resolution(res))
  }
})

test_that("validate_resolution errors on negative resolution", {
  expect_error(
    hexify:::validate_resolution(-1),
    "between"
  )
})

test_that("validate_resolution errors on resolution > 30", {
  expect_error(
    hexify:::validate_resolution(31),
    "between"
  )
})

test_that("validate_resolution errors on non-numeric input", {
  expect_error(
    hexify:::validate_resolution("10"),
    "must be a single numeric"
  )
})

test_that("validate_resolution errors on vector input", {
  expect_error(
    hexify:::validate_resolution(c(1, 2)),
    "must be a single numeric"
  )
})

test_that("validate_resolution errors on NA input", {
  expect_error(
    hexify:::validate_resolution(NA),
    "single numeric"
  )
})

# =============================================================================
# VALIDATE_APERTURE
# =============================================================================

test_that("validate_aperture accepts valid apertures", {
  expect_true(hexify:::validate_aperture(3))
  expect_true(hexify:::validate_aperture(4))
  expect_true(hexify:::validate_aperture(7))
})

test_that("validate_aperture errors on invalid aperture", {
  expect_error(
    hexify:::validate_aperture(5),
    "must be one of"
  )
  expect_error(
    hexify:::validate_aperture(1),
    "must be one of"
  )
})

test_that("validate_aperture errors on non-numeric input", {
  expect_error(
    hexify:::validate_aperture("3"),
    "must be a single numeric"
  )
})

test_that("validate_aperture errors on vector input", {
  expect_error(
    hexify:::validate_aperture(c(3, 4)),
    "must be a single numeric"
  )
})

# =============================================================================
# MAX_CELL_ID
# =============================================================================

test_that("max_cell_id returns 20 for resolution 0", {
  expect_equal(hexify:::max_cell_id(0, 3), 20)
  expect_equal(hexify:::max_cell_id(0, 4), 20)
  expect_equal(hexify:::max_cell_id(0, 7), 20)
})

test_that("max_cell_id increases with resolution", {
  max_res1 <- hexify:::max_cell_id(1, 3)
  max_res2 <- hexify:::max_cell_id(2, 3)
  max_res3 <- hexify:::max_cell_id(3, 3)

  expect_true(max_res2 > max_res1)
  expect_true(max_res3 > max_res2)
})

test_that("max_cell_id formula is correct", {
  # Formula: 10 * aperture^res + 2
  expect_equal(hexify:::max_cell_id(1, 3), 10 * 3^1 + 2)
  expect_equal(hexify:::max_cell_id(2, 4), 10 * 4^2 + 2)
  expect_equal(hexify:::max_cell_id(1, 7), 10 * 7^1 + 2)
})

# =============================================================================
# VALIDATE_CELL_ID
# =============================================================================

test_that("validate_cell_id returns TRUE for valid cell IDs", {
  result <- hexify:::validate_cell_id(
    c(1, 10, 20), resolution = 0, aperture = 3, warn = FALSE
  )
  expect_true(all(result))
})

test_that("validate_cell_id returns FALSE for invalid cell IDs", {
  result <- hexify:::validate_cell_id(
    c(1, 0, 100), resolution = 0, aperture = 3, warn = FALSE
  )
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("validate_cell_id handles NA values", {
  result <- hexify:::validate_cell_id(
    c(1, NA, 10), resolution = 0, aperture = 3, warn = FALSE
  )
  expect_equal(result, c(TRUE, TRUE, TRUE))
})

test_that("validate_cell_id warns on invalid values", {
  expect_warning(
    hexify:::validate_cell_id(c(0), resolution = 0, aperture = 3),
    "outside valid range"
  )
})

test_that("validate_cell_id errors on non-numeric input", {
  expect_error(
    hexify:::validate_cell_id("1", resolution = 0, aperture = 3),
    "must be numeric"
  )
})

# =============================================================================
# GET_GRID_RESOLUTION
# =============================================================================

test_that("get_grid_resolution extracts resolution from hexify_grid", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  res <- hexify:::get_grid_resolution(grid)
  expect_equal(res, grid$resolution)
})

test_that("get_grid_resolution handles res field name", {
  # Simulate dggridR-style grid with res field
  grid <- list(res = 10)
  res <- hexify:::get_grid_resolution(grid)
  expect_equal(res, 10)
})

test_that("get_grid_resolution returns NULL when field missing", {
  grid <- list(aperture = 3)
  res <- hexify:::get_grid_resolution(grid, require = FALSE)
  expect_null(res)
})

test_that("get_grid_resolution errors when require=TRUE and missing", {
  grid <- list(aperture = 3)
  expect_error(
    hexify:::get_grid_resolution(grid, require = TRUE),
    "missing resolution"
  )
})

# =============================================================================
# CONSTANTS VALUES
# =============================================================================

test_that("EARTH_SURFACE_KM2 is approximately correct", {
  # Earth surface area is approximately 510 million km²
  expect_true(hexify:::EARTH_SURFACE_KM2 > 500000000)
  expect_true(hexify:::EARTH_SURFACE_KM2 < 520000000)
})

test_that("EARTH_RADIUS_KM is approximately correct", {
  # Mean Earth radius is approximately 6371 km
  expect_equal(hexify:::EARTH_RADIUS_KM, 6371.0088, tolerance = 0.01)
})

test_that("VALID_APERTURES contains expected values", {
  expect_equal(hexify:::VALID_APERTURES, c(3L, 4L, 7L))
})

test_that("MIN_RESOLUTION and MAX_RESOLUTION are reasonable", {
  expect_equal(hexify:::MIN_RESOLUTION, 0L)
  expect_equal(hexify:::MAX_RESOLUTION, 30L)
})
