
# tests/testthat/test-assign.R
# Tests for hexify_assign function

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("hexify_assign returns correct structure", {
  setup_icosa()

  lon <- c(0, 45, -120)
  lat <- c(0, 30, -45)

  result <- hexify_assign(lon, lat, effective_res = 3)

  expect_s3_class(result, "data.frame")
  expect_named(result, c(
    "id", "face", "effective_res", "center_lon", "center_lat"
  ))
  expect_equal(nrow(result), 3)
})

test_that("hexify_assign produces valid cell IDs", {
  setup_icosa()

  lon <- c(16.37, 2.35, -3.70)
  lat <- c(48.21, 48.86, 40.42)

  result <- hexify_assign(lon, lat, effective_res = 5)

  # Cell IDs should be formatted strings
  expect_true(all(grepl("^F\\d+:Z3:", result$id)))

  # Face should be valid (0-19)
  expect_true(all(result$face >= 0))
  expect_true(all(result$face <= 19))
})

test_that("hexify_assign centers are valid coordinates", {
  setup_icosa()

  lon <- c(0, 90, -90, 180, -180)
  lat <- c(0, 45, -45, 60, -60)

  result <- hexify_assign(lon, lat, effective_res = 4)

  # Centers should be valid lon/lat
  expect_true(all(result$center_lon >= -180 & result$center_lon <= 180))
  expect_true(all(result$center_lat >= -90 & result$center_lat <= 90))
})

# =============================================================================
# RESOLUTION HANDLING
# =============================================================================

test_that("hexify_assign handles different resolutions", {
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  for (res in c(1, 3, 5, 7)) {
    result <- hexify_assign(lon, lat, effective_res = res)

    expect_equal(result$effective_res, res)
    expect_equal(nrow(result), 1)
  }
})

test_that("hexify_assign resolution affects cell ID length", {
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  result_low <- hexify_assign(lon, lat, effective_res = 2)
  result_high <- hexify_assign(lon, lat, effective_res = 6)

  # Higher resolution should produce longer digit sequences
  # Extract digit part of ID
  digits_low <- sub(".*:Z3:", "", result_low$id)
  digits_high <- sub(".*:Z3:", "", result_high$id)

  expect_lt(nchar(digits_low), nchar(digits_high))
})

# =============================================================================
# PARITY MATCHING
# =============================================================================

test_that("hexify_assign respects match_dggrid_parity", {
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  result_true <- hexify_assign(
    lon, lat, effective_res = 4, match_dggrid_parity = TRUE
  )
  result_false <- hexify_assign(
    lon, lat, effective_res = 4, match_dggrid_parity = FALSE
  )

  # Both should return valid results
  expect_equal(nrow(result_true), 1)
  expect_equal(nrow(result_false), 1)

  # Results may differ based on parity setting
  expect_s3_class(result_true, "data.frame")
  expect_s3_class(result_false, "data.frame")
})

# =============================================================================
# POLYGON GENERATION
# =============================================================================

test_that("hexify_assign with make_polygons=TRUE returns sf", {
  skip_if_not_installed("sf")
  setup_icosa()


  lon <- c(10, 20)
  lat <- c(50, 55)

  result <- hexify_assign(lon, lat, effective_res = 3, make_polygons = TRUE)

  expect_s3_class(result, "sf")
  expect_true("geometry" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify_assign polygons have 6 corners", {
  skip_if_not_installed("sf")
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  result <- hexify_assign(lon, lat, effective_res = 3, make_polygons = TRUE)

  # Get coordinates of polygon
  coords <- sf::st_coordinates(result$geometry[[1]])

  # Polygon should have 7 points (6 corners + closing point)
  expect_equal(nrow(coords), 7)
})

# =============================================================================
# INPUT VALIDATION
# =============================================================================

test_that("hexify_assign validates input lengths", {
  setup_icosa()

  expect_error(
    hexify_assign(c(0, 1), c(0), effective_res = 3),
    "length"
  )
})

test_that("hexify_assign validates resolution", {
  setup_icosa()

  expect_error(
    hexify_assign(0, 0, effective_res = 0),
    "effective_res"
  )

  expect_error(
    hexify_assign(0, 0, effective_res = c(1, 2)),
    "effective_res"
  )
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("hexify_assign handles poles", {
  setup_icosa()

  # North pole
  result_north <- hexify_assign(0, 90, effective_res = 3)
  expect_equal(nrow(result_north), 1)
  expect_true(is.finite(result_north$center_lat))

  # South pole
  result_south <- hexify_assign(0, -90, effective_res = 3)
  expect_equal(nrow(result_south), 1)
  expect_true(is.finite(result_south$center_lat))
})

test_that("hexify_assign handles date line", {
  setup_icosa()

  # Points near date line
  lon <- c(179, -179, 180, -180)
  lat <- c(0, 0, 45, -45)

  result <- hexify_assign(lon, lat, effective_res = 3)
  expect_equal(nrow(result), 4)
})
