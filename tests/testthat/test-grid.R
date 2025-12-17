# tests/testthat/test-grid.R
# Tests for hexify_grid.R - grid construction and validation

# =============================================================================
# GRID CONSTRUCTION
# =============================================================================

test_that("hexify_grid creates valid grid object", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  expect_s3_class(grid, "hexify_grid")
  expect_s3_class(grid, "dggs")
  expect_true("aperture" %in% names(grid))
  expect_true("resolution" %in% names(grid))
  expect_true("area" %in% names(grid))
})

test_that("hexify_grid sets correct aperture", {
  for (ap in c(3, 4, 7)) {
    grid <- hexify_grid(area = 1000, aperture = ap)
    expect_equal(grid$aperture, ap, info = sprintf("aperture %d", ap))
  }
})

test_that("hexify_grid sets correct index_type based on aperture", {
  grid3 <- hexify_grid(area = 1000, aperture = 3)
  expect_equal(grid3$index_type, "z3")

  grid4 <- hexify_grid(area = 1000, aperture = 4)
  expect_equal(grid4$index_type, "zorder")

  grid7 <- hexify_grid(area = 1000, aperture = 7)
  expect_equal(grid7$index_type, "z7")
})

test_that("hexify_grid resolution rounding works", {
  area <- 1000

  grid_nearest <- hexify_grid(area = area, resround = "nearest")
  grid_up <- hexify_grid(area = area, resround = "up")
  grid_down <- hexify_grid(area = area, resround = "down")

  expect_true(grid_up$resolution >= grid_down$resolution)
  expect_true(grid_nearest$resolution >= grid_down$resolution)
  expect_true(grid_nearest$resolution <= grid_up$resolution)
})

test_that("hexify_grid validates topology", {
  expect_error(
    hexify_grid(area = 1000, topology = "TRIANGLE"),
    "HEXAGON"
  )
})

test_that("hexify_grid validates projection", {
  expect_error(
    hexify_grid(area = 1000, projection = "FULLER"),
    "ISEA"
  )
})

test_that("hexify_grid validates aperture", {
  expect_error(
    hexify_grid(area = 1000, aperture = 5),
    "Aperture"
  )
})

test_that("hexify_grid validates resround", {
  expect_error(
    hexify_grid(area = 1000, resround = "invalid"),
    "resround"
  )
})

test_that("hexify_grid resolution stays in valid range", {
  # Very small area -> high resolution (capped at MAX)
  grid_small <- hexify_grid(area = 0.001, aperture = 3)
  expect_lte(grid_small$resolution, MAX_RESOLUTION)

  # Large area -> low resolution (capped at MIN)
  # Use area that still gives valid (non-NaN) resolution
  grid_large <- hexify_grid(area = 1e8, aperture = 3)
  expect_gte(grid_large$resolution, MIN_RESOLUTION)
})

test_that("hexify_grid includes dggridR-compatible fields", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # Check for dggridR compatibility fields
  expect_true("res" %in% names(grid))
  expect_true("topology_family" %in% names(grid))
  expect_true("pole_lon_deg" %in% names(grid))
  expect_true("pole_lat_deg" %in% names(grid))
  expect_true("azimuth_deg" %in% names(grid))
})

# =============================================================================
# GRID VERIFICATION
# =============================================================================

test_that("dgverify accepts valid grid", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  expect_true(dgverify(grid))
})

test_that("dgverify rejects non-grid objects", {
  expect_error(dgverify("not a grid"), "hexify_grid")
  expect_error(dgverify(list(x = 1)), "hexify_grid")
})

test_that("dgverify checks required fields", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # Remove required field
  bad_grid <- grid
  bad_grid$aperture <- NULL
  class(bad_grid) <- c("hexify_grid", "dggs", "list")

  expect_error(dgverify(bad_grid), "missing")
})

test_that("dgverify warns on non-HEXAGON topology", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  grid$topology <- "OTHER"

  expect_warning(dgverify(grid), "HEXAGON")
})

test_that("dgverify warns on non-ISEA projection", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  grid$projection <- "OTHER"

  expect_warning(dgverify(grid), "ISEA")
})

# =============================================================================
# PRINT METHOD (uses default list printing now)
# =============================================================================

test_that("print outputs grid information", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # Default list print shows fields
  expect_output(print(grid), "area")
  expect_output(print(grid), "resolution")
  expect_output(print(grid), "aperture")
})

# =============================================================================
# RESOLUTION CALCULATION
# =============================================================================

test_that("calculate_resolution_for_area is monotonic", {
  areas <- c(10000, 1000, 100, 10, 1)
  resolutions <- sapply(areas, calculate_resolution_for_area, aperture = 3)

  # Smaller area -> higher resolution
  expect_true(all(diff(resolutions) > 0))
})

test_that("calculate_resolution_for_area works for all apertures", {
  area <- 1000

  for (ap in c(3, 4, 7)) {
    res <- calculate_resolution_for_area(area, aperture = ap)
    expect_true(is.finite(res))
    expect_true(res > 0)
  }
})
