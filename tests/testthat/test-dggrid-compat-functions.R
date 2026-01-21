
# tests/testthat/test-dggrid-compat-functions.R
# Tests for dggridR compatibility layer functions

test_that("as_dggrid converts hexify_grid to dggridR format", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  dggs <- as_dggrid(grid)

  expect_type(dggs, "list")
  expect_equal(dggs$aperture, 3)
  expect_equal(dggs$res, grid$resolution)
  expect_equal(dggs$topology, "HEXAGON")
  expect_equal(dggs$projection, "ISEA")
  expect_equal(dggs$precision, 7L)
  expect_equal(dggs$pole_lon_deg, 11.25)
})

test_that("as_dggrid rejects non-hexify_grid objects", {
  expect_error(as_dggrid(list(a = 1)), "must be a hexify_grid object")
  expect_error(as_dggrid(data.frame()), "must be a hexify_grid object")
})

test_that("from_dggrid converts dggridR format to hexify_grid", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA"
  )

  grid <- from_dggrid(dggs)

  expect_s3_class(grid, "hexify_grid")
  expect_equal(grid$resolution, 5L)
  expect_equal(grid$aperture, 3L)
})

test_that("from_dggrid rejects invalid input", {
  expect_error(from_dggrid("not a list"), "must be a list")
  expect_error(from_dggrid(list()), "missing required fields")
  expect_error(from_dggrid(list(res = 5)), "missing required fields")
})

test_that("from_dggrid warns on unsupported projection", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "HEXAGON",
    projection = "FULLER"
  )
  expect_warning(from_dggrid(dggs), "Only ISEA projection")
})

test_that("from_dggrid warns on unsupported topology", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "DIAMOND",
    projection = "ISEA"
  )
  expect_warning(from_dggrid(dggs), "Only HEXAGON topology")
})

test_that("from_dggrid rejects unsupported aperture", {
  dggs <- list(
    res = 5L,
    aperture = 5L,
    topology = "HEXAGON",
    projection = "ISEA"
  )
  expect_error(from_dggrid(dggs), "Aperture 5 not supported")
})

test_that("from_dggrid warns on non-default orientation", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA",
    pole_lon_deg = 0
  )
  expect_warning(from_dggrid(dggs), "Non-default pole_lon_deg")

  dggs$pole_lon_deg <- 11.25
  dggs$pole_lat_deg <- 0
  expect_warning(from_dggrid(dggs), "Non-default pole_lat_deg")

  dggs$pole_lat_deg <- 58.28252559
  dggs$azimuth_deg <- 45
  expect_warning(from_dggrid(dggs), "Non-default azimuth_deg")
})

test_that("dggrid_is_compatible validates compatible grids", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA"
  )
  expect_true(dggrid_is_compatible(dggs))
})

test_that("dggrid_is_compatible rejects incompatible grids (strict=TRUE)", {
  # Not a list
  expect_error(dggrid_is_compatible("not a list"), "not compatible")

  # Wrong projection
  dggs <- list(aperture = 3L, topology = "HEXAGON", projection = "FULLER")
  expect_error(dggrid_is_compatible(dggs), "ISEA projection")

  # Wrong topology
  dggs <- list(aperture = 3L, topology = "DIAMOND", projection = "ISEA")
  expect_error(dggrid_is_compatible(dggs), "HEXAGON topology")

  # Wrong aperture
  dggs <- list(aperture = 5L, topology = "HEXAGON", projection = "ISEA")
  expect_error(dggrid_is_compatible(dggs), "Aperture must be")
})

test_that("dggrid_is_compatible returns FALSE for incompatible grids (strict=FALSE)", {
  dggs <- list(aperture = 3L, topology = "HEXAGON", projection = "FULLER")
  expect_false(dggrid_is_compatible(dggs, strict = FALSE))

  dggs <- list(aperture = 3L, topology = "DIAMOND", projection = "ISEA")
  expect_false(dggrid_is_compatible(dggs, strict = FALSE))
})

test_that("dggrid_is_compatible detects non-default orientation", {
  dggs <- list(
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA",
    pole_lon_deg = 0
  )
  expect_error(dggrid_is_compatible(dggs), "pole_lon_deg")
})

test_that("round-trip: hexify_grid -> dggridR -> hexify_grid", {
  original <- hexify_grid(area = 1000, aperture = 3)
  dggs <- as_dggrid(original)
  recovered <- from_dggrid(dggs)

  expect_equal(original$resolution, recovered$resolution)
  expect_equal(original$aperture, recovered$aperture)
})
