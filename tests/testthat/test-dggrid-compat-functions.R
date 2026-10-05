
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
    projection = "GNOMONIC"
  )
  expect_warning(from_dggrid(dggs), "ISEA and FULLER")
})

test_that("from_dggrid and as_dggrid carry the FULLER projection", {
  dggs <- list(res = 5L, aperture = 4L, topology = "HEXAGON", projection = "FULLER")
  grid <- from_dggrid(dggs)
  expect_equal(grid$projection, "FULLER")
  g <- extract_grid(grid)
  expect_equal(g@projection, "fuller")
  expect_equal(as_dggrid(g)$projection, "FULLER")
  expect_true(dggrid_is_compatible(dggs))
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

test_that("from_dggrid carries the orientation", {
  dggs <- list(
    res = 5L,
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA",
    pole_lon_deg = -20,
    pole_lat_deg = 10,
    azimuth_deg = 45
  )
  grid <- extract_grid(from_dggrid(dggs))
  expect_equal(unname(grid@orientation), c(-20, 10, 45))

  # DGGRID writes the standard latitude to eight decimals
  dggs$pole_lon_deg <- 11.25
  dggs$pole_lat_deg <- 58.28252559
  dggs$azimuth_deg <- 0
  expect_identical(extract_grid(from_dggrid(dggs))@orientation, ISEA_ORIENTATION)
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
  dggs <- list(aperture = 3L, topology = "HEXAGON", projection = "GNOMONIC")
  expect_error(dggrid_is_compatible(dggs), "ISEA or FULLER")

  # Wrong topology
  dggs <- list(aperture = 3L, topology = "DIAMOND", projection = "ISEA")
  expect_error(dggrid_is_compatible(dggs), "HEXAGON topology")

  # Wrong aperture
  dggs <- list(aperture = 5L, topology = "HEXAGON", projection = "ISEA")
  expect_error(dggrid_is_compatible(dggs), "Aperture must be")
})

test_that("dggrid_is_compatible returns FALSE for incompatible grids (strict=FALSE)", {
  dggs <- list(aperture = 3L, topology = "HEXAGON", projection = "GNOMONIC")
  expect_false(dggrid_is_compatible(dggs, strict = FALSE))

  dggs <- list(aperture = 3L, topology = "DIAMOND", projection = "ISEA")
  expect_false(dggrid_is_compatible(dggs, strict = FALSE))
})

test_that("dggrid_is_compatible accepts any orientation on the sphere", {
  dggs <- list(
    aperture = 3L,
    topology = "HEXAGON",
    projection = "ISEA",
    pole_lon_deg = 0,
    pole_lat_deg = 90,
    azimuth_deg = 30
  )
  expect_true(dggrid_is_compatible(dggs))

  dggs$pole_lat_deg <- 95
  expect_error(dggrid_is_compatible(dggs), "pole_lat_deg")
})

test_that("round-trip: hex_grid -> dggridR -> hex_grid keeps the orientation", {
  grid <- hex_grid(resolution = 6, orientation = c(30, -45, 100))
  dggs <- as_dggrid(grid)
  expect_equal(c(dggs$pole_lon_deg, dggs$pole_lat_deg, dggs$azimuth_deg),
               c(30, -45, 100))
  expect_equal(extract_grid(from_dggrid(dggs))@orientation, grid@orientation)
})

test_that("round-trip: hexify_grid -> dggridR -> hexify_grid", {
  original <- hexify_grid(area = 1000, aperture = 3)
  dggs <- as_dggrid(original)
  recovered <- from_dggrid(dggs)

  expect_equal(original$resolution, recovered$resolution)
  expect_equal(original$aperture, recovered$aperture)
})

test_that("as_dggrid warns that a dggs cannot carry a body radius", {
  grid <- hex_grid(area_km2 = 100000, radius_km = "mars")

  expect_warning(dggs <- as_dggrid(grid), "carries no body radius")
  expect_warning(as_dggrid(grid), "3389.5")

  # The rest of the conversion still happens
  expect_s3_class(dggs, "dggs")
  expect_equal(dggs$res, grid@resolution)
})

test_that("as_dggrid converts an Earth grid without warning", {
  expect_no_warning(as_dggrid(hex_grid(area_km2 = 100000)))
  expect_no_warning(as_dggrid(hexify_grid(area = 1000, aperture = 3)))
})

test_that("from_dggrid reads Earth unless told otherwise", {
  dggs <- suppressWarnings(as_dggrid(hex_grid(area_km2 = 100000, radius_km = "mars")))

  earth <- from_dggrid(dggs)
  expect_equal(hexify:::grid_radius_km(earth), hexify:::EARTH_RADIUS_KM)

  mars <- from_dggrid(dggs, radius_km = "mars")
  expect_equal(hexify:::grid_radius_km(mars), 3389.5)

  # Cell area follows the radius given, so the grid round-trips
  original <- hex_grid(area_km2 = 100000, radius_km = "mars")
  expect_equal(mars$area, original@area_km2, tolerance = 1e-8)
  expect_gt(earth$area, mars$area)
})
