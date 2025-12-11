# tests/testthat/test-hexify.R
# Tests for the main hexify() convenience function

test_that("hexify() works with data.frame and area parameter", {
  df <- data.frame(
    site = c("Vienna", "Paris", "Madrid"),
    lon = c(16.37, 2.35, -3.70),
    lat = c(48.21, 48.86, 40.42)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

  # Check that original columns are preserved
  expect_true("site" %in% names(result))
  expect_true("lon" %in% names(result))
  expect_true("lat" %in% names(result))

  # Check that hex columns are added
  expect_true("hex_id" %in% names(result))
  expect_true("hex_cen_lon" %in% names(result))
  expect_true("hex_cen_lat" %in% names(result))

  # Check types (hex_id is now integer SEQNUM)
  expect_type(result$hex_id, "integer")
  expect_type(result$hex_cen_lon, "double")
  expect_type(result$hex_cen_lat, "double")

  # Check that centers are valid coordinates
  expect_true(all(result$hex_cen_lon >= -180 & result$hex_cen_lon <= 180))
  expect_true(all(result$hex_cen_lat >= -90 & result$hex_cen_lat <= 90))

  # Check row count preserved
  expect_equal(nrow(result), 3)
})

test_that("hexify() works with spacing parameter", {
  df <- data.frame(lon = c(0, 10), lat = c(0, 45))

  result <- hexify(df, lon = "lon", lat = "lat", spacing = 50)

  expect_true("hex_id" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify() requires either area or spacing", {
  df <- data.frame(lon = 0, lat = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat"),
               "Either 'area'.*or 'spacing'.*must be provided")

  expect_error(hexify(df, lon = "lon", lat = "lat", area = 1000, spacing = 50),
               "Provide either 'area' or 'spacing', not both")
})

test_that("hexify() validates column names", {
  df <- data.frame(x = 0, y = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat", area = 1000),
               "Column 'lon' not found")

  expect_error(hexify(df, lon = "x", lat = "lat", area = 1000),
               "Column 'lat' not found")
})

test_that("hexify() works with sf objects", {
  skip_if_not_installed("sf")

  df <- data.frame(
    site = c("Vienna", "Paris"),
    lon = c(16.37, 2.35),
    lat = c(48.21, 48.86)
  )
  pts <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

  result <- hexify(pts, area = 1000)

  # Check sf class preserved
  expect_s3_class(result, "sf")

  # Check hex columns added
  expect_true("hex_id" %in% names(result))
  expect_true("hex_cen_lon" %in% names(result))
  expect_true("hex_cen_lat" %in% names(result))

  # Check original columns preserved
  expect_true("site" %in% names(result))
})

test_that("hexify() handles aperture 3 (ISEA3H)", {
  df <- data.frame(lon = 0, lat = 45)

  result_ap3 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)

  # Should work and return integer SEQNUM
  expect_true("hex_id" %in% names(result_ap3))
  expect_type(result_ap3$hex_id, "integer")
  expect_true(result_ap3$hex_id > 0)
})

test_that("hexify() rejects unsupported apertures with clear error", {
  df <- data.frame(lon = 0, lat = 45)

  # Aperture 4 and 7 not yet supported for DGGRID-compatible SEQNUM
  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4),
    "only supports aperture 3"
  )
  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7),
    "only supports aperture 3"
  )
})
