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

test_that("hexify() works with diagonal parameter", {
  df <- data.frame(lon = c(0, 10), lat = c(0, 45))

  result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)

  expect_true("hex_id" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify() requires either area or diagonal", {
  df <- data.frame(lon = 0, lat = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat"),
               "Either 'area'.*or 'diagonal'.*must be provided")

  expect_error(hexify(df, lon = "lon", lat = "lat", area = 1000, diagonal = 50),
               "Provide either 'area' or 'diagonal', not both")
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

test_that("hexify() handles aperture 4 (ISEA4H)", {
  df <- data.frame(lon = 0, lat = 45)

  result_ap4 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)

  # Should work and return integer SEQNUM
  expect_true("hex_id" %in% names(result_ap4))
  expect_type(result_ap4$hex_id, "integer")
  expect_true(result_ap4$hex_id > 0)
})

test_that("hexify() handles aperture 7 (ISEA7H)", {
  df <- data.frame(lon = 0, lat = 45)

  result_ap7 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)

  # Should work and return integer SEQNUM
  expect_true("hex_id" %in% names(result_ap7))
  expect_type(result_ap7$hex_id, "integer")
  expect_true(result_ap7$hex_id > 0)
})

test_that("hexify() rejects unsupported apertures with clear error", {
  df <- data.frame(lon = 0, lat = 45)

  # Apertures other than 3, 4, 7 and "4/3" are not supported
  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 5),
    "Aperture must be 3, 4, or 7"
  )
  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 2),
    "Aperture must be 3, 4, or 7"
  )
  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "invalid"),
    "Aperture must be 3, 4, 7, or '4/3'"
  )
})

test_that("hexify() handles mixed aperture 4/3 (ISEA43H)", {
  df <- data.frame(lon = 0, lat = 45)

  # Test with string "4/3"
  result_43 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "4/3")

  # Should work and return integer SEQNUM
  expect_true("hex_id" %in% names(result_43))
  expect_type(result_43$hex_id, "integer")
  expect_true(result_43$hex_id > 0)

  # Test with explicit mixed_aperture_level
  result_43_level <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                            aperture = "4/3", mixed_aperture_level = 4)
  expect_true("hex_id" %in% names(result_43_level))
  expect_type(result_43_level$hex_id, "integer")
})

test_that("hexify() returns hex_area and hex_diag columns", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)

  # Check that hex_area and hex_diag columns are present

  expect_true("hex_area" %in% names(result))
  expect_true("hex_diag" %in% names(result))

  # Check types
  expect_type(result$hex_area, "double")
  expect_type(result$hex_diag, "double")

  # Check values are reasonable (area should be within an order of magnitude of target)
  expect_true(result$hex_area > 100 && result$hex_area < 10000)
  expect_true(result$hex_diag > 10 && result$hex_diag < 200)
})

test_that("hexify() hex_area and hex_diag are consistent across apertures", {
  df <- data.frame(lon = c(0, 10, -5), lat = c(45, 30, -20))

  result_ap3 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)
  result_ap4 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)
  result_ap7 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)

  # All rows within same result should have same area/diag
  expect_equal(length(unique(result_ap3$hex_area)), 1)
  expect_equal(length(unique(result_ap3$hex_diag)), 1)
  expect_equal(length(unique(result_ap4$hex_area)), 1)
  expect_equal(length(unique(result_ap7$hex_area)), 1)

  # Area and diagonal should be positively related: diag ~ sqrt(area)
  expect_true(result_ap3$hex_diag[1] > 0)
  expect_true(result_ap4$hex_diag[1] > 0)
  expect_true(result_ap7$hex_diag[1] > 0)
})
