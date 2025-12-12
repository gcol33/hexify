# tests/testthat/test-hexify-function.R
# Tests for the main hexify() convenience function
#
# The hexify() function is the primary user-facing API that converts
# a data frame with lon/lat to hexagonal grid cells.

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("hexify works with data.frame and area parameter", {
  df <- data.frame(
    site = c("Vienna", "Paris", "Madrid"),
    lon = c(16.37, 2.35, -3.70),
    lat = c(48.21, 48.86, 40.42)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

  # Original columns preserved
  expect_true("site" %in% names(result))
  expect_true("lon" %in% names(result))
  expect_true("lat" %in% names(result))

  # Hex columns added
  expect_true("hex_id" %in% names(result))
  expect_true("hex_cen_lon" %in% names(result))
  expect_true("hex_cen_lat" %in% names(result))

  # Types
  expect_type(result$hex_id, "double")
  expect_type(result$hex_cen_lon, "double")
  expect_type(result$hex_cen_lat, "double")

  # Valid coordinates
  expect_true(all(result$hex_cen_lon >= -180 & result$hex_cen_lon <= 180))
  expect_true(all(result$hex_cen_lat >= -90 & result$hex_cen_lat <= 90))

  # Row count preserved
  expect_equal(nrow(result), 3)
})

test_that("hexify works with diagonal parameter", {
  df <- data.frame(lon = c(0, 10), lat = c(0, 45))

  result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)

  expect_true("hex_id" %in% names(result))
  expect_equal(nrow(result), 2)
})

# =============================================================================
# PARAMETER VALIDATION
# =============================================================================

test_that("hexify requires either area or diagonal", {
  df <- data.frame(lon = 0, lat = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat"),
               "Either 'area'.*or 'diagonal'.*must be provided")

  expect_error(hexify(df, lon = "lon", lat = "lat", area = 1000, diagonal = 50),
               "Provide either 'area' or 'diagonal', not both")
})

test_that("hexify validates column names", {
  df <- data.frame(x = 0, y = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat", area = 1000),
               "Column 'lon' not found")

  expect_error(hexify(df, lon = "x", lat = "lat", area = 1000),
               "Column 'lat' not found")
})

# =============================================================================
# SF INTEGRATION
# =============================================================================

test_that("hexify works with sf objects", {
  skip_if_not_installed("sf")

  df <- data.frame(
    site = c("Vienna", "Paris"),
    lon = c(16.37, 2.35),
    lat = c(48.21, 48.86)
  )
  pts <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

  result <- hexify(pts, area = 1000)

  expect_s3_class(result, "sf")
  expect_true("hex_id" %in% names(result))
  expect_true("hex_cen_lon" %in% names(result))
  expect_true("hex_cen_lat" %in% names(result))
  expect_true("site" %in% names(result))
})

# =============================================================================
# APERTURE SUPPORT
# =============================================================================

test_that("hexify handles aperture 3 (ISEA3H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)

  expect_true("hex_id" %in% names(result))
  expect_type(result$hex_id, "double")
  expect_true(result$hex_id > 0)
})

test_that("hexify handles aperture 4 (ISEA4H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)

  expect_true("hex_id" %in% names(result))
  expect_type(result$hex_id, "double")
  expect_true(result$hex_id > 0)
})

test_that("hexify handles aperture 7 (ISEA7H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)

  expect_true("hex_id" %in% names(result))
  expect_type(result$hex_id, "double")
  expect_true(result$hex_id > 0)
})

test_that("hexify handles mixed aperture 4/3 (ISEA43H)", {
  df <- data.frame(lon = 0, lat = 45)

  # Test with string "4/3"
  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "4/3")

  expect_true("hex_id" %in% names(result))
  expect_type(result$hex_id, "double")
  expect_true(result$hex_id > 0)

  # Test with explicit mixed_aperture_level
  result2 <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                    aperture = "4/3", mixed_aperture_level = 4)
  expect_true("hex_id" %in% names(result2))
})

test_that("hexify rejects unsupported apertures with clear error", {
  df <- data.frame(lon = 0, lat = 45)

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

# =============================================================================
# OUTPUT COLUMNS
# =============================================================================

test_that("hexify returns hex_area and hex_diag columns", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)

  expect_true("hex_area" %in% names(result))
  expect_true("hex_diag" %in% names(result))
  expect_type(result$hex_area, "double")
  expect_type(result$hex_diag, "double")

  # Values should be reasonable
  expect_true(result$hex_area > 100 && result$hex_area < 10000)
  expect_true(result$hex_diag > 10 && result$hex_diag < 200)
})

test_that("hexify hex_area and hex_diag are consistent across rows", {
  df <- data.frame(lon = c(0, 10, -5), lat = c(45, 30, -20))

  result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 3)

  # All rows should have same area/diag
  expect_equal(length(unique(result$hex_area)), 1)
  expect_equal(length(unique(result$hex_diag)), 1)
})

test_that("hexify hex_area and hex_diag work for all apertures", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 30))

  for (ap in c(3, 4, 7)) {
    result <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = ap)

    expect_true("hex_area" %in% names(result), info = sprintf("aperture %d", ap))
    expect_true("hex_diag" %in% names(result), info = sprintf("aperture %d", ap))
    expect_true(result$hex_diag[1] > 0)
  }
})

# =============================================================================
# NA HANDLING
# =============================================================================

test_that("hexify handles NA coordinate values", {
  df <- data.frame(
    lon = c(0, NA, 10),
    lat = c(45, 46, NA)
  )

  expect_warning(
    result <- hexify(df, lon = "lon", lat = "lat", area = 1000),
    "coordinate pairs contain NA"
  )

  expect_equal(nrow(result), 3)
})

test_that("hexify errors on all NA coordinates", {
  df <- data.frame(lon = as.numeric(c(NA, NA)), lat = as.numeric(c(NA, NA)))

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000),
    "All coordinates are NA"
  )
})

test_that("hexify validates numeric coordinates", {
  df <- data.frame(lon = c("a", "b"), lat = c("c", "d"))

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000),
    "Coordinates must be numeric"
  )
})

# =============================================================================
# SF INTEGRATION - ADVANCED
# =============================================================================

test_that("hexify with sf handles non-4326 CRS", {
  skip_if_not_installed("sf")

  df <- data.frame(
    site = c("Vienna", "Paris"),
    lon = c(16.37, 2.35),
    lat = c(48.21, 48.86)
  )

  pts_4326 <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  pts_3857 <- sf::st_transform(pts_4326, crs = 3857)

  result <- hexify(pts_3857, area = 1000)

  expect_s3_class(result, "sf")
  expect_true("hex_id" %in% names(result))
})

test_that("hexify with sf handles NA CRS", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(16.37, 2.35),
    lat = c(48.21, 48.86)
  )

  pts <- sf::st_as_sf(df, coords = c("lon", "lat"))
  sf::st_crs(pts) <- NA

  result <- hexify(pts, area = 1000)
  expect_s3_class(result, "sf")
})

# =============================================================================
# MIXED APERTURE DETAILED
# =============================================================================

test_that("hexify mixed_aperture_level validation", {
  df <- data.frame(lon = 0, lat = 45)

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area = 1000,
           aperture = "4/3", mixed_aperture_level = -1),
    "mixed_aperture_level must be between"
  )
})

test_that("hexify mixed aperture produces different results with different levels", {
  df <- data.frame(lon = 0, lat = 45)

  result_level2 <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                          aperture = "4/3", mixed_aperture_level = 2)
  result_level4 <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                          aperture = "4/3", mixed_aperture_level = 4)

  expect_true(result_level2$hex_id != result_level4$hex_id ||
              result_level2$hex_area != result_level4$hex_area)
})

# =============================================================================
# RESOLUTION ROUNDING
# =============================================================================

test_that("hexify respects resround parameter", {
  df <- data.frame(lon = 0, lat = 45)

  result_up <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                      aperture = 3, resround = "up")
  result_down <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                        aperture = 3, resround = "down")
  result_nearest <- hexify(df, lon = "lon", lat = "lat", area = 1000,
                           aperture = 3, resround = "nearest")

  expect_true("hex_area" %in% names(result_up))
  expect_true("hex_area" %in% names(result_down))
  expect_true("hex_area" %in% names(result_nearest))
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("hexify handles single point", {
  df <- data.frame(lon = 0, lat = 0)
  result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

  expect_equal(nrow(result), 1)
  expect_true(is.finite(result$hex_id))
})

test_that("hexify handles poles", {
  df <- data.frame(lon = c(0, 0), lat = c(90, -90))

  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$hex_id)))
})

test_that("hexify handles date line", {
  df <- data.frame(
    lon = c(179, -179, 180, -180),
    lat = c(0, 0, 45, -45)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_equal(nrow(result), 4)
  expect_true(all(is.finite(result$hex_id)))
})

test_that("hexify handles many points", {
  set.seed(42)
  df <- data.frame(
    lon = runif(100, -180, 180),
    lat = runif(100, -90, 90)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area = 10000)

  expect_equal(nrow(result), 100)
  expect_true(all(is.finite(result$hex_id)))
})
