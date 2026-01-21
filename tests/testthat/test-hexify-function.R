
# tests/testthat/test-hexify-function.R
# Tests for the main hexify() convenience function
#
# The hexify() function is the primary user-facing API that converts
# a data frame with lon/lat to hexagonal grid cells.
#
# Updated for S4 HexData return type.

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("hexify works with data.frame and area_km2 parameter", {
  df <- data.frame(
    site = c("Vienna", "Paris", "Madrid"),
    lon = c(16.37, 2.35, -3.70),
    lat = c(48.21, 48.86, 40.42)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Returns HexData object
  expect_s4_class(result, "HexData")

  # Original columns preserved in @data
  expect_true("site" %in% names(result@data))
  expect_true("lon" %in% names(result@data))
  expect_true("lat" %in% names(result@data))

  # Cell columns accessible via $ (virtual columns)
  expect_true("cell_id" %in% names(result))
  expect_true("cell_cen_lon" %in% names(result))
  expect_true("cell_cen_lat" %in% names(result))

  # Types
  expect_type(result$cell_id, "double")
  expect_type(result$cell_cen_lon, "double")
  expect_type(result$cell_cen_lat, "double")

  # Valid coordinates
  expect_true(all(result$cell_cen_lon >= -180 & result$cell_cen_lon <= 180))
  expect_true(all(result$cell_cen_lat >= -90 & result$cell_cen_lat <= 90))

  # Row count preserved
  expect_equal(nrow(result), 3)
})

test_that("hexify works with diagonal parameter", {
  df <- data.frame(lon = c(0, 10), lat = c(0, 45))

  result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)

  expect_s4_class(result, "HexData")
  expect_true("cell_id" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify works with grid parameter", {
  df <- data.frame(lon = c(0, 10), lat = c(0, 45))

  grid <- hex_grid(area_km2 = 1000)
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  expect_s4_class(result, "HexData")
  expect_identical(grid_info(result), grid)
})

# =============================================================================
# PARAMETER VALIDATION
# =============================================================================

test_that("hexify requires grid, area_km2, diagonal, or resolution", {
  df <- data.frame(lon = 0, lat = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat"),
               "Either 'grid', 'area_km2', 'diagonal', or 'resolution' must be provided")

  expect_error(hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, diagonal = 50),
               "Provide either 'area_km2' or 'diagonal', not both")
})

test_that("hexify validates column names", {
  df <- data.frame(x = 0, y = 0)

  expect_error(hexify(df, lon = "lon", lat = "lat", area_km2 = 1000),
               "Column 'lon' not found")

  expect_error(hexify(df, lon = "x", lat = "lat", area_km2 = 1000),
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

  result <- hexify(pts, area_km2 = 1000)

  expect_s4_class(result, "HexData")
  # Data slot should be sf
  expect_s3_class(result@data, "sf")
  expect_true("cell_id" %in% names(result))
  expect_true("cell_cen_lon" %in% names(result))
  expect_true("cell_cen_lat" %in% names(result))
  expect_true("site" %in% names(result))
})

# =============================================================================
# APERTURE SUPPORT
# =============================================================================

test_that("hexify handles aperture 3 (ISEA3H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 3)

  expect_s4_class(result, "HexData")
  expect_true("cell_id" %in% names(result))
  expect_type(result$cell_id, "double")
  expect_true(result$cell_id > 0)
  expect_equal(grid_info(result)@aperture, "3")
})

test_that("hexify handles aperture 4 (ISEA4H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 4)

  expect_s4_class(result, "HexData")
  expect_true("cell_id" %in% names(result))
  expect_type(result$cell_id, "double")
  expect_true(result$cell_id > 0)
  expect_equal(grid_info(result)@aperture, "4")
})

test_that("hexify handles aperture 7 (ISEA7H)", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 7)

  expect_s4_class(result, "HexData")
  expect_true("cell_id" %in% names(result))
  expect_type(result$cell_id, "double")
  expect_true(result$cell_id > 0)
  expect_equal(grid_info(result)@aperture, "7")
})

test_that("hexify handles mixed aperture 4/3 (ISEA43H)", {
  df <- data.frame(lon = 0, lat = 45)

  # Test with string "4/3"
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = "4/3")

  expect_s4_class(result, "HexData")
  expect_true("cell_id" %in% names(result))
  expect_type(result$cell_id, "double")
  expect_true(result$cell_id > 0)
  expect_equal(grid_info(result)@aperture, "4/3")
})

test_that("hexify rejects unsupported apertures with clear error", {
  df <- data.frame(lon = 0, lat = 45)

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 5),
    "Aperture must be 3, 4, 7"
  )

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 2),
    "Aperture must be 3, 4, 7"
  )

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = "invalid"),
    "Aperture must be 3, 4, 7"
  )
})

# =============================================================================
# OUTPUT COLUMNS
# =============================================================================

test_that("hexify returns cell_area_km2 and cell_diag_km columns", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 3)

  expect_true("cell_area_km2" %in% names(result))
  expect_true("cell_diag_km" %in% names(result))
  expect_type(result$cell_area_km2, "double")
  expect_type(result$cell_diag_km, "double")

  # Values should be reasonable
  expect_true(result$cell_area_km2 > 100 && result$cell_area_km2 < 10000)
  expect_true(result$cell_diag_km > 10 && result$cell_diag_km < 200)
})

test_that("hexify cell_area_km2 and cell_diag_km are consistent across rows", {
  df <- data.frame(lon = c(0, 10, -5), lat = c(45, 30, -20))

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = 3)

  # All rows should have same area/diag
  expect_equal(length(unique(result$cell_area_km2)), 1)
  expect_equal(length(unique(result$cell_diag_km)), 1)
})

test_that("hexify cell_area_km2 and cell_diag_km work for all apertures", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 30))

  for (ap in c(3, 4, 7)) {
    result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = ap)

    expect_true("cell_area_km2" %in% names(result), info = sprintf("aperture %d", ap))
    expect_true("cell_diag_km" %in% names(result), info = sprintf("aperture %d", ap))
    expect_true(result$cell_diag_km[1] > 0)
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
    result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000),
    "coordinate pairs contain NA"
  )

  expect_equal(nrow(result), 3)
})

test_that("hexify errors on all NA coordinates", {
  df <- data.frame(lon = as.numeric(c(NA, NA)), lat = as.numeric(c(NA, NA)))

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area_km2 = 1000),
    "All coordinates are NA"
  )
})

test_that("hexify validates numeric coordinates", {
  df <- data.frame(lon = c("a", "b"), lat = c("c", "d"))

  expect_error(
    hexify(df, lon = "lon", lat = "lat", area_km2 = 1000),
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

  result <- hexify(pts_3857, area_km2 = 1000)

  expect_s4_class(result, "HexData")
  expect_s3_class(result@data, "sf")
  expect_true("cell_id" %in% names(result))
})

test_that("hexify with sf handles NA CRS", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(16.37, 2.35),
    lat = c(48.21, 48.86)
  )

  pts <- sf::st_as_sf(df, coords = c("lon", "lat"))
  sf::st_crs(pts) <- NA

  result <- hexify(pts, area_km2 = 1000)
  expect_s4_class(result, "HexData")
})

# =============================================================================
# MIXED APERTURE DETAILED
# =============================================================================

test_that("hexify works with mixed aperture 4/3", {
  df <- data.frame(lon = 0, lat = 45)

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000, aperture = "4/3")
  expect_true("cell_id" %in% names(result))
  expect_equal(grid_info(result)@aperture, "4/3")
})

# =============================================================================
# RESOLUTION ROUNDING
# =============================================================================

test_that("hexify respects resround parameter", {
  df <- data.frame(lon = 0, lat = 45)

  result_up <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000,
                      aperture = 3, resround = "up")
  result_down <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000,
                        aperture = 3, resround = "down")
  result_nearest <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000,
                           aperture = 3, resround = "nearest")

  expect_true("cell_area_km2" %in% names(result_up))
  expect_true("cell_area_km2" %in% names(result_down))
  expect_true("cell_area_km2" %in% names(result_nearest))
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("hexify handles single point", {
  df <- data.frame(lon = 0, lat = 0)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_equal(nrow(result), 1)
  expect_true(is.finite(result$cell_id))
})

test_that("hexify handles poles", {
  df <- data.frame(lon = c(0, 0), lat = c(90, -90))

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$cell_id)))
})

test_that("hexify handles date line", {
  df <- data.frame(
    lon = c(179, -179, 180, -180),
    lat = c(0, 0, 45, -45)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_equal(nrow(result), 4)
  expect_true(all(is.finite(result$cell_id)))
})

test_that("hexify handles many points", {
  set.seed(42)
  df <- data.frame(
    lon = runif(100, -180, 180),
    lat = runif(100, -90, 90)
  )

  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_equal(nrow(result), 100)
  expect_true(all(is.finite(result$cell_id)))
})

# =============================================================================
# HEXDATA ACCESSORS
# =============================================================================

test_that("HexData accessors work correctly", {
  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # grid() accessor
  g <- grid_info(result)
  expect_s4_class(g, "HexGridInfo")

  # cells() accessor
  c <- cells(result)
  expect_type(c, "double")
  expect_equal(length(c), n_cells(result))

  # n_cells() accessor
  expect_type(n_cells(result), "integer")

  # nrow/ncol
  expect_equal(nrow(result), 3)
  expect_true(ncol(result) >= 6)

  # names
  expect_true("cell_id" %in% names(result))
})

test_that("HexData subsetting works", {
  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55), val = c(1, 2, 3))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Subset rows
  subset <- result[1:2, ]
  expect_s4_class(subset, "HexData")
  expect_equal(nrow(subset), 2)

  # Column access
  expect_equal(result$val, c(1, 2, 3))
  expect_equal(result[["val"]], c(1, 2, 3))
})

test_that("as.data.frame works on HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  df_out <- as.data.frame(result)
  expect_s3_class(df_out, "data.frame")
  expect_true("cell_id" %in% names(df_out))
})

# =============================================================================
# HEX_GRID CONSTRUCTOR
# =============================================================================

test_that("hex_grid constructor works with area_km2", {
  grid <- hex_grid(area_km2 = 1000)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "3")
  expect_true(grid@resolution >= 0)
  expect_true(grid@resolution <= 30)
})

test_that("hex_grid constructor works with resolution", {
  grid <- hex_grid(resolution = 8)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@resolution, 8L)
})

test_that("hex_grid rejects both area_km2 and resolution", {
  expect_error(
    hex_grid(area_km2 = 1000, resolution = 8),
    "Provide either 'area_km2' or 'resolution', not both"
  )
})

test_that("hex_grid supports different apertures", {
  grid3 <- hex_grid(area_km2 = 1000, aperture = 3)
  grid4 <- hex_grid(area_km2 = 1000, aperture = 4)
  grid7 <- hex_grid(area_km2 = 1000, aperture = 7)

  expect_equal(grid3@aperture, "3")
  expect_equal(grid4@aperture, "4")
  expect_equal(grid7@aperture, "7")
})

test_that("hex_grid supports mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")

  expect_equal(grid@aperture, "4/3")
})
