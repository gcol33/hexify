# tests/testthat/test-raster.R
# Tests for hex_extract() and hex_zonal()

test_that("hex_extract works with synthetic SpatRaster", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 10, ncols = 10,
                    xmin = -10, xmax = 10, ymin = 40, ymax = 55)
  terra::values(r) <- seq_len(100)
  names(r) <- "elevation"

  df <- data.frame(lon = c(0, 5), lat = c(45, 50))
  g <- hex_grid(area_km2 = 500)
  hd <- hexify(df, lon = "lon", lat = "lat", grid = g)

  result <- hex_extract(r, hd)
  expect_true(is.data.frame(result))
  expect_true("cell_id" %in% names(result))
  expect_true("elevation" %in% names(result))
})

test_that("hex_extract errors without terra", {
  skip_if(requireNamespace("terra", quietly = TRUE))
  expect_error(hex_extract(NULL, hex_grid(area_km2 = 1000)), "terra")
})

test_that("hex_zonal works with synthetic SpatRaster", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 100, ncols = 100,
                    xmin = -10, xmax = 10, ymin = 40, ymax = 55)
  terra::values(r) <- runif(10000)
  names(r) <- "temperature"

  df <- data.frame(lon = c(0, 5), lat = c(45, 50))
  g <- hex_grid(area_km2 = 2000)
  hd <- hexify(df, lon = "lon", lat = "lat", grid = g)

  result <- hex_zonal(r, hd, fun = "mean")
  expect_true(is.data.frame(result))
  expect_true("cell_id" %in% names(result))
  expect_true("temperature" %in% names(result))
})
