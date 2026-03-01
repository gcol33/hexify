# tests/testthat/test-interactive.R
# Smoke tests for hex_browse()

test_that("hex_browse creates a leaflet map", {
  skip_if_not_installed("leaflet")

  df <- data.frame(
    lon = runif(10, 5, 15),
    lat = runif(10, 45, 55),
    value = rnorm(10)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  m <- hex_browse(hd)
  expect_true(inherits(m, "leaflet"))
})

test_that("hex_browse with value column works", {
  skip_if_not_installed("leaflet")

  df <- data.frame(
    lon = runif(10, 5, 15),
    lat = runif(10, 45, 55),
    temp = rnorm(10, 15)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  m <- hex_browse(hd, value = "temp")
  expect_true(inherits(m, "leaflet"))
})

test_that("hex_browse errors without leaflet", {
  skip_if(requireNamespace("leaflet", quietly = TRUE))
  df <- data.frame(lon = 10, lat = 50)
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  expect_error(hex_browse(hd), "leaflet")
})
