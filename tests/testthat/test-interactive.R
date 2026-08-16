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

  set.seed(1)
  df <- data.frame(
    lon = runif(10, 5, 15),
    lat = runif(10, 45, 55),
    temp = rnorm(10, 15)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  m <- hex_browse(hd, value = "temp")
  expect_true(inherits(m, "leaflet"))

  # The value= mapping must actually reach the widget's polygon fill colors
  # (not just build *a* map): fillColor should vary across cells, and a map
  # built without value= should use the flat default fill instead. The map
  # carries one polygon per cell, while the hexified data carries one row per
  # point, and several points can share a cell.
  poly_call <- Filter(function(c) c$method == "addPolygons", m$x$calls)[[1]]
  fill_colors <- poly_call$args[[4]]$fillColor
  expect_length(fill_colors, length(unique(as.data.frame(hd)$cell_id)))
  expect_gt(length(unique(fill_colors)), 1)

  m_default <- hex_browse(hd)
  poly_call_default <- Filter(function(c) c$method == "addPolygons", m_default$x$calls)[[1]]
  expect_identical(poly_call_default$args[[4]]$fillColor, "#3388ff")
})

test_that("hex_browse errors without leaflet", {
  skip_if(requireNamespace("leaflet", quietly = TRUE))
  df <- data.frame(lon = 10, lat = 50)
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  expect_error(hex_browse(hd), "leaflet")
})
