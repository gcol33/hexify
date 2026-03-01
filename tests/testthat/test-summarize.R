# tests/testthat/test-summarize.R
# Tests for hex_summarize()

test_that("hex_summarize returns cell counts with no expressions", {
  df <- data.frame(lon = c(10, 10, 10.5, 20), lat = c(50, 50, 50.1, 55))
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  result <- hex_summarize(hd)

  expect_true(is.data.frame(result))
  expect_true("cell_id" %in% names(result))
  expect_true("n_points" %in% names(result))
  expect_true("cell_cen_lon" %in% names(result))
  expect_true("cell_cen_lat" %in% names(result))
  expect_true("cell_area_km2" %in% names(result))
  expect_equal(sum(result$n_points), nrow(df))
})

test_that("hex_summarize works with tidyeval expressions", {
  df <- data.frame(
    lon = runif(50, 5, 15),
    lat = runif(50, 45, 55),
    temp = rnorm(50, 15, 5)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  result <- hex_summarize(hd, mean_temp = mean(temp))

  expect_true("mean_temp" %in% names(result))
  expect_true(is.numeric(result$mean_temp))
})

test_that("hex_summarize works with .fns argument", {
  df <- data.frame(
    lon = runif(50, 5, 15),
    lat = runif(50, 45, 55),
    temp = rnorm(50, 15, 5)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  result <- hex_summarize(hd, .fns = list(mean_temp = ~mean(.x$temp)))

  expect_true("mean_temp" %in% names(result))
})

test_that("hex_summarize with geometry=TRUE returns sf", {
  skip_if_not_installed("sf")
  df <- data.frame(
    lon = runif(20, 5, 15),
    lat = runif(20, 45, 55)
  )
  hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
  result <- hex_summarize(hd, geometry = TRUE)

  expect_true(inherits(result, "sf"))
})

test_that("hex_summarize rejects non-HexData input", {
  expect_error(hex_summarize(data.frame(x = 1)), "HexData")
})
