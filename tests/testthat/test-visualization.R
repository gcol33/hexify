# tests/testthat/test-visualization.R
# Tests for visualization functions
#
# Functions tested:
# - hexify_map()
# - hexify_heatmap()
# - plot_world()

# =============================================================================
# HEXIFY_MAP
# =============================================================================

test_that("hexify_world data is available and valid", {
  expect_true(exists("hexify_world"))
  expect_s3_class(hexify_world, "sf")
  expect_true(nrow(hexify_world) > 100)
  expect_true("name" %in% names(hexify_world))
  expect_true("continent" %in% names(hexify_world))
  expect_equal(sf::st_crs(hexify_world)$epsg, 4326)
})

test_that("hexify_map works with hexify output", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(0, 5, 10),
    lat = c(45, 46, 45)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_silent(hexify_map(result))
})

test_that("hexify_map works with sf polygon input", {
  skip_if_not_installed("sf")

  hex_ids <- c(12847, 12532)
  polys <- hexify_polygons(hex_ids, resolution = 10, aperture = 3)

  expect_silent(hexify_map(polys))
})

test_that("hexify_map works with built-in world basemap", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(0, 5, 10),
    lat = c(45, 46, 45)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_no_error(hexify_map(result, basemap = "world"))
})

test_that("hexify_map works with custom sf basemap", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(0, 5, 10),
    lat = c(45, 46, 45)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_no_error(hexify_map(result, basemap = hexify_world))
})

test_that("hexify_map respects xlim and ylim", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(0, 5, 10),
    lat = c(45, 46, 45)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

  expect_silent(hexify_map(result, xlim = c(-5, 15), ylim = c(40, 50)))
})

test_that("hexify_map validates input", {
  expect_error(hexify_map(data.frame(x = 1)), "data.frame from hexify")
  expect_error(hexify_map(data.frame(hex_id = 1)), "hex_area")
  expect_error(hexify_map(list()), "data.frame from hexify")
})

test_that("hexify_map rejects invalid basemap types", {
  skip_if_not_installed("sf")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94)
  )

  expect_error(hexify_map(df, basemap = 123), "basemap must be")
  expect_error(hexify_map(df, basemap = "invalid"), "basemap must be")
})

test_that("hexify_map custom colors work", {
  skip_if_not_installed("sf")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94)
  )

  expect_no_error(hexify_map(df,
                              fill = "red",
                              border = "darkred",
                              alpha = 0.5,
                              basemap = "world",
                              basemap_fill = "ivory",
                              basemap_border = "gray"))
})

test_that("hexify_map respects lwd parameters", {
  skip_if_not_installed("sf")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94)
  )

  expect_silent(hexify_map(df, lwd = 2))
  expect_no_error(hexify_map(df, basemap = "world", lwd = 0.5, basemap_lwd = 2))
})

# =============================================================================
# PLOT_WORLD
# =============================================================================

test_that("plot_world works", {
  skip_if_not_installed("sf")

  expect_silent(plot_world())
  expect_silent(plot_world(fill = "lightblue", border = "navy"))
})

# =============================================================================
# HEXIFY_HEATMAP
# =============================================================================

test_that("hexify_heatmap requires ggplot2", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap auto-detects count column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap auto-detects n column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    n = c(10, 20)
  )

  result <- hexify_heatmap(df)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap errors without value column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94)
  )

  expect_error(hexify_heatmap(df), "No 'value' column specified")
})

test_that("hexify_heatmap works with basemap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", basemap = "world")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with breaks", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 200)
  )

  result <- hexify_heatmap(df, value = "count",
                            breaks = c(-Inf, 50, 100, Inf),
                            labels = c("Low", "Medium", "High"))
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with custom colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", colors = c("white", "red"))
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap respects styling parameters", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count",
                            hex_border = "darkblue",
                            hex_lwd = 0.5,
                            hex_alpha = 0.8,
                            title = "Test Plot",
                            legend_title = "Count")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with projection", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", crs = 3035)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap validates input", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  expect_error(hexify_heatmap(list()), "data must be")
  expect_error(hexify_heatmap(data.frame(x = 1)), "data must be")

  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  expect_error(hexify_heatmap(df, value = "nonexistent"), "not found in data")
  expect_error(hexify_heatmap(df, basemap = 123), "basemap must be")
})
