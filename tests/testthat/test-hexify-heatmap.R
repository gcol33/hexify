
# tests/testthat/test-hexify-heatmap.R
# Tests for hexify_heatmap function

test_that("hexify_heatmap creates ggplot object", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    value = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p$layers[[1]]$geom, "GeomSf")
  expect_true(nrow(p$layers[[1]]$data) > 0)
})

test_that("hexify_heatmap works with value column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    count = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, value = "count")

  expect_s3_class(p, "ggplot")
  # The fill aesthetic must actually reference the requested column, not just
  # produce *a* plot.
  expect_identical(rlang::as_label(p$layers[[1]]$mapping$fill), "count")
})

test_that("hexify_heatmap works with basemap = 'world'", {
  skip_on_cran()  # Slow sf intersection operations
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- suppressMessages(hexify_heatmap(result, basemap = "world"))

  expect_s3_class(p, "ggplot")
})

test_that("hexify_heatmap works with custom sf basemap", {
  skip_on_cran()  # Slow sf operations
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  custom_basemap <- hexify_world[hexify_world$continent == "Europe", ]
  p <- suppressMessages(hexify_heatmap(result, basemap = custom_basemap))

  expect_s3_class(p, "ggplot")
})

test_that("hexify_heatmap works with discrete values", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    category = c("A", "B", "A")
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, value = "category")

  expect_s3_class(p, "ggplot")
  expect_identical(rlang::as_label(p$layers[[1]]$mapping$fill), "category")
  # A discrete value column must produce a discrete (not continuous) scale.
  fill_scale <- Filter(function(s) identical(s$aesthetics, "fill"), p$scales$scales)[[1]]
  expect_s3_class(fill_scale, "ScaleDiscrete")
})

test_that("hexify_heatmap works with custom colors vector", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    value = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, value = "value",
                      colors = c("blue", "white", "red"))

  expect_s3_class(p, "ggplot")
  fill_scale <- Filter(function(s) identical(s$aesthetics, "fill"), p$scales$scales)[[1]]
  # The gradient's endpoints must actually be the requested colors, not the
  # default viridis palette.
  ends <- fill_scale$palette(c(0, 1))
  expect_equal(toupper(ends), c("#0000FF", "#FF0000"))
})

test_that("hexify_heatmap works with breaks", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    value = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, value = "value",
                      breaks = c(-Inf, 15, 25, Inf),
                      labels = c("Low", "Medium", "High"))

  expect_s3_class(p, "ggplot")
  expect_identical(rlang::as_label(p$layers[[1]]$mapping$fill), "value_bin")
  expect_identical(levels(p$layers[[1]]$data[["value_bin"]]), c("Low", "Medium", "High"))
  # Bin membership must match the supplied breaks, not just have 3 levels.
  expect_identical(
    as.character(p$layers[[1]]$data[["value_bin"]][order(p$layers[[1]]$data$value)]),
    c("Low", "Medium", "High")
  )
})

test_that("hexify_heatmap works with breaks (auto labels)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    value = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Auto-generated labels
  p <- hexify_heatmap(result, value = "value",
                      breaks = c(-Inf, 15, 25, Inf))

  expect_s3_class(p, "ggplot")
  bin_levels <- levels(p$layers[[1]]$data[["value_bin"]])
  expect_length(bin_levels, 3)
  # Auto-labels must be distinct and non-empty, not e.g. all "NA" or blank.
  expect_false(anyNA(bin_levels))
  expect_length(unique(bin_levels), 3)
})

test_that("hexify_heatmap works with xlim and ylim", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result,
                      xlim = c(-10, 20),
                      ylim = c(40, 60))

  expect_s3_class(p, "ggplot")
  expect_equal(p$coordinates$limits$x, c(-10, 20))
  expect_equal(p$coordinates$limits$y, c(40, 60))
})

test_that("hexify_heatmap works with title", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, title = "Test Title")

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Test Title")
})

test_that("hexify_heatmap works with legend_title", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37),
    value = c(10, 20)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, value = "value", legend_title = "Custom Legend")

  expect_s3_class(p, "ggplot")
  fill_scale <- Filter(function(s) identical(s$aesthetics, "fill"), p$scales$scales)[[1]]
  expect_identical(fill_scale$name, "Custom Legend")

  # Default (no legend_title) falls back to the value column name.
  p_default <- hexify_heatmap(result, value = "value")
  fill_scale_default <- Filter(function(s) identical(s$aesthetics, "fill"),
                                p_default$scales$scales)[[1]]
  expect_identical(fill_scale_default$name, "value")
})

test_that("hexify_heatmap works with mask_outside", {
  skip_on_cran()  # Slow sf intersection operations
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- suppressMessages(hexify_heatmap(result,
                                        basemap = "world",
                                        mask_outside = TRUE))
  p_unmasked <- suppressMessages(hexify_heatmap(result, basemap = "world"))

  expect_s3_class(p, "ggplot")
  # Masking to land must actually change the hex layer's geometry, not just
  # produce a plot with the same unmasked data.
  expect_false(isTRUE(all.equal(
    sf::st_area(p$layers[[1]]$data),
    sf::st_area(p_unmasked$layers[[1]]$data)
  )))
})

test_that("hexify_heatmap works with theme_void = FALSE", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result, theme_void = FALSE)
  p_void <- hexify_heatmap(result, theme_void = TRUE)

  expect_s3_class(p, "ggplot")
  # The two themes must actually differ, not both fall back to one default.
  expect_false(identical(p$theme$axis.text.x, p_void$theme$axis.text.x))
})

test_that("hexify_heatmap errors on invalid value column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_error(hexify_heatmap(result, value = "nonexistent"), "not found")
})

test_that("hexify_heatmap errors on invalid basemap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_error(hexify_heatmap(result, basemap = "invalid"), "must be")
})

test_that("hexify_heatmap works with custom styling", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(lon = c(2.35, 4.90), lat = c(48.86, 52.37))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  p <- hexify_heatmap(result,
                      hex_border = "blue",
                      hex_lwd = 1.5,
                      hex_alpha = 0.5)

  expect_s3_class(p, "ggplot")
  params <- p$layers[[1]]$aes_params
  expect_identical(params$colour, "blue")
  expect_equal(params$linewidth, 1.5)
  expect_equal(params$alpha, 0.5)
})

test_that("hexify_heatmap works with CRS transformation", {
  skip_on_cran()  # Slow CRS transformation
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # LAEA Europe projection
  p <- suppressWarnings(hexify_heatmap(result, crs = 3035))

  expect_s3_class(p, "ggplot")
  expect_equal(sf::st_crs(p$layers[[1]]$data)$epsg, 3035)
})

test_that("hexify_heatmap auto-detects count column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37),
    count = c(10, 20)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Should auto-detect 'count' column
  p <- hexify_heatmap(result)

  expect_s3_class(p, "ggplot")
  expect_identical(rlang::as_label(p$layers[[1]]$mapping$fill), "count")
})

test_that("hexify_heatmap auto-detects n column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37),
    n = c(10, 20)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Should auto-detect 'n' column
  p <- hexify_heatmap(result)

  expect_s3_class(p, "ggplot")
  expect_identical(rlang::as_label(p$layers[[1]]$mapping$fill), "n")
})

# =============================================================================
# plot_world Tests
# =============================================================================

test_that("plot_world creates plot", {
  skip_if_not_installed("sf")

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    result <- plot_world()
  })
  expect_null(result)
})

test_that("plot_world accepts custom colors", {
  skip_if_not_installed("sf")

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot_world(fill = "lightblue", border = "darkblue")
  })
})
