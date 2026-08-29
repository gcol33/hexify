
# tests/testthat/test-plot-methods.R
# Tests for plot methods

test_that("plot.HexData works with default settings", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Should not error
  expect_silent({
    pdf(NULL)  # Null device
    on.exit(dev.off())
    plot(result, basemap = FALSE)
  })
})

test_that("jitter_points_in_cells(jitter = FALSE) snaps points to the cell centroid", {
  skip_if_not_installed("sf")

  # Many points in one cell
  df <- data.frame(lon = rep(2.35, 20), lat = rep(48.86, 20))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 200000)
  hex_sf <- cell_to_sf(unique(result@cell_id), result@grid)

  jittered <- hexify:::jitter_points_in_cells(result@cell_id, hex_sf, jitter = TRUE)
  unjittered <- hexify:::jitter_points_in_cells(result@cell_id, hex_sf, jitter = FALSE)

  # With jitter, 20 random points in the same cell are extremely unlikely
  # to all land on the exact same coordinate
  expect_true(length(unique(jittered$lon)) > 1 || length(unique(jittered$lat)) > 1)

  # Without jitter, every point in the cell collapses to the same centroid
  expect_equal(length(unique(unjittered$lon)), 1)
  expect_equal(length(unique(unjittered$lat)), 1)
})

test_that("plot.HexData accepts jitter and point_alpha arguments", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE, jitter = FALSE,
         point_alpha = 0.5, point_color = "darkblue")
  })
})

test_that("plot.HexData works with basemap", {
  skip_on_cran()  # Slow sf operations
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # May produce sf messages about S2, so just check it runs
  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    suppressMessages(plot(result, basemap = TRUE))
  })
})

test_that("plot.HexData works with custom sf basemap", {
  skip_on_cran()  # Slow sf operations
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Use subset of world as custom basemap
  custom_basemap <- hexify_world[hexify_world$continent == "Europe", ]

  # May produce sf messages about S2, so just check it runs
  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    suppressMessages(plot(result, basemap = custom_basemap))
  })
})

test_that("plot.HexData works with an sfc basemap", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(2.35, 4.90), lat = c(48.86, 52.37))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  # Geometry column only: no rows to count when clipping
  custom_basemap <- sf::st_geometry(hexify_world)

  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    suppressMessages(plot(result, basemap = custom_basemap))
  })
})

test_that("plot.HexData works with a SpatRaster basemap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  df <- data.frame(lon = c(2.35, 4.90), lat = c(48.86, 52.37))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  r <- terra::rast(nrows = 20, ncols = 20,
                   xmin = -10, xmax = 15, ymin = 40, ymax = 55)
  terra::values(r) <- seq_len(400)

  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    suppressMessages(plot(result, basemap = r))
  })
})

test_that("plot.HexData works with the world_hires basemap", {
  skip_if_not_installed("sf")
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("rnaturalearthdata")

  df <- data.frame(lon = c(2.35, 4.90), lat = c(48.86, 52.37))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_no_error({
    pdf(NULL)
    on.exit(dev.off())
    suppressMessages(plot(result, basemap = "world_hires"))
  })
})

test_that("plot.HexData works with show_points", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 2.36, 2.37),
    lat = c(48.86, 48.87, 48.88)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE)
  })
})

test_that("plot.HexData works with point_size presets", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 2.36),
    lat = c(48.86, 48.87)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  for (size in c("tiny", "small", "normal", "large", "auto")) {
    expect_silent({
      pdf(NULL)
      on.exit(dev.off(), add = TRUE)
      plot(result, basemap = FALSE, show_points = TRUE, point_size = size)
    })
  }
})

test_that("plot.HexData works with numeric point_size", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 2.36),
    lat = c(48.86, 48.87)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE, point_size = 0.5)
  })
})

test_that("plot.HexData works with fill mapping", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    value = c(10, 20, 30)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, fill = "value")
  })
})

test_that("plot.HexData works with discrete fill", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90, -3.70),
    lat = c(48.86, 52.37, 40.42),
    category = c("A", "B", "A")
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, fill = "category")
  })
})

test_that("plot.HexData errors on invalid fill column", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_error({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, fill = "nonexistent")
  }, "not found")
})

test_that("plot.HexData works with crop = FALSE", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, crop = FALSE)
  })
})

test_that("plot.HexData works with clip_basemap = FALSE", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = TRUE, clip_basemap = FALSE)
  })
})

test_that("plot.HexData works with custom colors", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result,
         basemap = FALSE,
         grid_fill = "lightblue",
         grid_border = "darkblue",
         grid_alpha = 0.5,
         grid_lwd = 1.5)
  })
})

test_that("plot.HexData returns HexData invisibly", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 4.90),
    lat = c(48.86, 52.37)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  pdf(NULL)
  on.exit(dev.off())
  returned <- plot(result, basemap = FALSE)
  expect_s4_class(returned, "HexData")
})

test_that("plot.HexData works with title", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, main = "Test Title")
  })
})

test_that("plot.HexData warns on unknown point_size", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(2.35), lat = c(48.86))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  expect_warning({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE, point_size = "unknown")
  }, "Unknown point_size")
})

test_that("plot_grid creates ggplot", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]

  p <- suppressMessages(plot_grid(france, grid))
  expect_s3_class(p, "ggplot")
})

test_that("plot_grid works with custom colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]

  p <- suppressMessages(plot_grid(france, grid,
                                   grid_fill = "coral",
                                   boundary_fill = "lightyellow"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_grid works with crop = FALSE", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]

  p <- suppressMessages(plot_grid(france, grid, crop = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("plot_grid auto-generates title", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]

  p <- suppressMessages(plot_grid(france, grid, title = NULL))
  expect_s3_class(p, "ggplot")
})

test_that("plot_grid accepts custom title", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]

  p <- suppressMessages(plot_grid(france, grid, title = "Custom Title"))
  expect_s3_class(p, "ggplot")
})

test_that("plot.HexData works with 'very large' point_size", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 2.36),
    lat = c(48.86, 48.87)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE, point_size = "very large")
  })
})

test_that("plot.HexData works with 'verylarge' point_size variant", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(2.35, 2.36),
    lat = c(48.86, 48.87)
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  expect_silent({
    pdf(NULL)
    on.exit(dev.off())
    plot(result, basemap = FALSE, show_points = TRUE, point_size = "verylarge")
  })
})
