
# tests/testthat/test-visualization.R
# Tests for visualization functions
#
# Functions tested:
# - plot() method for HexData
# - hexify_heatmap()
# - plot_world()
# - Internal helper functions

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

test_that("calculate_view_buffer returns valid buffer", {
  bbox <- c(xmin = 0, xmax = 10, ymin = 40, ymax = 50)

  buffer <- hexify:::calculate_view_buffer(bbox)

  expect_type(buffer, "list")
  expect_true("x" %in% names(buffer))
  expect_true("y" %in% names(buffer))
  expect_true(buffer$x > 0)
  expect_true(buffer$y > 0)
})

test_that("calculate_view_buffer respects factor parameter", {
  bbox <- c(xmin = 0, xmax = 100, ymin = 0, ymax = 100)

  buffer_10 <- hexify:::calculate_view_buffer(bbox, factor = 0.1)
  buffer_20 <- hexify:::calculate_view_buffer(bbox, factor = 0.2)

  expect_true(buffer_20$x > buffer_10$x)
  expect_true(buffer_20$y > buffer_10$y)
})

test_that("calculate_view_buffer respects min_buffer", {
  bbox <- c(xmin = 0, xmax = 1, ymin = 0, ymax = 1)

  buffer <- hexify:::calculate_view_buffer(bbox, factor = 0.1, min_buffer = 5)

  expect_gte(buffer$x, 5)
  expect_gte(buffer$y, 5)
})

test_that("is_palette_name identifies palette names correctly", {
  expect_true(hexify:::is_palette_name("YlOrRd"))
  expect_true(hexify:::is_palette_name("viridis"))

  expect_false(hexify:::is_palette_name("#FF0000"))
  expect_false(hexify:::is_palette_name("red"))
  expect_false(hexify:::is_palette_name(c("red", "blue")))
})

test_that("is_brewer_palette identifies valid palettes", {
  skip_if_not_installed("RColorBrewer")

  expect_true(hexify:::is_brewer_palette("YlOrRd"))
  expect_true(hexify:::is_brewer_palette("Blues"))
  expect_true(hexify:::is_brewer_palette("Set1"))

  expect_false(hexify:::is_brewer_palette("NotAPalette"))
  expect_false(hexify:::is_brewer_palette("viridis"))
})

test_that("generate_bin_labels creates correct labels", {
  breaks <- c(0, 10, 20, 30)
  labels <- hexify:::generate_bin_labels(breaks)

  expect_length(labels, 3)
  expect_equal(labels, c("0-10", "10-20", "20-30"))
})

test_that("generate_bin_labels handles infinite bounds", {
  breaks <- c(-Inf, 10, 100, Inf)
  labels <- hexify:::generate_bin_labels(breaks)

  expect_length(labels, 3)
  expect_equal(labels[1], "<10")
  expect_equal(labels[3], ">100")
})

test_that("resolve_value_column auto-detects count column", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$count <- c(10, 20)

  col <- hexify:::resolve_value_column(hex_sf, NULL)
  expect_equal(col, "count")
})

test_that("resolve_value_column auto-detects n column", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$n <- c(10, 20)

  col <- hexify:::resolve_value_column(hex_sf, NULL)
  expect_equal(col, "n")
})

test_that("resolve_value_column returns specified column", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  col <- hexify:::resolve_value_column(hex_sf, "value")
  expect_equal(col, "value")
})

test_that("resolve_value_column errors on missing column", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)

  expect_error(
    hexify:::resolve_value_column(hex_sf, "nonexistent"),
    "not found in data"
  )
})

test_that("resolve_value_column returns NULL when no suitable column found", {
  skip_if_not_installed("sf")


  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)

  # With require = FALSE (default), returns NULL for uniform fill
  result <- hexify:::resolve_value_column(hex_sf, NULL)
  expect_null(result)

  # With require = TRUE, throws error
  expect_error(
    hexify:::resolve_value_column(hex_sf, NULL, require = TRUE),
    "No 'value' column specified"
  )
})

test_that("prepare_fill_column handles discrete data", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$category <- factor(c("A", "B"))

  result <- hexify:::prepare_fill_column(hex_sf, "category", NULL, NULL)

  expect_equal(result$fill_col, "category")
  expect_true(result$is_discrete)
})

test_that("prepare_fill_column handles continuous data without breaks", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  result <- hexify:::prepare_fill_column(hex_sf, "value", NULL, NULL)

  expect_equal(result$fill_col, "value")
  expect_false(result$is_discrete)
})

test_that("prepare_fill_column applies breaks to continuous data", {
  skip_if_not_installed("sf")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(5, 15)

  result <- hexify:::prepare_fill_column(
    hex_sf, "value", c(0, 10, 20), c("Low", "High")
  )

  expect_true("value_bin" %in% names(result$data))
  expect_equal(result$fill_col, "value_bin")
  expect_true(result$is_discrete)
})

test_that("resolve_basemap returns NULL for NULL input", {
  result <- hexify:::resolve_basemap(NULL)
  expect_null(result)
})

test_that("resolve_basemap returns hexify_world for 'world'", {
  skip_if_not_installed("sf")

  result <- hexify:::resolve_basemap("world")
  expect_s3_class(result, "sf")
})

test_that("resolve_basemap returns sf for sf input", {
  skip_if_not_installed("sf")

  result <- hexify:::resolve_basemap(hexify_world)
  expect_s3_class(result, "sf")
})

test_that("resolve_basemap errors on invalid input", {
  expect_error(hexify:::resolve_basemap(123), "basemap must be")
  expect_error(hexify:::resolve_basemap("invalid"), "basemap must be")
})

test_that("resolve_basemap_with_raster returns correct structure", {
  result <- hexify:::resolve_basemap_with_raster(NULL)

  expect_type(result, "list")
  expect_true("sf" %in% names(result))
  expect_true("raster" %in% names(result))
  expect_null(result$sf)
  expect_null(result$raster)
})

test_that("resolve_basemap_with_raster handles world basemap", {
  result <- hexify:::resolve_basemap_with_raster("world")

  expect_s3_class(result$sf, "sf")
  expect_null(result$raster)
})

test_that("resolve_basemap_with_raster handles sf input", {
  skip_if_not_installed("sf")

  result <- hexify:::resolve_basemap_with_raster(hexify_world)

  expect_s3_class(result$sf, "sf")
  expect_null(result$raster)
})

test_that("resolve_basemap_with_raster errors on invalid input", {
  expect_error(hexify:::resolve_basemap_with_raster(123), "basemap must be")
})

test_that("prepare_hex_sf_simple validates input", {
  expect_error(
    hexify:::prepare_hex_sf_simple(data.frame(x = 1), aperture = 3),
    "HexData object or an sf object"
  )

  expect_error(
    hexify:::prepare_hex_sf_simple(data.frame(cell_id = 1), aperture = 3),
    "cell_area"
  )
})

test_that("prepare_hex_sf_simple passes through sf objects", {
  skip_if_not_installed("sf")

  hex_sf <- hexify_cell_to_sf(c(12847), resolution = 10, aperture = 3)
  result <- hexify:::prepare_hex_sf_simple(hex_sf, aperture = 3)

  expect_s3_class(result, "sf")
})

test_that("prepare_hex_sf validates input", {
  expect_error(
    hexify:::prepare_hex_sf(data.frame(x = 1), aperture = 3),
    "HexData object or an sf object"
  )

  expect_error(
    hexify:::prepare_hex_sf(data.frame(cell_id = 1), aperture = 3),
    "cell_area"
  )
})

test_that("prepare_hex_sf passes through sf objects", {
  skip_if_not_installed("sf")

  hex_sf <- hexify_cell_to_sf(c(12847), resolution = 10, aperture = 3)
  result <- hexify:::prepare_hex_sf(hex_sf, aperture = 3)

  expect_s3_class(result, "sf")
})

test_that("prepare_hex_sf merges extra columns", {
  skip_if_not_installed("sf")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    extra_col = c("A", "B")
  )

  result <- hexify:::prepare_hex_sf(df, aperture = 3)

  expect_s3_class(result, "sf")
  expect_true("extra_col" %in% names(result))
})

# =============================================================================
# HEXIFY_WORLD DATA
# =============================================================================

test_that("hexify_world data is available and valid", {
  expect_true(exists("hexify_world"))
  expect_s3_class(hexify_world, "sf")
  expect_true(nrow(hexify_world) > 100)
  expect_true("name" %in% names(hexify_world))
  expect_true("continent" %in% names(hexify_world))
  expect_equal(sf::st_crs(hexify_world)$epsg, 4326)
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
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap auto-detects count column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap auto-detects n column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    n = c(10, 20)
  )

  result <- hexify_heatmap(df)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works without value column (uniform fill)", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94)
  )

  # Should work with uniform fill when no value column specified
  result <- hexify_heatmap(df)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with basemap", {
  skip_on_cran()  # Slow sf operations
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", basemap = "world")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with breaks", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
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
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", colors = c("white", "red"))
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap respects styling parameters", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
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
  skip_on_cran()  # Slow CRS transformation
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
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
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  expect_error(hexify_heatmap(df, value = "nonexistent"), "not found in data")
  expect_error(hexify_heatmap(df, basemap = 123), "basemap must be")
})

test_that("hexify_heatmap works with mask_outside", {
  skip_on_cran()  # Slow sf intersection operations
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", basemap = "world",
                           mask_outside = TRUE)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with discrete factor column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    category = factor(c("A", "B"))
  )

  result <- hexify_heatmap(df, value = "category")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with character column", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    type = c("type_A", "type_B"),
    stringsAsFactors = FALSE
  )

  result <- hexify_heatmap(df, value = "type")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with RColorBrewer palette", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", colors = "YlOrRd")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with RColorBrewer discrete scale", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    category = factor(c("A", "B"))
  )

  result <- hexify_heatmap(df, value = "category", colors = "Set1")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap works with viridis option", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", colors = "magma")
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap auto-generates break labels", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 200)
  )

  result <- hexify_heatmap(df, value = "count",
                           breaks = c(0, 50, 100, 200))
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap respects theme_void=FALSE", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count", theme_void = FALSE)
  expect_s3_class(result, "ggplot")
})

test_that("hexify_heatmap respects xlim and ylim", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    cell_id = c(12847, 12532),
    cell_area = c(863.94, 863.94),
    count = c(10, 20)
  )

  result <- hexify_heatmap(df, value = "count",
                           xlim = c(-10, 20), ylim = c(40, 60))
  expect_s3_class(result, "ggplot")
})

# =============================================================================
# COLOR SCALE HELPER FUNCTIONS
# =============================================================================

test_that("apply_discrete_scale works with NULL colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$category <- factor(c("A", "B"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = category))

  result <- hexify:::apply_discrete_scale(p, NULL, "Category", "gray90", 2)
  expect_s3_class(result, "ggplot")
})

test_that("apply_discrete_scale works with manual colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$category <- factor(c("A", "B"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = category))

  result <- hexify:::apply_discrete_scale(
    p, c("red", "blue"), "Category", "gray90", 2
  )
  expect_s3_class(result, "ggplot")
})

test_that("apply_continuous_scale works with NULL colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = value))

  result <- hexify:::apply_continuous_scale(p, NULL, "Value", "gray90")
  expect_s3_class(result, "ggplot")
})

test_that("apply_continuous_scale works with manual colors", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = value))

  result <- hexify:::apply_continuous_scale(
    p, c("white", "red"), "Value", "gray90"
  )
  expect_s3_class(result, "ggplot")
})

test_that("apply_continuous_scale works with Brewer palette", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = value))

  result <- hexify:::apply_continuous_scale(p, "YlOrRd", "Value", "gray90")
  expect_s3_class(result, "ggplot")
})

test_that("apply_discrete_scale falls back to viridis for unknown palette", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$category <- factor(c("A", "B"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = category))

  # Use an unknown palette name that's not a Brewer palette
  result <- hexify:::apply_discrete_scale(
    p, "inferno", "Category", "gray90", 2
  )
  expect_s3_class(result, "ggplot")
})

test_that("apply_continuous_scale falls back to viridis for unknown palette", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  cell_ids <- c(12847, 12532)
  hex_sf <- hexify_cell_to_sf(cell_ids, resolution = 10, aperture = 3)
  hex_sf$value <- c(10, 20)

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = hex_sf, ggplot2::aes(fill = value))

  # Use an unknown palette name
  result <- hexify:::apply_continuous_scale(p, "cividis", "Value", "gray90")
  expect_s3_class(result, "ggplot")
})

test_that("resolve_basemap handles world_hires without package", {
  # Skip if rnaturalearth is installed - we want to test the error path
  skip_if(requireNamespace("rnaturalearth", quietly = TRUE),
          "rnaturalearth is installed, skipping error test")

  expect_error(
    hexify:::resolve_basemap("world_hires"),
    "rnaturalearth"
  )
})

test_that("resolve_basemap_with_raster handles sfc input", {
  skip_if_not_installed("sf")

  sfc <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 4326)
  result <- hexify:::resolve_basemap_with_raster(sfc)

  expect_s3_class(result$sf, "sfc")
  expect_null(result$raster)
})

test_that("hexify_heatmap handles data with NA CRS", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")

  # Create sf object and remove CRS
  hex_sf <- hexify_cell_to_sf(c(12847, 12532), resolution = 10, aperture = 3)
  hex_sf$count <- c(10, 20)
  sf::st_crs(hex_sf) <- NA

  result <- hexify_heatmap(hex_sf, value = "count")
  expect_s3_class(result, "ggplot")
})
