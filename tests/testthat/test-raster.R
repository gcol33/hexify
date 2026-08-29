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

test_that("a boundary works on any body, and the cells come back on it", {
  skip_if_not_installed("terra")

  boundary <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = -10, ymin = 35, xmax = 30, ymax = 60), crs = 4326))
  raster <- terra::rast(nrows = 90, ncols = 180)
  terra::values(raster) <- seq_len(terra::ncell(raster))
  names(raster) <- "value"

  for (body in list(EARTH_RADIUS_KM, "mars", "moon")) {
    label <- paste("radius", body)
    grid <- hex_grid(area_km2 = 2e5, radius_km = body)

    clipped <- suppressMessages(grid_clip(boundary, grid))
    expect_gt(nrow(clipped), 0)
    expect_true(all(sf::st_is_valid(clipped)), info = label)
    expect_equal(sf::st_crs(clipped), st_crs(grid), info = label)

    zonal <- suppressMessages(hex_zonal(raster, grid, boundary = boundary))
    expect_equal(nrow(zonal), nrow(clipped), info = label)
    expect_true("value" %in% names(zonal), info = label)

    extracted <- suppressMessages(hex_extract(raster, grid, boundary = boundary))
    expect_equal(nrow(extracted), nrow(clipped), info = label)
    expect_true("value" %in% names(extracted), info = label)
  }
})

test_that("a lon/lat boundary on another body is read in the grid's CRS", {
  boundary <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = -10, ymin = 35, xmax = 30, ymax = 60), crs = 4326))
  grid <- hex_grid(area_km2 = 2e5, radius_km = "mars")

  expect_message(grid_clip(boundary, grid), "name the same region on both")
  expect_message(grid_clip(boundary, grid), "the boundary is in EPSG:4326")

  # Same angles, so the same cells as the grid's own CRS names
  on_body <- boundary
  on_body <- suppressWarnings(sf::st_set_crs(on_body, st_crs(grid)))
  expect_equal(suppressMessages(grid_clip(boundary, grid))$cell_id,
               grid_clip(on_body, grid)$cell_id)
})

test_that("a boundary that cannot reach the grid's body says so", {
  boundary <- sf::st_transform(
    sf::st_as_sfc(sf::st_bbox(c(xmin = -10, ymin = 35, xmax = 30, ymax = 60),
                              crs = 4326)),
    3035)
  grid <- hex_grid(area_km2 = 2e5, radius_km = "mars")

  expect_error(grid_clip(boundary, grid),
               "which are on different bodies", fixed = TRUE)

  # On Earth the same boundary reprojects and clips
  expect_gt(nrow(grid_clip(boundary, hex_grid(area_km2 = 2e5))), 0)
})

test_that("hexify() reads sf coordinates on the grid's own body", {
  points <- sf::st_as_sf(data.frame(lon = c(0, 10), lat = c(45, 50)),
                         coords = c("lon", "lat"), crs = 4326)

  on_mars <- suppressMessages(
    hexify(points, grid = hex_grid(area_km2 = 2e5, radius_km = "mars")))
  expect_equal(nrow(on_mars), 2L)
  expect_false(anyNA(on_mars@cell_id))

  expect_message(hexify(points, grid = hex_grid(area_km2 = 2e5, radius_km = "mars")),
                 "name the same region on both")

  # An Earth grid reprojects as before, without a word
  on_earth <- hex_grid(area_km2 = 2e5)
  expect_silent(hexify(sf::st_transform(points, 4258), grid = on_earth))
  expect_equal(nrow(hexify(sf::st_transform(points, 3035), grid = on_earth)), 2L)
})
