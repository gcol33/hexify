# tests/testthat/test-grid-helpers.R
# Tests for grid helper functions

test_that("lonlat_to_cell works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
  expect_true(all(cells > 0))
})

test_that("lonlat_to_cell works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Use the HexData object as grid source
  cells <- lonlat_to_cell(lon = 5, lat = 48, grid = result)

  expect_type(cells, "double")
  expect_length(cells, 1)
})

test_that("lonlat_to_cell works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
})

test_that("cell_to_lonlat returns cell centers", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  expect_length(coords$lon_deg, 2)
})

test_that("cell_to_lonlat works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_length(coords$lon_deg, 2)
})

test_that("lonlat_to_cell -> cell_to_lonlat round-trip lands in same cell", {
  grid <- hex_grid(area_km2 = 1000)

  original_lon <- c(0, 10, -5)
  original_lat <- c(45, 50, 48)

  cells1 <- lonlat_to_cell(original_lon, original_lat, grid)
  centers <- cell_to_lonlat(cells1, grid)
  cells2 <- lonlat_to_cell(centers$lon_deg, centers$lat_deg, grid)

  expect_equal(cells1, cells2)
})

test_that("cell_to_sf creates sf polygons", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  polys <- cell_to_sf(cells, grid)

  expect_s3_class(polys, "sf")
  expect_true("cell_id" %in% names(polys))
  expect_equal(nrow(polys), length(unique(cells)))
})

test_that("cell_to_sf works with HexData (no cell_id)", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  polys <- cell_to_sf(grid = result)

  expect_s3_class(polys, "sf")
  expect_equal(nrow(polys), length(unique(result@cell_id)))
})

test_that("cell_to_sf errors without cell_id for HexGridInfo", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(grid = grid), "cell_id required")
})

test_that("cell_to_sf errors on empty cell_id", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(cell_id = numeric(0), grid = grid), "No valid")
  expect_error(cell_to_sf(cell_id = c(NA, NA), grid = grid), "No valid")
})

test_that("grid_rect generates regional grid", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 50000)
  europe <- grid_rect(c(-10, 35, 30, 60), grid)

  expect_s3_class(europe, "sf")
  expect_gt(nrow(europe), 0)
})

test_that("grid_rect accepts sf bbox input", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 50000)
  france <- hexify_world[hexify_world$name == "France", ]
  hexes <- grid_rect(france, grid)

  expect_s3_class(hexes, "sf")
  expect_gt(nrow(hexes), 0)
})

test_that("grid_global generates global grid", {
  skip_if_not_installed("sf")

  # Use very coarse grid to keep test fast
  grid <- hex_grid(area_km2 = 500000)
  global <- grid_global(grid)

  expect_s3_class(global, "sf")
  expect_gt(nrow(global), 10)  # Should have at least some cells
})

test_that("grid_global warns for large grids", {
  skip_if_not_installed("sf")

  # Small cells = many cells = warning
  grid <- hex_grid(area_km2 = 100)  # Very small cells
  expect_warning(grid_global(grid), "cells")
})

test_that("grid_clip clips to boundary", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 20000)
  france <- hexify_world[hexify_world$name == "France", ]
  clipped <- grid_clip(france, grid)

  expect_s3_class(clipped, "sf")
  expect_gt(nrow(clipped), 0)
})

test_that("grid_clip with crop = FALSE keeps full hexagons", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 20000)
  france <- hexify_world[hexify_world$name == "France", ]
  clipped <- grid_clip(france, grid, crop = FALSE)

  expect_s3_class(clipped, "sf")
})

test_that("extract_grid works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  g <- hexify:::extract_grid(result)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with legacy hexify_grid", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid errors on invalid input", {
  expect_error(hexify:::extract_grid(list(a = 1)), "Cannot extract grid")
  expect_error(hexify:::extract_grid(data.frame()), "Cannot extract grid")
})
