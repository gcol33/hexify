# tests/testthat/test-coercion.R
# Tests for type coercion functions (as_tibble, st_as_sf)

# =============================================================================
# st_as_sf Tests
# =============================================================================

test_that("st_as_sf.HexData creates sf with point geometry by default", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  sf_obj <- st_as_sf(result)

  expect_s3_class(sf_obj, "sf")
  # Check geometry type is POINT
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true("POINT" %in% geom_types)
})

test_that("st_as_sf.HexData creates sf with polygon geometry", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55), value = 1:3)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  sf_obj <- st_as_sf(result, geometry = "polygon")

  expect_s3_class(sf_obj, "sf")
  # Check geometry type is POLYGON
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true(any(geom_types %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("st_as_sf.HexData preserves data columns", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(100, 200))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  sf_obj <- st_as_sf(result)

  expect_true("value" %in% names(sf_obj))
  expect_equal(sf_obj$value, c(100, 200))
})

test_that("st_as_sf.HexData works with sf input data", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  sf_obj <- st_as_sf(result)

  expect_s3_class(sf_obj, "sf")
})

test_that("st_as_sf.HexData polygon geometry works with sf input", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 10000)

  sf_obj <- st_as_sf(result, geometry = "polygon")

  expect_s3_class(sf_obj, "sf")
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true(any(geom_types %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("sf::st_as_sf() dispatches on hexify objects", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Reached through sf's own generic, without hexify attached to the call
  expect_s3_class(sf::st_as_sf(result), "sf")
  expect_s3_class(sf::st_as_sf(result, geometry = "polygon"), "sf")
})

test_that("st_as_sf.HexGridInfo returns the global cell set", {
  skip_if_not_installed("sf")

  grid <- hex_grid(resolution = 2)
  cells <- st_as_sf(grid)

  expect_s3_class(cells, "sf")
  expect_equal(nrow(cells), 2 + 10 * 3^2)
  expect_equal(nrow(cells), nrow(grid_global(grid)))
})

# =============================================================================
# as_tibble Tests
# =============================================================================

test_that("tibble::as_tibble() dispatches on HexData", {
  skip_if_not_installed("tibble")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  tbl <- tibble::as_tibble(result)

  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 2)
})

test_that("as_tibble.HexData preserves columns", {
  skip_if_not_installed("tibble")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(100, 200))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  tbl <- tibble::as_tibble(result)

  expect_true("value" %in% names(tbl))
  expect_equal(tbl$value, c(100, 200))
})

test_that("as_tibble.HexData drops sf geometry", {
  skip_if_not_installed("tibble")
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  tbl <- tibble::as_tibble(result)

  expect_s3_class(tbl, "tbl_df")
  expect_false(inherits(tbl, "sf"))
})
