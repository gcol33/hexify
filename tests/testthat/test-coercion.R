
# tests/testthat/test-coercion.R
# Tests for type coercion functions (as_tibble, as_sf)

# =============================================================================
# as_sf Tests
# =============================================================================

test_that("as_sf.HexData creates sf with point geometry by default", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  sf_obj <- as_sf(result)

  expect_s3_class(sf_obj, "sf")
  # Check geometry type is POINT
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true("POINT" %in% geom_types)
})

test_that("as_sf.HexData creates sf with polygon geometry", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55), value = 1:3)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  sf_obj <- as_sf(result, geometry = "polygon")

  expect_s3_class(sf_obj, "sf")
  # Check geometry type is POLYGON
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true(any(geom_types %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("as_sf.HexData preserves data columns", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(100, 200))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  sf_obj <- as_sf(result)

  expect_true("value" %in% names(sf_obj))
  expect_equal(sf_obj$value, c(100, 200))
})

test_that("as_sf.HexData works with sf input data", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  sf_obj <- as_sf(result)

  expect_s3_class(sf_obj, "sf")
})

test_that("as_sf.HexData polygon geometry works with sf input", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 10000)

  sf_obj <- as_sf(result, geometry = "polygon")

  expect_s3_class(sf_obj, "sf")
  geom_types <- unique(sf::st_geometry_type(sf_obj))
  expect_true(any(geom_types %in% c("POLYGON", "MULTIPOLYGON")))
})

test_that("as_sf.default errors for non-HexData", {
  expect_error(as_sf(data.frame()), "not defined for objects")
  expect_error(as_sf(list()), "not defined for objects")
})

test_that("as_sf requires sf package", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # This test just confirms the function runs - package is installed
  sf_obj <- as_sf(result)
  expect_s3_class(sf_obj, "sf")
})

# =============================================================================
# as_tibble Tests
# =============================================================================

test_that("as_tibble.HexData creates tibble", {
  skip_if_not_installed("tibble")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  tbl <- as_tibble.HexData(result)

  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 2)
})

test_that("as_tibble.HexData preserves columns", {
  skip_if_not_installed("tibble")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(100, 200))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  tbl <- as_tibble.HexData(result)

  expect_true("value" %in% names(tbl))
  expect_equal(tbl$value, c(100, 200))
})

test_that("as_tibble.HexData drops sf geometry", {
  skip_if_not_installed("tibble")
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  tbl <- as_tibble.HexData(result)

  expect_s3_class(tbl, "tbl_df")
  expect_false(inherits(tbl, "sf"))
})
