# tests/testthat/test-polygons.R
# Tests for polygon generation functions
#
# Functions tested:
# - hexify_polygons() / hexify_cell_to_sf()
# - hexify_to_polygons()
# - hexify_plot()
# - hexify_grid_rect()
# - hexify_grid_global()
# - hex_corners_to_sf()

# =============================================================================
# HEXIFY_POLYGONS / HEXIFY_CELL_TO_SF
# =============================================================================

test_that("hexify_polygons returns data frame with return_sf=FALSE", {
  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("hex_id", "lon", "lat", "order") %in% names(result)))
})

test_that("hexify_polygons returns 7 vertices per cell (closed polygon)", {
  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(nrow(result), length(hex_ids) * 7)
  expect_equal(unique(result$order), 1:7)
})

test_that("hexify_polygons returns valid coordinates", {
  hex_ids <- c(12847, 12532)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_true(all(result$lon >= -180 & result$lon <= 180))
  expect_true(all(result$lat >= -90 & result$lat <= 90))
})

test_that("hexify_polygons produces closed polygons", {
  hex_ids <- c(12847)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  first_vertex <- result[result$order == 1, c("lon", "lat")]
  last_vertex <- result[result$order == 7, c("lon", "lat")]

  expect_equal(first_vertex$lon, last_vertex$lon)
  expect_equal(first_vertex$lat, last_vertex$lat)
})

test_that("hexify_polygons returns sf object with return_sf=TRUE", {
  skip_if_not_installed("sf")

  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = TRUE)

  expect_s3_class(result, "sf")
  expect_true("hex_id" %in% names(result))
  expect_true("geometry" %in% names(result))
  expect_equal(sf::st_crs(result)$epsg, 4326)

  geom_types <- sf::st_geometry_type(result)
  expect_true(all(geom_types == "POLYGON"))
})

test_that("hexify_polygons removes duplicates", {
  hex_ids <- c(12847, 12847, 12532)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(length(unique(result$hex_id)), 2)
  expect_equal(nrow(result), 2 * 7)
})

test_that("hexify_polygons handles NA values", {
  hex_ids <- c(12847, NA, 12532)

  result <- hexify_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(length(unique(result$hex_id)), 2)
})

test_that("hexify_polygons validates aperture", {
  expect_error(
    hexify_polygons(c(12847), resolution = 10, aperture = 5),
    "aperture must be 3, 4, or 7"
  )
})

test_that("hexify_polygons validates resolution", {
  expect_error(
    hexify_polygons(c(12847), resolution = -1, aperture = 3),
    "resolution must be between 0 and 30"
  )

  expect_error(
    hexify_polygons(c(12847), resolution = 31, aperture = 3),
    "resolution must be between 0 and 30"
  )
})

test_that("hexify_polygons works with aperture 4", {
  hex_ids <- c(100, 200, 300)

  result <- hexify_polygons(hex_ids, resolution = 8, aperture = 4, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

test_that("hexify_polygons works with aperture 7", {
  hex_ids <- c(100, 200, 300)

  result <- hexify_polygons(hex_ids, resolution = 5, aperture = 7, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

# =============================================================================
# HEXIFY_TO_POLYGONS
# =============================================================================

test_that("hexify_to_polygons auto-detects resolution from hex_area", {
  skip_if_not_installed("sf")

  df <- data.frame(
    name = c("A", "B"),
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    hex_diag = c(44.67, 44.67)
  )

  result <- hexify_to_polygons(df, aperture = 3)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)
})

test_that("hexify_to_polygons requires hex_id column", {
  df <- data.frame(name = c("A", "B"), hex_area = c(863, 863))

  expect_error(
    hexify_to_polygons(df),
    "must contain 'hex_id' column"
  )
})

test_that("hexify_to_polygons requires hex_area column", {
  df <- data.frame(hex_id = c(12847, 12532))

  expect_error(
    hexify_to_polygons(df),
    "must contain 'hex_area' column"
  )
})

# =============================================================================
# HEXIFY_PLOT
# =============================================================================

test_that("hexify_plot creates plot without error", {
  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    hex_diag = c(44.67, 44.67)
  )

  expect_silent(hexify_plot(df, col = "lightblue"))
})

test_that("hexify_plot validates input", {
  df1 <- data.frame(name = c("A", "B"))
  df2 <- data.frame(hex_id = c(1, 2))

  expect_error(hexify_plot(df1), "must contain 'hex_id' column")
  expect_error(hexify_plot(df2), "must contain 'hex_area' column")
})

# =============================================================================
# HEXIFY_GRID_RECT
# =============================================================================

test_that("hexify_grid_rect generates grid for rectangular region", {
  skip_if_not_installed("sf")

  grid <- hexify_grid_rect(
    minlon = 0, maxlon = 5,
    minlat = 45, maxlat = 48,
    area = 5000
  )

  expect_s3_class(grid, "sf")
  expect_true(nrow(grid) > 0)

  geom_types <- sf::st_geometry_type(grid)
  expect_true(all(geom_types == "POLYGON"))
})

# =============================================================================
# HEX_CORNERS_TO_SF
# =============================================================================

test_that("hex_corners_to_sf builds valid polygon", {
  skip_if_not_installed("sf")

  lon <- c(0, 1, 1, 0, -1, -1)
  lat <- c(0, 0.5, 1, 1, 0.5, 0)

  poly <- hex_corners_to_sf(lon, lat)

  expect_s3_class(poly, "sf")
  expect_true(sf::st_is_valid(poly))
  expect_equal(nrow(poly), 1L)
  expect_identical(as.character(sf::st_geometry_type(poly)), "POLYGON")
})

test_that("hex_corners_to_sf closes polygon correctly", {
  skip_if_not_installed("sf")

  lon <- c(0, 1, 1, 0, -1, -1)
  lat <- c(0, 0.5, 1, 1, 0.5, 0)

  poly <- hex_corners_to_sf(lon, lat)

  coords <- sf::st_coordinates(poly)
  expect_true(all(c("X", "Y") %in% colnames(coords)))

  # XY must match provided points + closing vertex
  xy_expected <- rbind(cbind(lon, lat), c(lon[1], lat[1]))
  actual_xy <- unname(as.matrix(coords[, c("X", "Y"), drop = FALSE]))
  expected_xy <- unname(as.matrix(xy_expected))

  expect_equal(actual_xy, expected_xy, tolerance = 0)
})
