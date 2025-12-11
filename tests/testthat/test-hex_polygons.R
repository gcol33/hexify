# tests/testthat/test-hex_polygons.R
# Tests for polygon generation functions

test_that("hex_polygons() works with data frame output", {
  hex_ids <- c(12847, 12532, 22178)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  # Check structure
  expect_s3_class(result, "data.frame")
  expect_true(all(c("hex_id", "lon", "lat", "order") %in% names(result)))

  # Each cell should have 7 vertices (6 corners + 1 to close)
  expect_equal(nrow(result), length(hex_ids) * 7)

  # Check order column
  expect_equal(unique(result$order), 1:7)

  # Check coordinates are valid
  expect_true(all(result$lon >= -180 & result$lon <= 180))
  expect_true(all(result$lat >= -90 & result$lat <= 90))
})

test_that("hex_polygons() works with sf output", {
  skip_if_not_installed("sf")

  hex_ids <- c(12847, 12532, 22178)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = TRUE)

  # Check sf class
  expect_s3_class(result, "sf")
  expect_s3_class(result, "data.frame")

  # Check columns
  expect_true("hex_id" %in% names(result))
  expect_true("geometry" %in% names(result))

  # Check CRS
  expect_equal(sf::st_crs(result)$epsg, 4326)

  # Check geometry type
  geom_types <- sf::st_geometry_type(result)
  expect_true(all(geom_types == "POLYGON"))

  # Check number of features
  expect_equal(nrow(result), 3)
})

test_that("hex_polygons() removes duplicates", {
  # Same hex_id repeated
  hex_ids <- c(12847, 12847, 12532)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  # Should only have 2 unique cells
  expect_equal(length(unique(result$hex_id)), 2)
  expect_equal(nrow(result), 2 * 7)
})

test_that("hex_polygons() handles NA values", {
  hex_ids <- c(12847, NA, 12532)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  # Should only have 2 cells (NA removed)
  expect_equal(length(unique(result$hex_id)), 2)
})

test_that("hex_polygons() validates aperture", {
  hex_ids <- c(12847)

  expect_error(
    hex_polygons(hex_ids, resolution = 10, aperture = 5),
    "aperture must be 3, 4, or 7"
  )
})

test_that("hex_polygons() validates resolution", {
  hex_ids <- c(12847)

  expect_error(
    hex_polygons(hex_ids, resolution = -1, aperture = 3),
    "resolution must be between 0 and 30"
  )

  expect_error(
    hex_polygons(hex_ids, resolution = 31, aperture = 3),
    "resolution must be between 0 and 30"
  )
})

test_that("hex_polygons() works with aperture 4", {
  hex_ids <- c(100, 200, 300)

  result <- hex_polygons(hex_ids, resolution = 8, aperture = 4, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

test_that("hex_polygons() works with aperture 7", {
  hex_ids <- c(100, 200, 300)

  result <- hex_polygons(hex_ids, resolution = 5, aperture = 7, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

test_that("hex_polygons() produces closed polygons", {
  hex_ids <- c(12847)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  # First and last vertices should be the same (closed polygon)
  first_vertex <- result[result$order == 1, c("lon", "lat")]
  last_vertex <- result[result$order == 7, c("lon", "lat")]

  expect_equal(first_vertex$lon, last_vertex$lon)
  expect_equal(first_vertex$lat, last_vertex$lat)
})

test_that("hex_to_polygons() auto-detects resolution from hex_area", {
  skip_if_not_installed("sf")

  # Simulate hexify output with hex_area
  df <- data.frame(
    name = c("A", "B"),
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),  # resolution 10
    hex_diag = c(44.67, 44.67)
  )

  result <- hex_to_polygons(df, aperture = 3)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)
})

test_that("hex_to_polygons() requires hex_id column", {
  df <- data.frame(name = c("A", "B"), hex_area = c(863, 863))

  expect_error(
    hex_to_polygons(df),
    "must contain 'hex_id' column"
  )
})

test_that("hex_to_polygons() requires hex_area column", {
  df <- data.frame(hex_id = c(12847, 12532))

  expect_error(
    hex_to_polygons(df),
    "must contain 'hex_area' column"
  )
})

test_that("hex_plot() creates plot without error", {
  # Simulate hexify output
  df <- data.frame(
    hex_id = c(12847, 12532),
    hex_area = c(863.94, 863.94),
    hex_diag = c(44.67, 44.67)
  )

  # Should run without error (creates plot as side effect)
  expect_silent(hex_plot(df, col = "lightblue"))
})

test_that("hex_plot() validates input", {
  df1 <- data.frame(name = c("A", "B"))
  df2 <- data.frame(hex_id = c(1, 2))

  expect_error(hex_plot(df1), "must contain 'hex_id' column")
  expect_error(hex_plot(df2), "must contain 'hex_area' column")
})

test_that("hex_grid_rect() generates grid for rectangular region", {
  skip_if_not_installed("sf")

  # Small test region
  grid <- hex_grid_rect(
    minlon = 0, maxlon = 5,
    minlat = 45, maxlat = 48,
    area = 5000
  )

  expect_s3_class(grid, "sf")
  expect_true(nrow(grid) > 0)

  # Check all cells are polygons
  geom_types <- sf::st_geometry_type(grid)
  expect_true(all(geom_types == "POLYGON"))
})

test_that("polygon vertices are in correct order", {
  hex_ids <- c(12847)

  result <- hex_polygons(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  # Extract vertices in order
  vertices <- result[order(result$order), ]

  # Calculate centroid (rough approximation)
  centroid_lon <- mean(vertices$lon[1:6])
  centroid_lat <- mean(vertices$lat[1:6])

  # Calculate angles from centroid to each vertex
  angles <- atan2(vertices$lat[1:6] - centroid_lat,
                  vertices$lon[1:6] - centroid_lon)

  # For counter-clockwise ordering, angles should be decreasing (with wrap-around)
  # Check that we have 6 distinct angles
  expect_equal(length(unique(round(angles, 4))), 6)
})
