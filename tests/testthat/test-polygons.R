# tests/testthat/test-polygons.R
# Tests for polygon generation functions
#
# Functions tested:
# - hexify_cell_to_sf()
# - cell_to_sf()
# - hexify_grid_rect()
# - hexify_grid_global()
# - hex_corners_to_sf()

# =============================================================================
# HEXIFY_CELL_TO_SF
# =============================================================================

test_that("hexify_cell_to_sf returns data frame with return_sf=FALSE", {
  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("cell_id", "lon", "lat", "order") %in% names(result)))
})

test_that("hexify_cell_to_sf returns 7 vertices per cell (closed polygon)", {
  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(nrow(result), length(hex_ids) * 7)
  expect_equal(unique(result$order), 1:7)
})

test_that("hexify_cell_to_sf returns valid coordinates", {
  hex_ids <- c(12847, 12532)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_true(all(result$lon >= -180 & result$lon <= 180))
  expect_true(all(result$lat >= -90 & result$lat <= 90))
})

test_that("hexify_cell_to_sf produces closed polygons", {
  hex_ids <- c(12847)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  first_vertex <- result[result$order == 1, c("lon", "lat")]
  last_vertex <- result[result$order == 7, c("lon", "lat")]

  expect_equal(first_vertex$lon, last_vertex$lon)
  expect_equal(first_vertex$lat, last_vertex$lat)
})

test_that("hexify_cell_to_sf returns sf object with return_sf=TRUE", {
  skip_if_not_installed("sf")

  hex_ids <- c(12847, 12532, 22178)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = TRUE)

  expect_s3_class(result, "sf")
  expect_true("cell_id" %in% names(result))
  expect_true("geometry" %in% names(result))
  expect_equal(sf::st_crs(result)$epsg, 4326)

  geom_types <- sf::st_geometry_type(result)
  expect_true(all(geom_types == "POLYGON"))
})

test_that("hexify_cell_to_sf removes duplicates", {
  hex_ids <- c(12847, 12847, 12532)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(length(unique(result$cell_id)), 2)
  expect_equal(nrow(result), 2 * 7)
})

test_that("hexify_cell_to_sf handles NA values", {
  hex_ids <- c(12847, NA, 12532)

  result <- hexify_cell_to_sf(hex_ids, resolution = 10, aperture = 3, return_sf = FALSE)

  expect_equal(length(unique(result$cell_id)), 2)
})

test_that("hexify_cell_to_sf validates aperture", {
  expect_error(
    hexify_cell_to_sf(c(12847), resolution = 10, aperture = 5),
    "aperture must be 3, 4, or 7"
  )
})

test_that("hexify_cell_to_sf validates resolution", {
  expect_error(
    hexify_cell_to_sf(c(12847), resolution = -1, aperture = 3),
    "resolution must be between 0 and 30"
  )

  expect_error(
    hexify_cell_to_sf(c(12847), resolution = 31, aperture = 3),
    "resolution must be between 0 and 30"
  )
})

test_that("hexify_cell_to_sf works with aperture 4", {
  hex_ids <- c(100, 200, 300)

  result <- hexify_cell_to_sf(hex_ids, resolution = 8, aperture = 4, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

test_that("hexify_cell_to_sf works with aperture 7", {
  hex_ids <- c(100, 200, 300)

  result <- hexify_cell_to_sf(hex_ids, resolution = 5, aperture = 7, return_sf = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3 * 7)
})

# =============================================================================
# CELL_TO_SF (grid helper)
# =============================================================================

test_that("cell_to_sf works with HexGrid object", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 1000)
  hex_ids <- c(12847, 12532, 22178)

  result <- cell_to_sf(hex_ids, grid)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
})

test_that("cell_to_sf works with HexData object (extracts unique cells)", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 47))
  hex_data <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  # When called with HexData and no cell_id, it uses the cells from data
  result <- cell_to_sf(grid = hex_data)

  expect_s3_class(result, "sf")
  expect_true(nrow(result) >= 1)
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

test_that("hex_corners_to_sf validates input lengths", {
  expect_error(
    hex_corners_to_sf(c(1, 2, 3), c(1, 2, 3, 4, 5, 6)),
    ""
  )
})

# =============================================================================
# ADDITIONAL INPUT VALIDATION
# =============================================================================

test_that("hexify_cell_to_sf validates cell_id is numeric", {
  expect_error(
    hexify_cell_to_sf("not_numeric", resolution = 10, aperture = 3),
    "cell_id must be numeric"
  )
})

test_that("hexify_cell_to_sf errors on empty input", {
  expect_error(
    hexify_cell_to_sf(numeric(0), resolution = 10, aperture = 3),
    "No valid cell_id values"
  )
})

test_that("hexify_cell_to_sf errors on all NA input", {
  expect_error(
    hexify_cell_to_sf(as.numeric(c()), resolution = 10, aperture = 3),
    "No valid cell_id values"
  )
})

# =============================================================================
# AS_SF METHOD
# =============================================================================

test_that("as_sf creates point geometry from HexData object", {
  skip_if_not_installed("sf")

  df <- data.frame(
    name = c("A", "B"),
    lon = c(10.5, 11.2),
    lat = c(48.5, 49.1)
  )

  hex_data <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  result <- as_sf(hex_data, geometry = "point")

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)
  expect_true(all(sf::st_geometry_type(result) == "POINT"))
})

test_that("as_sf creates polygon geometry from HexData object", {
  skip_if_not_installed("sf")

  df <- data.frame(
    name = c("A", "B"),
    lon = c(10.5, 11.2),
    lat = c(48.5, 49.1)
  )

  hex_data <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  result <- as_sf(hex_data, geometry = "polygon")

  expect_s3_class(result, "sf")
  expect_true(all(sf::st_geometry_type(result) == "POLYGON"))
})

# =============================================================================
# RESOLUTION FROM AREA HELPER
# =============================================================================

test_that(".resolution_from_area returns valid resolution", {
  res <- hexify:::.resolution_from_area(1000, aperture = 3)
  expect_true(is.numeric(res))
  expect_true(res >= 0 && res <= 30)
})

# =============================================================================
# HEXIFY_GRID_RECT
# =============================================================================

test_that("hexify_grid_rect validates parameters", {
  skip_if_not_installed("sf")

  # Test basic functionality was already done, test different apertures
  grid_ap4 <- hexify_grid_rect(
    minlon = 0, maxlon = 5,
    minlat = 45, maxlat = 48,
    area = 5000, aperture = 4
  )

  expect_s3_class(grid_ap4, "sf")
  expect_true(nrow(grid_ap4) > 0)
})

test_that("hexify_grid_rect works with resround", {
  skip_if_not_installed("sf")

  grid_up <- hexify_grid_rect(
    minlon = 0, maxlon = 5,
    minlat = 45, maxlat = 48,
    area = 5000, resround = "up"
  )

  grid_down <- hexify_grid_rect(
    minlon = 0, maxlon = 5,
    minlat = 45, maxlat = 48,
    area = 5000, resround = "down"
  )

  expect_s3_class(grid_up, "sf")
  expect_s3_class(grid_down, "sf")
})

# =============================================================================
# HEXIFY_GRID_GLOBAL
# =============================================================================

test_that("hexify_grid_global works with large area", {
  skip_if_not_installed("sf")

  # Use very large area to avoid warning
  grid <- hexify_grid_global(area = 10000000)

  expect_s3_class(grid, "sf")
  expect_true(nrow(grid) > 0)
})

test_that("hexify_grid_global warns on small area", {
  skip_if_not_installed("sf")

  expect_warning(
    hexify_grid_global(area = 1000),
    "approximately.*cells"
  )
})

# =============================================================================
# EDGE CASES FOR HEXIFY_CELL_TO_SF
# =============================================================================

test_that("hexify_cell_to_sf works with large cell IDs", {
  skip_if_not_installed("sf")

  # Higher resolution -> larger cell IDs
  hex_ids <- c(100000, 200000, 300000)

  result <- hexify_cell_to_sf(
    hex_ids, resolution = 12, aperture = 3, return_sf = TRUE
  )

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 3)
})

test_that("hexify_cell_to_sf handles duplicate IDs correctly", {
  skip_if_not_installed("sf")

  hex_ids <- c(100, 100, 200, 200, 200)

  result <- hexify_cell_to_sf(
    hex_ids, resolution = 5, aperture = 3, return_sf = TRUE
  )

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)  # Should be deduplicated
})

test_that("hexify_cell_to_sf with return_sf=FALSE produces correct structure", {
  hex_ids <- c(100, 200)

  result <- hexify_cell_to_sf(
    hex_ids, resolution = 5, aperture = 3, return_sf = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("cell_id", "lon", "lat", "order") %in% names(result)))
  expect_equal(nrow(result), 2 * 7)  # 7 vertices per hex
})

# =============================================================================
# AS_SF ADDITIONAL TESTS
# =============================================================================

test_that("as_sf preserves all columns", {
  skip_if_not_installed("sf")

  df <- data.frame(
    name = c("A", "B"),
    lon = c(10.5, 11.2),
    lat = c(48.5, 49.1),
    custom_col = c("x", "y")
  )

  hex_data <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  result <- as_sf(hex_data, geometry = "point")

  expect_true("name" %in% names(result))
  expect_true("custom_col" %in% names(result))
})

test_that("as_sf polygon geometry preserves attributes", {
  skip_if_not_installed("sf")

  df <- data.frame(
    lon = c(10.5, 11.2),
    lat = c(48.5, 49.1),
    value = c(100, 200)
  )

  hex_data <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  result <- as_sf(hex_data, geometry = "polygon")

  expect_s3_class(result, "sf")
  expect_true("value" %in% names(result))
})

# =============================================================================
# CELL_TO_SF DEDUPLICATION
# =============================================================================

test_that("cell_to_sf deduplicates cell_ids", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 1000)
  hex_ids <- c(12847, 12847, 12532)  # duplicate

  result <- cell_to_sf(hex_ids, grid)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)  # Should deduplicate
})

# =============================================================================
# HEXIFY_CELL_TO_SF WITH GRID PARAMETER
# =============================================================================

test_that("hexify_cell_to_sf extracts parameters from grid object", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 1000)
  hex_ids <- c(12847, 12532)

  # Pass grid instead of resolution/aperture

  result <- hexify_cell_to_sf(hex_ids, grid = grid)

  expect_s3_class(result, "sf")
  expect_equal(nrow(result), 2)
})

test_that("hexify_cell_to_sf errors when resolution/aperture missing and no grid", {
  expect_error(
    hexify_cell_to_sf(c(12847), resolution = 10),
    "resolution and aperture must be provided"
  )

  expect_error(
    hexify_cell_to_sf(c(12847), aperture = 3),
    "resolution and aperture must be provided"
  )
})
