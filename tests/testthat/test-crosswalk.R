# test-crosswalk.R
# Tests for h3_crosswalk()

skip_if_not_installed("sf")

# =============================================================================
# ISEA -> H3
# =============================================================================

test_that("h3_crosswalk() maps ISEA to H3", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
  xwalk <- h3_crosswalk(cells, grid)

  expect_true(is.data.frame(xwalk))
  expect_true(all(c("isea_cell_id", "h3_cell_id", "isea_area_km2",
                     "h3_area_km2", "area_ratio") %in% names(xwalk)))
  expect_equal(nrow(xwalk), length(unique(cells)))
  expect_true(is.character(xwalk$h3_cell_id))
  expect_true(all(xwalk$area_ratio > 0))
})

test_that("h3_crosswalk() auto-selects H3 resolution", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0), c(45), grid)
  xwalk <- h3_crosswalk(cells, grid)

  # H3 area should be in the same ballpark as ISEA area
  expect_true(xwalk$area_ratio[1] > 0.1)
  expect_true(xwalk$area_ratio[1] < 10)
})

test_that("h3_crosswalk() accepts explicit H3 resolution", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0), c(45), grid)
  xwalk <- h3_crosswalk(cells, grid, h3_resolution = 3)

  expect_equal(nrow(xwalk), 1)
  expect_true(is.character(xwalk$h3_cell_id))
})

test_that("h3_crosswalk() deduplicates cell IDs", {
  grid <- hex_grid(area_km2 = 1000)
  # Two points that map to the same cell
  cells <- lonlat_to_cell(c(0, 0.001), c(45, 45.001), grid)
  xwalk <- h3_crosswalk(cells, grid)

  # Result should have unique ISEA cell IDs
  expect_equal(nrow(xwalk), length(unique(cells)))
})

test_that("h3_crosswalk() works with HexData", {
  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  xwalk <- h3_crosswalk(grid = result)

  expect_true(is.data.frame(xwalk))
  expect_true(nrow(xwalk) > 0)
})

# =============================================================================
# H3 -> ISEA
# =============================================================================

test_that("h3_crosswalk() maps H3 to ISEA", {
  h3 <- hex_grid(resolution = 5, type = "h3")
  isea <- hex_grid(area_km2 = 1000)
  h3_cells <- lonlat_to_cell(c(0, 10), c(45, 50), h3)

  xwalk <- h3_crosswalk(h3_cells, h3, isea_grid = isea, direction = "h3_to_isea")

  expect_true(is.data.frame(xwalk))
  expect_equal(nrow(xwalk), length(unique(h3_cells)))
  expect_true(is.numeric(xwalk$isea_cell_id))
  expect_true(is.character(xwalk$h3_cell_id))
})

test_that("h3_crosswalk() errors without isea_grid for h3_to_isea", {
  h3 <- hex_grid(resolution = 5, type = "h3")
  h3_cells <- lonlat_to_cell(c(0), c(45), h3)

  expect_error(
    h3_crosswalk(h3_cells, h3, direction = "h3_to_isea"),
    "isea_grid is required"
  )
})

# =============================================================================
# Validation
# =============================================================================

test_that("h3_crosswalk() errors on wrong grid type for direction", {
  h3 <- hex_grid(resolution = 5, type = "h3")
  isea <- hex_grid(area_km2 = 1000)

  h3_cells <- lonlat_to_cell(c(0), c(45), h3)
  isea_cells <- lonlat_to_cell(c(0), c(45), isea)

  expect_error(
    h3_crosswalk(h3_cells, h3, direction = "isea_to_h3"),
    "requires an ISEA grid"
  )
  expect_error(
    h3_crosswalk(isea_cells, isea, direction = "h3_to_isea"),
    "requires an H3 grid"
  )
})

test_that("h3_crosswalk() validates H3 resolution range", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0), c(45), grid)

  expect_error(
    h3_crosswalk(cells, grid, h3_resolution = 16),
    "between 0 and 15"
  )
})

test_that("h3_crosswalk() errors when cell_id is NULL and grid is HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  expect_error(h3_crosswalk(grid = grid), "cell_id required")
})
