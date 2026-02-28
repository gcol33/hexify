# test-cell-area.R
# Tests for cell_area() and the H3 virtual column fix

# =============================================================================
# ISEA cell_area()
# =============================================================================

test_that("cell_area() returns constant area for ISEA grids", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0, 10, 20), c(0, 45, 80), grid)
  areas <- cell_area(cells, grid)

  expect_true(is.numeric(areas))
  expect_equal(length(areas), 3)
  # All areas should be the same (equal-area property)
  expect_true(all(areas == areas[1]))
  # Area should match grid spec
  expect_equal(unname(areas[1]), grid@area_km2)
})

test_that("cell_area() returns named vector for ISEA", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
  areas <- cell_area(cells, grid)

  expect_true(!is.null(names(areas)))
  expect_equal(length(names(areas)), 2)
})

test_that("cell_area() works with HexData and NULL cell_id", {
  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  areas <- cell_area(grid = result)

  expect_equal(length(areas), nrow(df))
})

test_that("cell_area() errors when cell_id is NULL and grid is HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  expect_error(cell_area(grid = grid), "cell_id required")
})

# =============================================================================
# H3 cell_area()
# =============================================================================

skip_if_not_installed("sf")

test_that("cell_area() returns per-cell areas for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  # Equator vs high latitude
  cells <- lonlat_to_cell(c(0, 0), c(0, 80), grid)
  areas <- cell_area(cells, grid)

  expect_true(is.numeric(areas))
  expect_equal(length(areas), 2)
  # Areas should differ (H3 is NOT equal-area)
  expect_false(areas[1] == areas[2])
  # Both should be positive
  expect_true(all(areas > 0))
})

test_that("cell_area() caching works for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(c(0), c(45), grid)

  # First call
  a1 <- cell_area(cells, grid)
  # Second call (should hit cache)
  a2 <- cell_area(cells, grid)

  expect_equal(a1, a2)
})

test_that("cell_area() deduplicates H3 cell IDs", {
  grid <- hex_grid(resolution = 5, type = "h3")
  # Same cell repeated
  cells <- lonlat_to_cell(c(0, 0), c(45, 45), grid)
  areas <- cell_area(cells, grid)

  # Should still return 2 values (one per input)
  expect_equal(length(areas), 2)
  expect_equal(unname(areas[1]), unname(areas[2]))
})

test_that("cell_area() with HexData for H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  areas <- cell_area(grid = result)

  expect_equal(length(areas), 2)
  expect_true(all(areas > 0))
})

# =============================================================================
# H3 virtual column fix
# =============================================================================

test_that("$cell_area_km2 returns per-cell areas for H3 HexData", {
  grid <- hex_grid(resolution = 5, type = "h3")
  df <- data.frame(lon = c(0, 0), lat = c(0, 80))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  areas <- result$cell_area_km2
  expect_equal(length(areas), 2)
  # Areas should differ for equator vs polar
  expect_false(areas[1] == areas[2])
})

test_that("[[cell_area_km2]] returns per-cell areas for H3 HexData", {
  grid <- hex_grid(resolution = 5, type = "h3")
  df <- data.frame(lon = c(0, 0), lat = c(0, 80))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  areas <- result[["cell_area_km2"]]
  expect_equal(length(areas), 2)
  expect_false(areas[1] == areas[2])
})

test_that("as.data.frame() returns per-cell areas for H3 HexData", {
  grid <- hex_grid(resolution = 5, type = "h3")
  df <- data.frame(lon = c(0, 0), lat = c(0, 80))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  out <- as.data.frame(result)

  expect_true("cell_area_km2" %in% names(out))
  # Should not be a single repeated value
  expect_false(out$cell_area_km2[1] == out$cell_area_km2[2])
})

test_that("$cell_area_km2 still works for ISEA HexData", {
  grid <- hex_grid(area_km2 = 1000)
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  areas <- result$cell_area_km2
  expect_equal(length(areas), 2)
  # ISEA: all areas should be the same
  expect_equal(areas[1], areas[2])
  expect_equal(areas[1], grid@area_km2)
})
