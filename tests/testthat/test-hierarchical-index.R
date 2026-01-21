
# tests/testthat/test-hierarchical-index.R
# Tests for hierarchical index functions

test_that("cell_to_index converts cell IDs to index strings", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  indices <- cell_to_index(cells, grid)

  expect_type(indices, "character")
  expect_length(indices, 2)
  expect_true(all(nchar(indices) > 0))
})

test_that("cell_to_index works with aperture 4", {
  grid <- hex_grid(area_km2 = 1000, aperture = 4)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  indices <- cell_to_index(cells, grid)

  expect_type(indices, "character")
  expect_length(indices, 2)
})

test_that("cell_to_index works with aperture 7", {
  grid <- hex_grid(area_km2 = 10000, aperture = 7)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  indices <- cell_to_index(cells, grid)

  expect_type(indices, "character")
  expect_length(indices, 2)
})

test_that("get_parent returns parent cell for aperture 3", {
  grid <- hex_grid(area_km2 = 1000, aperture = 3)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  parents <- get_parent(cells, grid)

  expect_type(parents, "double")
  expect_length(parents, 2)
  # Parents should be different from children
  expect_false(all(parents == cells))
})

test_that("get_parent returns parent cell for aperture 4", {
  grid <- hex_grid(area_km2 = 1000, aperture = 4)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  parents <- get_parent(cells, grid)

  expect_type(parents, "double")
  expect_length(parents, 2)
})

test_that("get_parent returns parent cell for aperture 7", {
  grid <- hex_grid(area_km2 = 10000, aperture = 7)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  parents <- get_parent(cells, grid)

  expect_type(parents, "double")
  expect_length(parents, 2)
})

test_that("get_children returns child cells for aperture 3", {
  grid <- hex_grid(area_km2 = 100000, aperture = 3)
  cells <- lonlat_to_cell(lon = 0, lat = 45, grid = grid)

  children <- get_children(cells, grid)

  expect_type(children, "list")
  expect_length(children, 1)
  # Aperture 3 should produce 3 children
  expect_length(children[[1]], 3)
})

test_that("get_children returns child cells for aperture 4", {
  grid <- hex_grid(area_km2 = 100000, aperture = 4)
  cells <- lonlat_to_cell(lon = 0, lat = 45, grid = grid)

  children <- get_children(cells, grid)

  expect_type(children, "list")
  expect_length(children, 1)
  # Aperture 4 should produce 4 children
  expect_length(children[[1]], 4)
})

test_that("get_children returns child cells for aperture 7", {
  grid <- hex_grid(area_km2 = 100000, aperture = 7)
  cells <- lonlat_to_cell(lon = 0, lat = 45, grid = grid)

  children <- get_children(cells, grid)

  expect_type(children, "list")
  expect_length(children, 1)
  # Aperture 7 should produce 7 children
  expect_length(children[[1]], 7)
})

test_that("get_parent errors at minimum resolution", {
  grid <- hex_grid(resolution = 0, aperture = 3)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_error(get_parent(cells, grid), "Cannot get parent")
})

test_that("get_children errors at maximum resolution", {
  # Resolution 30 is maximum
  expect_error(
    {
      grid <- hex_grid(resolution = 30, aperture = 3)
      cells <- lonlat_to_cell(lon = 0, lat = 45, grid = grid)
      get_children(cells, grid)
    },
    "Cannot get children"
  )
})
