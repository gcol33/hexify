# tests/testthat/test-topology.R
# Tests for is_pentagon() and hex_distance()

test_that("is_pentagon detects 12 pentagons for H3", {
  skip_on_cran()

  g <- hex_grid(resolution = 1, type = "h3")
  cells_data <- grid_global(g)
  pent <- is_pentagon(cells_data$cell_id, g)
  expect_equal(sum(pent), 12L)
})

test_that("is_pentagon returns logical vector", {
  g <- hex_grid(resolution = 3, type = "h3")
  cells <- lonlat_to_cell(c(10, 20), c(50, 55), g)
  result <- is_pentagon(cells, g)
  expect_type(result, "logical")
  expect_length(result, 2)
})

test_that("hex_distance returns 0 for same cell (H3)", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  expect_equal(hex_distance(cell, cell, g), 0L)
})

test_that("hex_distance returns 1 for adjacent H3 cells", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  nbrs <- get_neighbors(cell, g)
  dist <- hex_distance(cell, nbrs[[1]][1], g)
  expect_equal(dist, 1L)
})

test_that("hex_distance is symmetric", {
  g <- hex_grid(resolution = 5, type = "h3")
  a <- lonlat_to_cell(10, 50, g)
  b <- lonlat_to_cell(10.5, 50.5, g)
  expect_equal(hex_distance(a, b, g), hex_distance(b, a, g))
})

test_that("hex_distance handles vectorized input", {
  g <- hex_grid(resolution = 5, type = "h3")
  a <- lonlat_to_cell(c(10, 20), c(50, 55), g)
  b <- lonlat_to_cell(c(10.1, 20.1), c(50.1, 55.1), g)
  result <- hex_distance(a, b, g)
  expect_length(result, 2)
  expect_type(result, "integer")
})

test_that("hex_distance works for same-quad ISEA cells", {
  g <- hex_grid(area_km2 = 1000)
  a <- lonlat_to_cell(10, 50, g)
  b <- lonlat_to_cell(10, 50, g)  # Same location = same cell
  expect_equal(hex_distance(a, b, g), 0L)
})
