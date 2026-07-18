# tests/testthat/test-neighbors.R
# Tests for get_neighbors() / k_ring

test_that("get_neighbors returns 6 neighbors for interior ISEA cell (ap4)", {
  g <- hex_grid(area_km2 = 1000)
  cell <- lonlat_to_cell(10, 50, g)
  expect_false(is_pentagon(cell, g))
  nbrs <- get_neighbors(cell, g)
  expect_type(nbrs, "list")
  expect_length(nbrs, 1)
  # Interior (non-pentagon) hex cells always have exactly 6 neighbors
  expect_equal(length(nbrs[[1]]), 6)
  # Self should not be in neighbors
  expect_false(cell %in% nbrs[[1]])
})

test_that("get_neighbors returns 6 neighbors for interior ISEA cell (ap3)", {
  g <- hex_grid(area_km2 = 1000, aperture = 3)
  cell <- lonlat_to_cell(10, 50, g)
  expect_false(is_pentagon(cell, g))
  nbrs <- get_neighbors(cell, g)
  expect_equal(length(nbrs[[1]]), 6)
})

test_that("get_neighbors returns 6 neighbors for interior ISEA cell (ap7)", {
  g <- hex_grid(area_km2 = 1000, aperture = 7)
  cell <- lonlat_to_cell(10, 50, g)
  expect_false(is_pentagon(cell, g))
  nbrs <- get_neighbors(cell, g)
  expect_equal(length(nbrs[[1]]), 6)
})

test_that("get_neighbors works for H3 cells", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  nbrs <- get_neighbors(cell, g)
  expect_type(nbrs, "list")
  expect_length(nbrs, 1)
  # H3 pentagons (12 per resolution) have 5 neighbors; this point isn't one
  expect_equal(length(nbrs[[1]]), 6)
})

test_that("k_ring with k=2 returns more cells than k=1", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  nbrs_k1 <- get_neighbors(cell, g, k = 1)
  nbrs_k2 <- get_neighbors(cell, g, k = 2)
  expect_true(length(nbrs_k2[[1]]) > length(nbrs_k1[[1]]))
})

test_that("include_self adds origin to result", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  nbrs <- get_neighbors(cell, g, include_self = TRUE)
  expect_true(cell %in% nbrs[[1]])
})

test_that("distances=TRUE returns data.frame", {
  g <- hex_grid(resolution = 5, type = "h3")
  cell <- lonlat_to_cell(10, 50, g)
  result <- get_neighbors(cell, g, k = 2, distances = TRUE)
  expect_true(is.data.frame(result[[1]]))
  expect_true("cell_id" %in% names(result[[1]]))
  expect_true("ring_distance" %in% names(result[[1]]))
  expect_true(all(result[[1]]$ring_distance %in% 1:2))
})

test_that("get_neighbors handles vectorized input", {
  g <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(c(10, 20), c(50, 55), g)
  nbrs <- get_neighbors(cells, g)
  expect_length(nbrs, 2)
})
