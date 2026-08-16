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

test_that("ISEA neighbours stay inside the grid's cell ID range", {
  for (ap in c(3, 4, 7)) {
    for (res in 0:4) {
      n_cells <- 10 * ap^res + 2
      nbrs <- get_neighbors(seq_len(n_cells), hex_grid(resolution = res, aperture = ap))
      ids <- unlist(nbrs)

      expect_true(all(ids >= 1 & ids <= n_cells),
                  info = sprintf("aperture %d res %d neighbour ID range", ap, res))
    }
  }
})

test_that("ISEA grids hold 12 pentagons and hexagons everywhere else", {
  # The icosahedron's 12 vertices carry a pentagon each, so exactly 12 cells
  # have five neighbours and the rest have six.
  for (ap in c(3, 4, 7)) {
    for (res in 0:4) {
      n_cells <- 10 * ap^res + 2
      cells <- seq_len(n_cells)
      nbrs <- get_neighbors(cells, hex_grid(resolution = res, aperture = ap))
      label <- sprintf("aperture %d res %d", ap, res)

      expect_equal(sum(lengths(nbrs) == 5), 12L,
                   info = sprintf("%s pentagon count", label))
      expect_equal(sum(lengths(nbrs) == 6), n_cells - 12L,
                   info = sprintf("%s hexagon count", label))
      expect_true(all(vapply(nbrs, function(x) !anyDuplicated(x), TRUE)),
                  info = sprintf("%s duplicate neighbours", label))
      expect_true(all(mapply(function(c, x) !(c %in% x), cells, nbrs)),
                  info = sprintf("%s self in neighbours", label))
    }
  }
})

test_that("ISEA adjacency is symmetric", {
  for (ap in c(3, 4, 7)) {
    for (res in 0:4) {
      n_cells <- 10 * ap^res + 2
      cells <- seq_len(n_cells)
      nbrs <- get_neighbors(cells, hex_grid(resolution = res, aperture = ap))

      from <- rep(cells, lengths(nbrs))
      to <- unlist(nbrs)

      expect_setequal(paste(from, to), paste(to, from))
    }
  }
})

test_that("ISEA neighbours are the adjacent cells, not the ring beyond", {
  # Every neighbour centre sits about one cell spacing away; a coordinate that
  # stepped into the wrong quad would land near twice that.
  for (ap in c(3, 4, 7)) {
    for (res in c(2, 3)) {
      g <- hex_grid(resolution = res, aperture = ap)
      n_cells <- 10 * ap^res + 2
      cells <- seq_len(n_cells)
      nbrs <- get_neighbors(cells, g)

      ll <- cell_to_lonlat(cells, g)
      from <- rep(cells, lengths(nbrs))
      to <- unlist(nbrs)
      p <- pi / 180
      h <- sin((ll$lat[to] - ll$lat[from]) * p / 2)^2 +
        cos(ll$lat[from] * p) * cos(ll$lat[to] * p) *
        sin((ll$lon[to] - ll$lon[from]) * p / 2)^2
      km <- 6371 * 2 * asin(pmin(1, sqrt(h)))

      spacing <- sqrt(2 * (5.100656e8 / n_cells) / sqrt(3))
      expect_lt(max(km / spacing), 1.3)
    }
  }
})
