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

test_that("adjacent ISEA cells are one hop apart, at every aperture", {
  for (aperture in c(3L, 4L, 7L)) {
    for (resolution in 2:3) {
      label <- sprintf("aperture %d resolution %d", aperture, resolution)
      grid <- hex_grid(resolution = resolution, aperture = aperture)
      n_cells <- 2 + 10 * aperture^resolution

      neighbours <- get_neighbors(seq_len(n_cells), grid)
      from <- rep(seq_len(n_cells), lengths(neighbours))
      to <- unlist(neighbours, use.names = FALSE)

      expect_equal(hex_distance(from, to, grid),
                   rep(1L, length(from)), info = label)
    }
  }
})

test_that("ISEA distances match breadth-first search over the neighbours", {
  bfs_hops <- function(source, resolution, aperture, depth) {
    hops <- c(0L)
    names(hops) <- as.character(source)
    frontier <- source

    for (step in seq_len(depth)) {
      reached <- unique(unlist(cpp_get_neighbors_isea(frontier, resolution,
                                                      aperture)))
      reached <- reached[!(as.character(reached) %in% names(hops))]
      if (length(reached) == 0) break
      hops[as.character(reached)] <- step
      frontier <- reached
    }

    hops
  }

  set.seed(11)
  for (aperture in c(3L, 4L, 7L)) {
    for (resolution in 2:3) {
      label <- sprintf("aperture %d resolution %d", aperture, resolution)
      grid <- hex_grid(resolution = resolution, aperture = aperture)
      n_cells <- 2 + 10 * aperture^resolution

      for (source in sample(n_cells, 3)) {
        truth <- bfs_hops(source, resolution, aperture, 4L)
        targets <- as.numeric(names(truth))
        expect_equal(hex_distance(rep(source, length(targets)), targets, grid),
                     as.integer(truth), info = label)
      }
    }
  }
})

test_that("aperture 3 odd resolutions read the rotated cell lattice", {
  # One substrate point in three is a cell there, so the six neighbours stand
  # sqrt(3) substrate units apart and a raw i/j reading triples the distance
  grid <- hex_grid(resolution = 3, aperture = 3)
  cell <- 229

  neighbours <- get_neighbors(cell, grid)[[1]]
  expect_equal(hex_distance(cell, neighbours, grid),
               rep(1L, length(neighbours)))

  expect_equal(cpp_cell_lattice_generator(3L, 3L), c(2, 1))
  expect_equal(cpp_cell_lattice_generator(3L, 4L), c(1, 0))
  expect_equal(cpp_cell_lattice_generator(4L, 3L), c(1, 0))
  expect_equal(cpp_cell_lattice_generator(7L, 3L), c(1, 0))
})

test_that("a mixed aperture sequence walks like a pure one", {
  for (spelling in c("4/3", "4/7", "7/3", "3/4")) {
    for (resolution in 2:3) {
      label <- sprintf("aperture %s resolution %d", spelling, resolution)
      grid <- hex_grid(resolution = resolution, aperture = spelling)
      n_cells <- n_cells(grid)

      neighbours <- get_neighbors(seq_len(n_cells), grid)

      # Twelve icosahedron vertices with five neighbours, six everywhere else
      expect_equal(sum(lengths(neighbours) == 5L), 12L, info = label)
      expect_equal(sum(lengths(neighbours) == 6L), n_cells - 12L, info = label)

      from <- rep(seq_len(n_cells), lengths(neighbours))
      to <- unlist(neighbours, use.names = FALSE)
      expect_equal(hex_distance(from, to, grid), rep(1L, length(from)),
                   info = label)
    }
  }
})

test_that("mixed-sequence distances match breadth-first search", {
  bfs_hops <- function(source, grid, depth) {
    hops <- c(0L)
    names(hops) <- as.character(source)
    frontier <- source

    for (step in seq_len(depth)) {
      reached <- unique(unlist(grid_neighbors_isea(frontier, grid)))
      reached <- reached[!(as.character(reached) %in% names(hops))]
      if (length(reached) == 0) break
      hops[as.character(reached)] <- step
      frontier <- reached
    }

    hops
  }

  set.seed(5)
  for (spelling in c("4/3", "4/7", "7/3")) {
    label <- paste("aperture", spelling)
    grid <- hex_grid(resolution = 3, aperture = spelling)

    for (source in sample(n_cells(grid), 3)) {
      truth <- bfs_hops(source, grid, 4L)
      targets <- as.numeric(names(truth))
      expect_equal(hex_distance(rep(source, length(targets)), targets, grid),
                   as.integer(truth), info = label)
    }
  }
})

test_that("a missing cell ID gives a missing distance", {
  grid <- hex_grid(resolution = 3, aperture = "4/3")
  expect_equal(hex_distance(c(5, NA, 20), c(5, 10, NA), grid),
               c(0L, NA_integer_, NA_integer_))
})
