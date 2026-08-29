
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

test_that("get_parent honours levels on pure ISEA grids", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 4, aperture = ap)
    coarse <- hex_grid(resolution = 3, aperture = ap)

    cells <- lonlat_to_cell(
      lon = c(0, 10, -60, 120, -140),
      lat = c(45, 50, -20, 5, -65),
      grid = g
    )

    expect_equal(
      get_parent(cells, g, levels = 2),
      get_parent(get_parent(cells, g), coarse)
    )
  }
})

test_that("get_children honours levels on pure ISEA grids", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 2, aperture = ap)
    mid <- hex_grid(resolution = 3, aperture = ap)

    cells <- lonlat_to_cell(lon = c(0, 120), lat = c(45, -20), grid = g)

    twice <- lapply(get_children(cells, g), function(kids) {
      sort(unique(unlist(get_children(kids, mid), use.names = FALSE)))
    })
    direct <- lapply(get_children(cells, g, levels = 2), sort)

    expect_equal(direct, twice)
    expect_length(direct[[1]], ap * ap)
  }
})

test_that("a cell is among the children of its two-level ancestor", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 3, aperture = ap)
    coarse <- hex_grid(resolution = 1, aperture = ap)

    cells <- lonlat_to_cell(
      lon = seq(-170, 170, by = 20),
      lat = rep(c(-55, -20, 15, 50), length.out = 18),
      grid = g
    )
    ancestors <- get_parent(cells, g, levels = 2)

    if (ap == 7) {
      # A pentagon has six aperture-7 children, so a Z7 index does not name
      # seven distinct descendants below one.
      keep <- !is_pentagon(ancestors, coarse)
      cells <- cells[keep]
      ancestors <- ancestors[keep]
    }

    expect_gt(length(cells), 0)

    kids <- get_children(ancestors, coarse, levels = 2)
    expect_true(all(mapply(function(cell, k) cell %in% k, cells, kids)))
  }
})

# =============================================================================
# Children at the twelve icosahedron vertices
# =============================================================================

test_that("get_children returns cells that exist in the child grid", {
  for (aperture in c(3L, 4L, 7L)) {
    for (resolution in 1:3) {
      grid <- hex_grid(resolution = resolution, aperture = aperture)
      n_parent <- 2 + 10 * aperture^resolution
      n_child <- 2 + 10 * aperture^(resolution + 1)
      label <- sprintf("aperture %d, resolution %d", aperture, resolution)

      kids <- unlist(get_children(seq_len(n_parent), grid))

      expect_false(anyNA(kids), info = label)
      expect_true(all(kids >= 1 & kids <= n_child), info = label)
    }
  }
})

test_that("get_children inverts get_parent, pentagons included", {
  for (aperture in c(3L, 4L, 7L)) {
    resolution <- 2L
    grid <- hex_grid(resolution = resolution, aperture = aperture)
    child_grid <- hex_grid(resolution = resolution + 1L, aperture = aperture)
    n_parent <- 2 + 10 * aperture^resolution
    n_child <- 2 + 10 * aperture^(resolution + 1)
    label <- sprintf("aperture %d", aperture)

    kids <- get_children(seq_len(n_parent), grid)

    # Every child names its parent back
    for (k in seq_len(n_parent)) {
      expect_equal(unique(get_parent(kids[[k]], child_grid)), k,
                   info = sprintf("%s, parent %d", label, k))
    }

    # And between them the parents claim the whole child grid, once each
    claimed <- unlist(kids)
    expect_equal(sort(claimed), seq_len(n_child), info = label)
  }
})

test_that("a pentagon has one child fewer than the aperture allows", {
  # Aperture 7 refines a hexagon into seven cells; the vertex deficit takes one
  grid <- hex_grid(resolution = 2, aperture = 7)
  n_parent <- 2 + 10 * 7^2
  counts <- lengths(get_children(seq_len(n_parent), grid))
  pentagons <- is_pentagon(seq_len(n_parent), grid)

  expect_equal(sum(pentagons), 12L)
  expect_true(all(counts[pentagons] == 6L))
  expect_true(all(counts[!pentagons] == 7L))
})

test_that("the round trip in both directions holds at a pentagon", {
  parent_grid <- hex_grid(resolution = 2, aperture = 3)
  child_grid <- hex_grid(resolution = 3, aperture = 3)

  kids <- get_children(92L, parent_grid)[[1]]
  expect_true(all(kids <= 272))
  expect_no_error(back <- get_parent(kids, child_grid))
  expect_equal(unique(back), 92)
})
