# tests/testthat/test-compact.R
# Tests for hex_compact() and hex_uncompact()

test_that("H3 compact/uncompact roundtrip", {
  g <- hex_grid(resolution = 3, type = "h3")

  # Get children of a parent cell
  parent <- lonlat_to_cell(10, 50, g)
  parent_res2 <- cpp_h3_cellToParent(parent, 2L)
  children <- cpp_h3_cellToChildren(parent_res2, 3L)[[1]]

  # Compact children should give parent

  compacted <- hex_compact(children, g)
  expect_true(parent_res2 %in% compacted)

  # Uncompact should give back children
  uncompacted <- hex_uncompact(compacted, g, target_resolution = 3L)
  expect_equal(sort(uncompacted), sort(children))
})

test_that("ISEA Z7 compact/uncompact roundtrip", {
  parent_grid <- hex_grid(resolution = 1, aperture = 7)
  child_grid <- hex_grid(resolution = 2, aperture = 7)

  parent <- 8L
  children <- cell_to_index(get_children(parent, parent_grid)[[1]], child_grid)

  compacted <- hex_compact(children, child_grid)
  expect_equal(compacted, cell_to_index(parent, parent_grid))

  uncompacted <- hex_uncompact(compacted, child_grid, target_resolution = 2L)
  expect_equal(sort(uncompacted), sort(children))
})

test_that("hex_compact with no full parents returns unchanged", {
  g <- hex_grid(resolution = 3, type = "h3")
  cells <- lonlat_to_cell(c(10, 20, 30), c(50, 55, 60), g)
  compacted <- hex_compact(cells, g)
  # Non-sibling cells should not be compacted
  expect_equal(length(compacted), length(cells))
})

test_that("compaction works on every pure ISEA aperture", {
  for (aperture in c(3L, 4L, 7L)) {
    label <- sprintf("aperture %d", aperture)
    resolution <- 3L
    grid <- hex_grid(resolution = resolution, aperture = aperture)
    n_cells <- 2 + 10 * aperture^resolution

    all_cells <- cell_to_index(seq_len(n_cells), grid)
    compacted <- hex_compact(all_cells, grid)

    expect_lt(length(compacted), n_cells)
    expect_equal(sort(hex_uncompact(compacted, grid, resolution)),
                 sort(all_cells), info = label)
  }
})

test_that("a full sibling set compacts to its parent, at any aperture", {
  for (aperture in c(3L, 4L, 7L)) {
    parent_grid <- hex_grid(resolution = 2, aperture = aperture)
    child_grid <- hex_grid(resolution = 3, aperture = aperture)
    label <- sprintf("aperture %d", aperture)

    # A cell away from the icosahedron vertices
    parent <- which(!is_pentagon(seq_len(2 + 10 * aperture^2), parent_grid))[10]
    children <- cell_to_index(get_children(parent, parent_grid)[[1]], child_grid)

    expect_equal(hex_compact(children, parent_grid),
                 cell_to_index(parent, parent_grid), info = label)

    # One child short is not a sibling set
    expect_equal(sort(hex_compact(children[-1], parent_grid)),
                 sort(children[-1]), info = label)
  }
})

test_that("uncompaction expands to the cells the hierarchy holds", {
  for (aperture in c(3L, 4L, 7L)) {
    parent_grid <- hex_grid(resolution = 2, aperture = aperture)
    child_grid <- hex_grid(resolution = 3, aperture = aperture)
    n_parent <- 2 + 10 * aperture^2
    label <- sprintf("aperture %d", aperture)

    parents <- seq_len(n_parent)
    expanded <- hex_uncompact(cell_to_index(parents, parent_grid),
                              parent_grid, target_resolution = 3L)

    # Every cell of the child grid, named once
    expect_equal(sort(expanded),
                 sort(cell_to_index(seq_len(2 + 10 * aperture^3), child_grid)),
                 info = label)
  }
})

test_that("an index naming no cell is refused by name", {
  grid <- hex_grid(resolution = 3, aperture = 3)

  # A digit appended to the vertex cell "1100" writes a string that names
  # cell 278 of a grid holding 272
  expect_error(hex_compact(c("11000", "11001"), grid),
               "name no cell of the grid")
  expect_error(hex_uncompact("11001", grid, target_resolution = 4L),
               "name no cell of the grid")
  expect_error(hex_compact("1", grid), "shorter than that")
})
