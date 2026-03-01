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
  g <- hex_grid(resolution = 3, aperture = 7)

  # Create 7 Z7 children of a parent
  parent <- "010"  # quad 01, digit 0 (res 1)
  children <- paste0(parent, 0:6)

  compacted <- hex_compact(children, g)
  expect_true(parent %in% compacted)

  uncompacted <- hex_uncompact(compacted, g, target_resolution = 3L)
  # Should have 7^1 = 7 cells at target res relative to the single res-2 parent
  expect_true(all(nchar(uncompacted) == 5))  # "01" + 3 digits
})

test_that("hex_compact with no full parents returns unchanged", {
  g <- hex_grid(resolution = 3, type = "h3")
  cells <- lonlat_to_cell(c(10, 20, 30), c(50, 55, 60), g)
  compacted <- hex_compact(cells, g)
  # Non-sibling cells should not be compacted
  expect_equal(length(compacted), length(cells))
})

test_that("hex_compact errors for non-aperture-7 ISEA", {
  g <- hex_grid(resolution = 3, aperture = 3)
  expect_error(hex_compact(c("0100", "0101", "0102"), g), "aperture 7")
})
