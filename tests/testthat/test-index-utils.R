# tests/testthat/test-index-utils.R
# Tests for hexify_index.R utility functions
#
# Functions tested:
# - hexify_get_parent()
# - hexify_get_children()
# - hexify_get_resolution()
# - hexify_compare_indices()
# - hexify_is_valid_index_type()
# - hexify_default_index_type()
# - hexify_eff_res_to_area()
# - hexify_area_to_eff_res()
# - hexify_eff_res_to_resolution()
# - hexify_resolution_to_eff_res()

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# GET PARENT INDEX
# =============================================================================

test_that("hexify_get_parent returns shorter index for aperture 3", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 5, aperture = 3)
  parent <- hexify_get_parent(index, aperture = 3)

  expect_true(nchar(parent) < nchar(index))
})

test_that("hexify_get_parent returns shorter index for aperture 4", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 5, aperture = 4)
  parent <- hexify_get_parent(index, aperture = 4, index_type = "zorder")

  expect_true(nchar(parent) < nchar(index))
})

test_that("hexify_get_parent returns shorter index for aperture 7", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 5, aperture = 7)
  parent <- hexify_get_parent(index, aperture = 7, index_type = "z7")

  expect_true(nchar(parent) < nchar(index))
})

test_that("hexify_get_parent of child produces shorter index", {
  setup_icosa()

  # Get a parent index
  parent <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 3)

  # Get children
  children <- hexify_get_children(parent, aperture = 3)

  # Verify each child's parent is shorter than the child
  for (child in children) {
    recovered_parent <- hexify_get_parent(child, aperture = 3)
    expect_true(nchar(recovered_parent) < nchar(child))
  }
})

# =============================================================================
# GET CHILDREN INDICES
# =============================================================================

test_that("hexify_get_children returns correct number for aperture 3", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 3)
  children <- hexify_get_children(index, aperture = 3)

  # Aperture 3 should have 3 children
  expect_equal(length(children), 3)
})

test_that("hexify_get_children returns correct number for aperture 4", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 4)
  children <- hexify_get_children(index, aperture = 4, index_type = "zorder")

  # Aperture 4 should have 4 children
  expect_equal(length(children), 4)
})

test_that("hexify_get_children returns correct number for aperture 7", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 7)
  children <- hexify_get_children(index, aperture = 7, index_type = "z7")

  # Aperture 7 should have 7 children
  expect_equal(length(children), 7)
})

test_that("hexify_get_children returns longer indices", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 3)
  children <- hexify_get_children(index, aperture = 3)

  for (child in children) {
    expect_true(nchar(child) > nchar(index))
  }
})

# =============================================================================
# GET RESOLUTION
# =============================================================================

test_that("hexify_get_resolution returns correct value for aperture 3", {
  setup_icosa()

  for (res in c(1, 3, 5, 7)) {
    index <- hexify_lonlat_to_index(10, 50, resolution = res, aperture = 3)
    recovered_res <- hexify_get_resolution(index, aperture = 3)
    expect_equal(recovered_res, res)
  }
})

test_that("hexify_get_resolution returns correct value for aperture 4", {
  setup_icosa()

  for (res in c(1, 3, 5, 7)) {
    index <- hexify_lonlat_to_index(10, 50, resolution = res, aperture = 4)
    recovered_res <- hexify_get_resolution(
      index, aperture = 4, index_type = "zorder"
    )
    expect_equal(recovered_res, res)
  }
})

test_that("hexify_get_resolution returns correct value for aperture 7", {
  setup_icosa()

  for (res in c(1, 3, 5)) {
    index <- hexify_lonlat_to_index(10, 50, resolution = res, aperture = 7)
    recovered_res <- hexify_get_resolution(
      index, aperture = 7, index_type = "z7"
    )
    expect_equal(recovered_res, res)
  }
})

# =============================================================================
# COMPARE INDICES
# =============================================================================

test_that("hexify_compare_indices returns 0 for equal indices", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 5, aperture = 3)
  result <- hexify_compare_indices(index, index)

  expect_equal(result, 0)
})

test_that("hexify_compare_indices returns -1 for smaller first index", {
  result <- hexify_compare_indices("010000", "020000")
  expect_equal(result, -1)
})

test_that("hexify_compare_indices returns 1 for larger first index", {
  result <- hexify_compare_indices("020000", "010000")
  expect_equal(result, 1)
})

test_that("hexify_compare_indices compares by length when equal prefix", {
  # Shorter index should be "less than" longer
  result <- hexify_compare_indices("010", "0100")
  expect_equal(result, -1)
})

# =============================================================================
# VALID INDEX TYPE
# =============================================================================

test_that("hexify_is_valid_index_type returns TRUE for matching types", {
  expect_true(hexify_is_valid_index_type(3, "z3"))
  expect_true(hexify_is_valid_index_type(4, "zorder"))
  expect_true(hexify_is_valid_index_type(7, "z7"))
})

test_that("hexify_is_valid_index_type returns TRUE for auto", {
  expect_true(hexify_is_valid_index_type(3, "auto"))
  expect_true(hexify_is_valid_index_type(4, "auto"))
  expect_true(hexify_is_valid_index_type(7, "auto"))
})

# =============================================================================
# DEFAULT INDEX TYPE
# =============================================================================

test_that("hexify_default_index_type returns correct types", {
  expect_equal(hexify_default_index_type(3), "z3")
  expect_equal(hexify_default_index_type(4), "zorder")
  expect_equal(hexify_default_index_type(7), "z7")
})

# =============================================================================
# EFFECTIVE RESOLUTION CONVERSIONS
# =============================================================================

test_that("hexify_eff_res_to_area returns positive values", {
  for (eff_res in c(1, 5, 10, 15)) {
    area <- hexify_eff_res_to_area(eff_res)
    expect_true(area > 0)
  }
})

test_that("hexify_eff_res_to_area decreases with resolution", {
  areas <- vapply(1:10, hexify_eff_res_to_area, numeric(1))

  # Each resolution should have smaller area than previous
  expect_true(all(diff(areas) < 0))
})

test_that("hexify_area_to_eff_res is inverse of hexify_eff_res_to_area", {
  for (eff_res in c(1, 5, 10, 15)) {
    area <- hexify_eff_res_to_area(eff_res)
    recovered_res <- hexify_area_to_eff_res(area)
    expect_equal(recovered_res, eff_res, tolerance = 1e-10)
  }
})

test_that("hexify_area_to_eff_res returns expected values", {
  # At eff_res 10, area should be ISEA3H_RES10_AREA_KM2 = 863.8006
  eff_res <- hexify_area_to_eff_res(863.8006)
  expect_equal(eff_res, 10, tolerance = 1e-3)
})

test_that("hexify_eff_res_to_resolution returns integer", {
  for (eff_res in c(1, 5, 10)) {
    res <- hexify_eff_res_to_resolution(eff_res)
    expect_true(is.integer(res))
  }
})

test_that("hexify_eff_res_to_resolution formula is correct", {
  # res_index = 2 * eff_res - 1
  expect_equal(hexify_eff_res_to_resolution(1), 1L)
  expect_equal(hexify_eff_res_to_resolution(5), 9L)
  expect_equal(hexify_eff_res_to_resolution(10), 19L)
})

test_that("hexify_resolution_to_eff_res is inverse", {
  for (res in c(1, 5, 9, 19)) {
    eff_res <- hexify_resolution_to_eff_res(res)
    recovered_res <- hexify_eff_res_to_resolution(eff_res)
    expect_equal(recovered_res, res)
  }
})

test_that("hexify_resolution_to_eff_res formula is correct", {
  # eff_res = (resolution + 1) / 2
  expect_equal(hexify_resolution_to_eff_res(1), 1)
  expect_equal(hexify_resolution_to_eff_res(9), 5)
  expect_equal(hexify_resolution_to_eff_res(19), 10)
})

# =============================================================================
# INDEX STRING OPERATIONS
# =============================================================================

test_that("hexify_cell_to_index and hexify_index_to_cell round-trip", {
  setup_icosa()

  face <- 5
  i <- 10
  j <- 15
  resolution <- 5
  aperture <- 3

  # cell -> index -> cell
  index <- hexify_cell_to_index(face, i, j, resolution, aperture)
  cell <- hexify_index_to_cell(index, aperture)

  # Verify structure is returned
  expect_true("face" %in% names(cell))
  expect_true("i" %in% names(cell))
  expect_true("j" %in% names(cell))
  expect_true("resolution" %in% names(cell))
  expect_equal(cell$resolution, resolution)
})

test_that("hexify_lonlat_to_index and hexify_index_to_lonlat are consistent", {
  setup_icosa()

  lon <- 10
  lat <- 50

  # lon/lat -> index -> lon/lat
  index <- hexify_lonlat_to_index(lon, lat, resolution = 5, aperture = 3)
  coords <- hexify_index_to_lonlat(index, aperture = 3)

  # Reconstructed coordinates should be within one cell diameter
  lon_diff <- abs(coords["lon"] - lon)
  lat_diff <- abs(coords["lat"] - lat)

  expect_true(lon_diff < 5)
  expect_true(lat_diff < 5)
})
