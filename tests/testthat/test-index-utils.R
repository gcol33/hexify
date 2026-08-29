
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
  children <- hexify_get_children(parent, aperture = 3)[[1]]

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
  children <- hexify_get_children(index, aperture = 3)[[1]]

  # Aperture 3 should have 3 children
  expect_equal(length(children), 3)
})

test_that("hexify_get_children returns correct number for aperture 4", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 4)
  children <- hexify_get_children(index, aperture = 4, index_type = "zorder")[[1]]

  # Aperture 4 should have 4 children
  expect_equal(length(children), 4)
})

test_that("hexify_get_children returns correct number for aperture 7", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 7)
  children <- hexify_get_children(index, aperture = 7, index_type = "z7")[[1]]

  # Aperture 7 should have 7 children
  expect_equal(length(children), 7)
})

test_that("hexify_get_children returns longer indices", {
  setup_icosa()

  index <- hexify_lonlat_to_index(10, 50, resolution = 3, aperture = 3)
  children <- hexify_get_children(index, aperture = 3)[[1]]

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

test_that("hexify_index_to_lonlat agrees with cell_to_lonlat for apertures 3 and 4", {
  setup_icosa()

  # Both are cell-centre lookups for the same cell, so they must agree exactly.
  # They used to differ at the polar pentagons, where folding (0, 0) through the
  # quad frame returns the icosahedron vertex rather than the pole.
  set.seed(311)
  lon <- c(runif(30, -170, 170), 0, 0)
  lat <- c(runif(30, -80, 80), 89.9, -89.9)

  for (ap in c(3, 4)) {
    for (res in c(1, 2, 3, 5)) {
      g <- hex_grid(resolution = res, aperture = ap)
      cells <- lonlat_to_cell(lon, lat, g)
      ctr <- cell_to_lonlat(cells, g)
      idx <- cell_to_index(cells, g)

      for (k in seq_along(idx)) {
        got <- hexify_index_to_lonlat(idx[k], aperture = ap)
        expect_equal(as.numeric(got[["lon"]]), ctr$lon[k], tolerance = 1e-4,
                     info = sprintf("ap=%d res=%d cell %d lon", ap, res, k))
        expect_equal(as.numeric(got[["lat"]]), ctr$lat[k], tolerance = 1e-4,
                     info = sprintf("ap=%d res=%d cell %d lat", ap, res, k))
      }
    }
  }
})

# =============================================================================
# VECTORISED INPUT
# =============================================================================

test_that("hexify_lonlat_to_index takes a vector of points", {
  setup_icosa()

  lon <- c(0, 10, 20, -175, 0)
  lat <- c(45, 50, 55, -60, 89.5)

  for (ap in c(3, 4, 7)) {
    label <- sprintf("aperture %d", ap)
    vec <- hexify_lonlat_to_index(lon, lat, resolution = 3, aperture = ap)

    expect_type(vec, "character")
    expect_length(vec, length(lon))

    one_by_one <- vapply(seq_along(lon), function(k) {
      hexify_lonlat_to_index(lon[k], lat[k], resolution = 3, aperture = ap)
    }, character(1))
    expect_identical(vec, one_by_one, info = label)
  }
})

test_that("the index consumers take the vector cell_to_index returns", {
  for (ap in c(3, 4, 7)) {
    label <- sprintf("aperture %d", ap)
    g <- hex_grid(resolution = 3, aperture = ap)
    ix <- cell_to_index(1:5, g)

    cells <- hexify_index_to_cell(ix, aperture = ap)
    expect_s3_class(cells, "data.frame")
    expect_equal(nrow(cells), 5L, info = label)
    expect_named(cells, c("face", "i", "j", "resolution"))

    coords <- hexify_index_to_lonlat(ix, aperture = ap)
    expect_s3_class(coords, "data.frame")
    expect_equal(nrow(coords), 5L, info = label)
    expect_equal(coords$lon, cell_to_lonlat(1:5, g)$lon_deg,
                 tolerance = 1e-6, info = label)

    expect_equal(hexify_get_resolution(ix, aperture = ap), rep(3L, 5),
                 info = label)
    expect_length(hexify_get_parent(ix, aperture = ap), 5L)

    kids <- hexify_get_children(ix, aperture = ap)
    expect_type(kids, "list")
    expect_length(kids, 5L)

    # Every index compares equal to itself and orders against a fixed one
    expect_equal(hexify_compare_indices(ix, ix), rep(0L, 5), info = label)
    expect_length(hexify_compare_indices(ix, ix[1]), 5L)
  }

  g7 <- hex_grid(resolution = 3, aperture = 7)
  ix7 <- cell_to_index(1:5, g7)
  expect_identical(hexify_z7_canonical(ix7), ix7)
})

test_that("hexify_cell_to_index reads face, i and j in step", {
  g <- hex_grid(resolution = 3, aperture = 3)
  qij <- cpp_cell_to_quad_ij(1:5, 3L, 3L)

  expect_identical(
    hexify_cell_to_index(qij$quad, qij$i, qij$j, resolution = 3, aperture = 3),
    cell_to_index(1:5, g)
  )
})

test_that("a missing coordinate or index carries through as missing", {
  setup_icosa()

  idx <- hexify_lonlat_to_index(c(0, NA, 20), c(45, 50, NA),
                                resolution = 3, aperture = 3)
  expect_equal(is.na(idx), c(FALSE, TRUE, TRUE))

  with_na <- c(idx[1], NA)
  expect_equal(is.na(hexify_get_resolution(with_na, aperture = 3)),
               c(FALSE, TRUE))
  expect_equal(is.na(hexify_index_to_lonlat(with_na, aperture = 3)$lon),
               c(FALSE, TRUE))
  expect_equal(is.na(hexify_index_to_cell(with_na, aperture = 3)$face),
               c(FALSE, TRUE))
  expect_equal(is.na(hexify_get_parent(with_na, aperture = 3)),
               c(FALSE, TRUE))
})

test_that("arguments read in step are named when their lengths differ", {
  lonlat <- tryCatch(
    hexify_lonlat_to_index(c(0, 10), 45, resolution = 3, aperture = 3),
    error = conditionMessage
  )
  expect_match(lonlat, "hexify_lonlat_to_index", fixed = TRUE)
  expect_match(lonlat, "`lon` and `lat` are read in step", fixed = TRUE)
  expect_match(lonlat, "lon = 2, lat = 1", fixed = TRUE)

  cell <- tryCatch(
    hexify_cell_to_index(c(0, 1), 0, 0, resolution = 3, aperture = 3),
    error = conditionMessage
  )
  expect_match(cell, "hexify_cell_to_index", fixed = TRUE)
  expect_match(cell, "`face` and `i` and `j` are read in step", fixed = TRUE)
})
