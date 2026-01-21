
# tests/testthat/test-cell-id.R
# Tests for cell ID functions
#
# Functions tested:
# - hexify_lonlat_to_cell()
# - hexify_cell_to_lonlat()
# - hexify_cell_to_quad_ij()
# - hexify_cell_to_icosa_tri()
# - hexify_cell_id_to_quad_ij()
# - hexify_grid_to_cell()
# - hexify_grid_cell_to_lonlat()

# =============================================================================
# LONLAT TO CELL (CELL ID)
# =============================================================================

test_that("lonlat_to_cell returns positive integers for aperture 3", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)

  cell_ids <- hexify_lonlat_to_cell(lons, lats, resolution = 10, aperture = 3)

  expect_true(all(cell_ids > 0))
  expect_true(is.numeric(cell_ids))
  expect_true(all(is.finite(cell_ids)))
})

test_that("lonlat_to_cell returns positive integers for aperture 4", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)

  cell_ids <- hexify_lonlat_to_cell(lons, lats, resolution = 10, aperture = 4)

  expect_true(all(cell_ids > 0))
  expect_true(is.numeric(cell_ids))
  expect_true(all(is.finite(cell_ids)))
})

test_that("lonlat_to_cell returns positive integers for aperture 7", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)

  cell_ids <- hexify_lonlat_to_cell(lons, lats, resolution = 5, aperture = 7)

  expect_true(all(cell_ids > 0))
  expect_true(is.numeric(cell_ids))
  expect_true(all(is.finite(cell_ids)))
})

test_that("lonlat_to_cell handles extreme coordinates", {
  extreme_coords <- data.frame(
    lon = c(-180, 180, 0, 0),
    lat = c(0, 0, -89, 89)
  )

  for (i in seq_len(nrow(extreme_coords))) {
    cell_id <- hexify_lonlat_to_cell(
      extreme_coords$lon[i],
      extreme_coords$lat[i],
      resolution = 10,
      aperture = 3
    )

    expect_true(is.finite(cell_id))
    expect_true(cell_id > 0)
  }
})

test_that("resolution 0 returns face numbers (1-20)", {
  lons <- c(0, 60, 120, -60, -120, 0)
  lats <- c(45, 45, 45, 45, 45, -45)

  cell_ids <- hexify_lonlat_to_cell(lons, lats, resolution = 0, aperture = 3)

  expect_true(all(cell_ids >= 1 & cell_ids <= 20))
})

# =============================================================================
# CELL TO LONLAT
# =============================================================================

test_that("cell_to_lonlat returns valid coordinates", {
  cell_ids <- c(100, 1000, 10000)

  coords <- hexify_cell_to_lonlat(cell_ids, resolution = 10, aperture = 3)

  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  expect_true(all(is.finite(coords$lon_deg)))
  expect_true(all(is.finite(coords$lat_deg)))
  expect_true(all(coords$lon_deg >= -180 & coords$lon_deg <= 180))
  expect_true(all(coords$lat_deg >= -90 & coords$lat_deg <= 90))
})

# =============================================================================
# ROUND-TRIP CONVERSION
# =============================================================================

test_that("round-trip conversion works for aperture 3", {
  test_lons <- c(0, 5, 10, -30, 45, 90, -90)
  test_lats <- c(0, 45, -30, 15, -60, 30, -15)
  resolution <- 10
  tolerance <- 1.0  # degrees

  cell_ids <- hexify_lonlat_to_cell(test_lons, test_lats, resolution, aperture = 3)
  coords <- hexify_cell_to_lonlat(cell_ids, resolution, aperture = 3)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])

    # Account for longitude wrapping
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

test_that("round-trip conversion works for aperture 4", {
  test_lons <- c(0, 5, 10, -30, 45)
  test_lats <- c(0, 45, -30, 15, -60)
  resolution <- 10
  tolerance <- 1.0

  cell_ids <- hexify_lonlat_to_cell(test_lons, test_lats, resolution, aperture = 4)
  coords <- hexify_cell_to_lonlat(cell_ids, resolution, aperture = 4)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

test_that("round-trip conversion works for aperture 7", {
  test_lons <- c(0, 5, 10, -30, 45)
  test_lats <- c(0, 45, -30, 15, -60)
  resolution <- 5
  # Aperture 7 uses surrogate-substrate conversion with rounding
  tolerance <- 4.0

  cell_ids <- hexify_lonlat_to_cell(test_lons, test_lats, resolution, aperture = 7)
  coords <- hexify_cell_to_lonlat(cell_ids, resolution, aperture = 7)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

# =============================================================================
# CELL TO QUAD IJ
# =============================================================================

test_that("cell_to_quad_ij returns valid structure", {
  cell_ids <- c(100, 1000, 10000)

  result <- hexify_cell_to_quad_ij(cell_ids, resolution = 10, aperture = 3)

  expect_true("quad" %in% names(result))
  expect_true("i" %in% names(result))
  expect_true("j" %in% names(result))
  expect_true(all(result$quad >= 0 & result$quad <= 11))
  expect_true(all(result$i >= 0))
  expect_true(all(result$j >= 0))
})

test_that("cell_to_quad_ij round-trips with quad_ij_to_cell", {
  cell_ids <- c(100, 1000, 5000, 10000)
  resolution <- 10

  for (ap in c(3, 4)) {
    result <- hexify_cell_to_quad_ij(cell_ids, resolution, aperture = ap)
    recovered <- hexify_quad_ij_to_cell(result$quad, result$i, result$j,
                                         resolution, aperture = ap)

    expect_equal(recovered, cell_ids,
                 info = sprintf("aperture %d quad_ij round-trip", ap))
  }
})

# =============================================================================
# CELL TO ICOSA TRIANGLE
# =============================================================================

test_that("cell_to_icosa_tri returns valid structure", {
  cell_ids <- c(100, 1000, 10000)

  result <- hexify_cell_to_icosa_tri(cell_ids, resolution = 10, aperture = 3)

  expect_true("icosa_triangle_face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
  expect_true(all(result$icosa_triangle_face >= 0 & result$icosa_triangle_face <= 19))
  expect_true(all(is.finite(result$icosa_triangle_x)))
  expect_true(all(is.finite(result$icosa_triangle_y)))
})

# =============================================================================
# CELL ID TO QUAD IJ INFO
# =============================================================================

test_that("cell_id_to_quad_ij returns valid cell info", {
  cell_ids <- c(1702, 1954, 100)

  info <- hexify_cell_id_to_quad_ij(cell_ids, resolution = 5, aperture = 3)

  expect_true(is.data.frame(info))
  expect_true("quad" %in% names(info))
  expect_true("i" %in% names(info))
  expect_true("j" %in% names(info))

  # Quads should be 0-11
  expect_true(all(info$quad >= 0 & info$quad <= 11))

  # i, j should be non-negative
  expect_true(all(info$i >= 0))
  expect_true(all(info$j >= 0))
})

# =============================================================================
# GRID-BASED WRAPPERS
# =============================================================================

test_that("grid_to_cell works with hexify_grid object", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))

  expect_true(all(cell_ids > 0))
  expect_length(cell_ids, 2)
})

test_that("grid_cell_to_lonlat works with hexify_grid object", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  cell_ids <- hexify_grid_to_cell(grid, lon = 5, lat = 45)
  coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)

  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
})

# =============================================================================
# BATCH PROCESSING
# =============================================================================

test_that("batch processing handles large datasets", {
  set.seed(123)
  n <- 1000
  lons <- runif(n, -180, 180)
  lats <- runif(n, -90, 90)

  cell_ids <- hexify_lonlat_to_cell(lons, lats, resolution = 10, aperture = 3)

  expect_length(cell_ids, n)
  expect_true(all(cell_ids > 0))

  # Check reasonable uniqueness
  uniqueness <- length(unique(cell_ids)) / n
  expect_gt(uniqueness, 0.8)
})
