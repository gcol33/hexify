# tests/testthat/test-stats.R
# Tests for grid statistics functions

# =============================================================================
# DGGRID HELPER FUNCTIONS
# =============================================================================

test_that("dggrid_43h_sequence creates correct aperture sequences", {
  # Basic 43H pattern
  seq1 <- dggrid_43h_sequence(2, 3)
  expect_equal(seq1, c(4L, 4L, 3L, 3L, 3L))

  # All aperture 4
  seq2 <- dggrid_43h_sequence(5, 0)
  expect_equal(seq2, c(4L, 4L, 4L, 4L, 4L))

  # All aperture 3
  seq3 <- dggrid_43h_sequence(0, 4)
  expect_equal(seq3, c(3L, 3L, 3L, 3L))

  # Single level each
  seq4 <- dggrid_43h_sequence(1, 1)
  expect_equal(seq4, c(4L, 3L))
})

# =============================================================================
# EARTH STATISTICS
# =============================================================================

test_that("dgearthstat returns correct structure", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  stats <- dgearthstat(grid)

  expect_type(stats, "list")
  expect_named(stats, c(
    "area_km", "n_cells", "cell_area_km2",
    "cell_spacing_km", "cls_km", "resolution", "aperture"
  ))

  # Check values are positive

  expect_gt(stats$area_km, 0)
  expect_gt(stats$n_cells, 0)
  expect_gt(stats$cell_area_km2, 0)
  expect_gt(stats$cell_spacing_km, 0)
  expect_gt(stats$cls_km, 0)
})

test_that("dgearthstat returns correct values for aperture 3", {
  grid <- hexify_grid(area = 10000, aperture = 3)
  stats <- dgearthstat(grid)

  # Earth surface area should be correct
  expect_equal(stats$area_km, EARTH_SURFACE_KM2)

  # Check aperture matches
  expect_equal(stats$aperture, 3)

  # n_cells should equal 10 * 3^resolution for aperture 3
  expected_cells <- 10 * (3 ^ stats$resolution) + 2
  expect_equal(stats$n_cells, expected_cells)

  # Cell area = Earth surface / n_cells
  expected_area <- EARTH_SURFACE_KM2 / expected_cells
  expect_equal(stats$cell_area_km2, expected_area)
})

test_that("dgearthstat works for different apertures", {
  for (ap in c(3, 4, 7)) {
    grid <- hexify_grid(area = 1000, aperture = ap)
    stats <- dgearthstat(grid)

    expect_equal(stats$aperture, ap)
    expect_gt(stats$n_cells, 0)
  }
})

test_that("dgearthstat validates input", {
  expect_error(dgearthstat("not a grid"), "hexify_grid")
  expect_error(dgearthstat(list(x = 1)), "hexify_grid")
})

# =============================================================================
# MAX CELL
# =============================================================================

test_that("max_cell_id returns correct number of cells", {
  grid <- hexify_grid(area = 10000, aperture = 3)
  max_cells <- max_cell_id(grid$resolution, grid$aperture)

  stats <- dgearthstat(grid)
  expect_equal(max_cells, stats$n_cells)
})

test_that("max_cell_id handles resolution 0", {
  # Resolution 0 returns 20 (20 icosahedron faces)
  expect_equal(max_cell_id(0, 3), 20)
  expect_equal(max_cell_id(0, 4), 20)
  expect_equal(max_cell_id(0, 7), 20)
})

# =============================================================================
# CLOSEST RESOLUTION FUNCTIONS
# =============================================================================

test_that("dg_closest_res_to_area finds correct resolution", {
  grid <- list(aperture = 3, topology = "HEXAGON")
  class(grid) <- c("hexify_grid", "dggs", "list")

  # Find resolution for ~1000 km² cells
  res <- dg_closest_res_to_area(grid, area = 1000, metric = TRUE)

  expect_type(res, "double")
  expect_gte(res, 0)
  expect_lte(res, 30)
})

test_that("dg_closest_res_to_area respects rounding", {
  grid <- list(aperture = 3, topology = "HEXAGON")
  class(grid) <- c("hexify_grid", "dggs", "list")

  res_up <- dg_closest_res_to_area(grid, area = 1000, round = "up")
  res_down <- dg_closest_res_to_area(grid, area = 1000, round = "down")
  res_nearest <- dg_closest_res_to_area(grid, area = 1000, round = "nearest")

  expect_gte(res_up, res_down)
})

test_that("dg_closest_res_to_area supports show_info", {
  grid <- list(aperture = 3, topology = "HEXAGON")
  class(grid) <- c("hexify_grid", "dggs", "list")

  expect_message(
    dg_closest_res_to_area(grid, area = 1000, show_info = TRUE),
    "Resolution"
  )
})


# =============================================================================
# RESOLUTION COMPARISON
# =============================================================================

test_that("hexify_compare_resolutions returns correct structure", {
  comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:5)

  expect_s3_class(comparison, "data.frame")
  expect_named(comparison, c(
    "resolution", "n_cells", "cell_area_km2",
    "cell_spacing_km", "cls_km"
  ))
  expect_equal(nrow(comparison), 6)
})

test_that("hexify_compare_resolutions values are monotonic", {
  comparison <- hexify_compare_resolutions(aperture = 3, res_range = 0:10)

  # Number of cells should increase with resolution
  expect_true(all(diff(comparison$n_cells) > 0))

  # Cell area should decrease with resolution
  expect_true(all(diff(comparison$cell_area_km2) < 0))
})


# =============================================================================
# ADDITIONAL COVERAGE FOR HEXIFY_STATS
# =============================================================================

test_that("dgearthstat handles aperture 7 differently", {
  grid <- hexify_grid(area = 10000, aperture = 7)
  stats <- dgearthstat(grid)

  # Aperture 7 uses 12 base faces
  expected_cells <- 10 * (7 ^ stats$resolution) + 2
  expect_equal(stats$n_cells, expected_cells)
})

test_that("dg_closest_res_to_area handles non-metric input", {
  grid <- list(aperture = 3, topology = "HEXAGON")
  class(grid) <- c("hexify_grid", "dggs", "list")

  # Square miles input
  res <- dg_closest_res_to_area(grid, area = 386, metric = FALSE)

  expect_type(res, "double")
  expect_gte(res, 0)
})


test_that("hexify_compare_resolutions works for different apertures", {
  comp_ap3 <- hexify_compare_resolutions(aperture = 3, res_range = 0:3)
  comp_ap4 <- hexify_compare_resolutions(aperture = 4, res_range = 0:3)
  comp_ap7 <- hexify_compare_resolutions(aperture = 7, res_range = 0:3)

  expect_equal(nrow(comp_ap3), 4)
  expect_equal(nrow(comp_ap4), 4)
  expect_equal(nrow(comp_ap7), 4)
})

test_that("hexify_compare_resolutions formats large cell counts", {
  # Resolution 11+ should have cells in millions
  expect_output(
    hexify_compare_resolutions(aperture = 3, res_range = 11:12, print = TRUE),
    "M"  # Million suffix
  )
})

test_that("hexify_compare_resolutions formats thousand cell counts", {
  # Resolution 3-5 should have cells in thousands
  expect_output(
    hexify_compare_resolutions(aperture = 3, res_range = 3:5, print = TRUE),
    "K"  # Thousand suffix
  )
})
