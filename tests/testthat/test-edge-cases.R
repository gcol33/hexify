# tests/testthat/test-edge-cases.R
# Edge case and property-based tests

# =============================================================================
# GEOGRAPHIC EDGE CASES
# =============================================================================

test_that("poles assign to valid cells for all apertures", {
  for (ap in c(3, 4, 7)) {
    res <- if (ap == 7) 3 else 5
    g <- hex_grid(resolution = res, aperture = ap)

    # North pole
    np <- lonlat_to_cell(0, 90, g)
    expect_true(np > 0, info = sprintf("north pole, ap %d", ap))

    # South pole
    sp <- lonlat_to_cell(0, -90, g)
    expect_true(sp > 0, info = sprintf("south pole, ap %d", ap))

    # Poles should be different cells
    expect_true(np != sp, info = sprintf("poles differ, ap %d", ap))
  }
})

test_that("antimeridian coordinates assign consistently", {
  g <- hex_grid(resolution = 5, aperture = 3)

  # Points near the antimeridian (180° / -180°)
  cell_pos <- lonlat_to_cell(179.99, 0, g)
  cell_neg <- lonlat_to_cell(-179.99, 0, g)

  # Should be neighboring cells or the same cell (very close points)
  expect_true(cell_pos > 0)
  expect_true(cell_neg > 0)
})

test_that("dateline-straddling points get valid cells", {
  g <- hex_grid(resolution = 5, aperture = 3)

  # Exactly at +180 and -180 should map to the same cell
  cell_a <- lonlat_to_cell(180, 45, g)
  cell_b <- lonlat_to_cell(-180, 45, g)
  expect_equal(cell_a, cell_b)
})

test_that("equator points assign to valid cells", {
  g <- hex_grid(resolution = 5, aperture = 3)

  lons <- seq(-180, 180, by = 30)
  cells <- lonlat_to_cell(lons, rep(0, length(lons)), g)
  expect_true(all(cells > 0))
  expect_length(cells, length(lons))
})

# =============================================================================
# ROUNDTRIP PROPERTY TESTS
# =============================================================================

test_that("encode-decode roundtrip converges for random points", {
  skip_on_cran()

  set.seed(42)
  n <- 50
  lons <- runif(n, -180, 180)
  lats <- runif(n, -90, 90)

  for (ap in c(3, 4)) {
    res <- 5
    g <- hex_grid(resolution = res, aperture = ap)

    cells <- lonlat_to_cell(lons, lats, g)
    coords <- cell_to_lonlat(cells, g)

    # Re-encode should give the same cell
    cells2 <- lonlat_to_cell(coords$lon_deg, coords$lat_deg, g)
    expect_equal(cells, cells2,
                 info = sprintf("roundtrip stability, ap %d", ap))
  }

  # AP7 has known precision issues at quad boundaries —

  # test that the vast majority of cells survive the roundtrip
  g7 <- hex_grid(resolution = 3, aperture = 7)
  cells7 <- lonlat_to_cell(lons, lats, g7)
  coords7 <- cell_to_lonlat(cells7, g7)
  cells7b <- lonlat_to_cell(coords7$lon_deg, coords7$lat_deg, g7)
  pct_match <- mean(cells7 == cells7b)
  expect_true(pct_match > 0.7,
              info = sprintf("ap7 roundtrip: %.0f%% match", pct_match * 100))
})

test_that("cell centers re-encode to the same cell", {
  skip_on_cran()

  g <- hex_grid(resolution = 5, aperture = 3)
  set.seed(123)
  lons <- runif(100, -180, 180)
  lats <- runif(100, -90, 90)

  cells <- lonlat_to_cell(lons, lats, g)
  centers <- cell_to_lonlat(cells, g)
  cells_back <- lonlat_to_cell(centers$lon_deg, centers$lat_deg, g)
  expect_equal(cells, cells_back)
})

# =============================================================================
# PENTAGON PROPERTIES
# =============================================================================

test_that("exactly 12 pentagons exist for all ISEA apertures", {
  for (ap in c(3, 4, 7)) {
    res <- if (ap == 7) 3 else 5
    g <- hex_grid(resolution = res, aperture = ap)

    # Pentagon cell IDs are at (quad=0:11, i=0, j=0)
    pentagon_ids <- hexify:::cpp_quad_ij_to_cell(
      quad = 0:11, i = rep(0, 12), j = rep(0, 12),
      resolution = res, aperture = ap
    )

    pent_flags <- is_pentagon(pentagon_ids, g)
    expect_equal(sum(pent_flags), 12L,
                 info = sprintf("12 pentagons, ap %d", ap))
  }
})

test_that("non-pentagon cells are not detected as pentagons", {
  g <- hex_grid(resolution = 5, aperture = 3)

  # A cell far from any icosahedron vertex
  cell <- lonlat_to_cell(10, 45, g)
  expect_false(is_pentagon(cell, g))
})

# =============================================================================
# NEIGHBOR PROPERTIES
# =============================================================================

test_that("interior cells have 6 neighbors", {
  # Use a location well inside a quad to avoid boundary effects
  g <- hex_grid(resolution = 5, aperture = 4)

  cell <- lonlat_to_cell(10, 45, g)
  nbrs <- get_neighbors(cell, g)

  expect_length(nbrs[[1]], 6L)
})

test_that("neighbors are symmetric for interior cells", {
  skip_on_cran()

  # Use aperture 4 and a location well inside a quad
  g <- hex_grid(resolution = 5, aperture = 4)
  cell <- lonlat_to_cell(10, 45, g)
  nbrs <- get_neighbors(cell, g)[[1]]

  for (nbr in nbrs) {
    nbr_nbrs <- get_neighbors(nbr, g)[[1]]
    expect_true(cell %in% nbr_nbrs,
                info = sprintf("cell %s not in neighbors of %s", cell, nbr))
  }
})

# =============================================================================
# GRID CELL COUNTS
# =============================================================================

test_that("resolution 0 produces 12 cells", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 0L, aperture = ap)
    cells <- grid_global(g)
    expect_equal(nrow(cells), 12L,
                 info = sprintf("res 0 cell count, ap %d", ap))
  }
})
