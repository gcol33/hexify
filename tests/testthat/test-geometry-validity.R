# tests/testthat/test-geometry-validity.R
# Tests for geometry validity at icosahedral vertices (pentagon cells)
#
# The ISEA3H grid has 12 icosahedral vertices where pentagons occur
# instead of hexagons. These tests verify that all geometries are valid
# and renderable.

test_that("grid_global produces all valid geometries", {
  skip_on_cran()
  skip_if_not_installed("sf")

  # Use a coarse grid for faster testing
  grid <- hex_grid(area_km2 = 100000)
  global <- grid_global(grid)

  # All geometries should be valid
  validity <- sf::st_is_valid(global)
  invalid_count <- sum(!validity)


  expect_true(all(validity),
    info = sprintf(
      "%d invalid geometries found out of %d total",
      invalid_count, nrow(global)
    )
  )
})

test_that("cell_to_sf produces valid geometries for polar pentagon cells", {
  skip_on_cran()
  skip_if_not_installed("sf")

  # Test at multiple apertures
  for (aperture in c(3L, 4L)) {
    grid <- hex_grid(area_km2 = 100000, aperture = aperture)

    # Get polar cells (north and south pole)
    n_cells <- 10 * (aperture^grid@resolution) + 2
    polar_cells <- c(1, n_cells)  # Quad 0 and quad 11

    polys <- cell_to_sf(polar_cells, grid)

    validity <- sf::st_is_valid(polys)
    expect_true(all(validity),
      info = sprintf("Aperture %d: polar pentagon cells have invalid geometry", aperture)
    )

    # Both are pentagons: 5 corners, the first repeated last
    rings <- cpp_cell_to_corners(numeric(0), as_cell_id(polar_cells), grid@resolution,
                                 aperture, integer(0))
    expect_true(all(vapply(rings, nrow, integer(1)) == 6),
      info = sprintf("Aperture %d: polar pentagon cells should have 6 coordinates", aperture)
    )
  }
})

test_that("grid_global geometries are non-empty", {
  skip_on_cran()
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 100000)
  global <- grid_global(grid)

  # No empty geometries
  empty <- sf::st_is_empty(global)
  expect_false(any(empty),
    info = sprintf("%d empty geometries found", sum(empty))
  )
})

test_that("exactly the 12 vertex cells are pentagons", {
  skip_on_cran()
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 100000)
  global <- grid_global(grid)

  # The 12 cells at icosahedral vertices are pentagons, every other a hexagon
  rings <- cpp_cell_to_corners(numeric(0), global$cell_id, grid@resolution,
                               as.integer(grid@aperture), integer(0))
  corner_counts <- vapply(rings, nrow, integer(1)) - 1L
  expect_equal(sum(corner_counts == 5L), 12L)
  expect_true(all(corner_counts %in% c(5L, 6L)))
})

test_that("grid_rect produces all valid geometries", {
  skip_on_cran()
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 50000)

  # Test a region that doesn't include icosahedral vertices
  europe <- grid_rect(c(-10, 35, 30, 60), grid)

  validity <- sf::st_is_valid(europe)
  expect_true(all(validity),
    info = sprintf("%d invalid geometries in Europe grid", sum(!validity))
  )
})

test_that("near-polar cells stay non-degenerate across apertures", {
  skip_on_cran()
  skip_if_not_installed("sf")

  # One resolution per aperture fine enough that a whole cell fits inside a
  # couple of degrees of the pole.
  cases <- list(
    list(aperture = 3L, resolution = 5L),
    list(aperture = 4L, resolution = 5L),
    list(aperture = 7L, resolution = 4L)
  )

  for (case in cases) {
    label <- sprintf("aperture %d, resolution %d", case$aperture, case$resolution)
    grid <- hex_grid(resolution = case$resolution, aperture = case$aperture)

    n_total <- 2 + 10 * case$aperture^case$resolution
    centers <- cell_to_lonlat(seq_len(n_total), grid)
    polar <- which(abs(centers$lat_deg) > 80)
    expect_gt(length(polar), 0)

    polys <- cell_to_sf(polar, grid)

    empty <- sf::st_is_empty(polys)
    expect_false(any(empty),
      info = sprintf("%s: %d empty polar geometries", label, sum(empty))
    )
    expect_true(all(sf::st_is_valid(polys)),
      info = sprintf("%s: invalid polar geometries", label)
    )

    unique_vertices <- vapply(sf::st_geometry(polys), function(g) {
      nrow(unique(sf::st_coordinates(g)[, 1:2, drop = FALSE]))
    }, integer(1))
    expect_true(all(unique_vertices >= 5),
      info = sprintf("%s: %d cells with fewer than 5 distinct corners",
                     label, sum(unique_vertices < 5))
    )

    # A pole falls inside a cell or on a cell edge, never on a corner
    corners <- hexify_cell_to_sf(polar, grid = grid, return_sf = FALSE)
    expect_lt(max(abs(corners$lat)), 90)
  }
})

test_that("grid_global covers the poles without warning", {
  skip_on_cran()
  skip_if_not_installed("sf")

  for (aperture in c(3L, 4L, 7L)) {
    label <- sprintf("aperture %d", aperture)
    expect_no_warning({
      global <- grid_global(hex_grid(resolution = 2, aperture = aperture))
    })

    expect_true(all(sf::st_is_valid(global)), info = label)
    expect_false(any(sf::st_is_empty(global)), info = label)

    # The cells holding a pole carry it as a corner
    reach <- vapply(sf::st_geometry(global), function(g) {
      max(abs(sf::st_coordinates(g)[, 2]))
    }, numeric(1))
    expect_equal(max(reach), 90, info = label)
  }
})

test_that("grid_global returns every cell once", {
  skip_on_cran()
  skip_if_not_installed("sf")

  global <- grid_global(hex_grid(resolution = 6, aperture = 4))
  expect_identical(global$cell_id, as_cell_id(seq_len(10 * 4^6 + 2)))
})

test_that("the dateline wrap keeps each cell's area", {
  skip_on_cran()
  skip_if_not_installed("sf")

  cases <- list(
    list(aperture = 3L, resolution = 3L),
    list(aperture = 4L, resolution = 3L),
    list(aperture = 7L, resolution = 2L)
  )

  for (case in cases) {
    label <- sprintf("aperture %d, resolution %d", case$aperture, case$resolution)
    grid <- hex_grid(resolution = case$resolution, aperture = case$aperture)
    cells <- seq_len(2 + 10 * case$aperture^case$resolution)

    whole <- cell_to_sf(cells, grid, wrap_dateline = FALSE)
    split <- cell_to_sf(cells, grid, wrap_dateline = TRUE)

    expect_true(all(sf::st_is_valid(split)), info = label)

    bb <- sf::st_bbox(whole)
    expect_gte(bb[["xmin"]], -180, label = label)
    expect_lte(bb[["xmax"]], 180, label = label)
    expect_no_warning(whole_area <- as.numeric(sf::st_area(whole)))
    expect_equal(as.numeric(sf::st_area(split)), whole_area,
                 tolerance = 1e-6, info = label)
  }
})
