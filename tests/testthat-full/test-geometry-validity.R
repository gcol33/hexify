# tests/testthat-full/test-geometry-validity.R
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

    # Both should have 5 vertices (6 coords with closing)
    vertex_counts <- sapply(sf::st_geometry(polys), function(g) nrow(sf::st_coordinates(g)))
    expect_true(all(vertex_counts == 6),
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

test_that("pentagon cells have approximately 5 unique vertices", {
  skip_on_cran()
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 100000)
  global <- grid_global(grid)

  # Count unique vertices for each cell
  vertex_counts <- sapply(sf::st_geometry(global), function(g) {
    coords <- sf::st_coordinates(g)
    nrow(unique(coords[, 1:2]))
  })

  # Most cells should have 6 vertices (hexagons)
  # 12 cells should have ~5 vertices (pentagons at icosahedral vertices)
  n_pentagons <- sum(vertex_counts <= 5)
  n_hexagons <- sum(vertex_counts >= 6)

  expect_true(n_pentagons >= 1 && n_pentagons <= 20,
    info = sprintf("Expected 1-20 pentagons, found %d", n_pentagons)
  )
  expect_true(n_hexagons > n_pentagons,
    info = "Hexagons should outnumber pentagons"
  )
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
