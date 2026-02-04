# tests/testthat/test-grid-helpers.R
# Tests for grid helper functions
# Note: grid_global, grid_rect, grid_clip tests removed for CRAN speed

test_that("lonlat_to_cell works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
  expect_true(all(cells > 0))
})

test_that("lonlat_to_cell works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Use the HexData object as grid source
  cells <- lonlat_to_cell(lon = 5, lat = 48, grid = result)

  expect_type(cells, "double")
  expect_length(cells, 1)
})

test_that("lonlat_to_cell works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
})

test_that("cell_to_lonlat returns cell centers", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  expect_length(coords$lon_deg, 2)
})

test_that("cell_to_lonlat works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_length(coords$lon_deg, 2)
})

test_that("lonlat_to_cell -> cell_to_lonlat round-trip lands in same cell", {
  grid <- hex_grid(area_km2 = 1000)

  original_lon <- c(0, 10, -5)
  original_lat <- c(45, 50, 48)

  cells1 <- lonlat_to_cell(original_lon, original_lat, grid)
  centers <- cell_to_lonlat(cells1, grid)
  cells2 <- lonlat_to_cell(centers$lon_deg, centers$lat_deg, grid)

  expect_equal(cells1, cells2)
})

test_that("cell_to_sf creates sf polygons", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  polys <- cell_to_sf(cells, grid)

  expect_s3_class(polys, "sf")
  expect_true("cell_id" %in% names(polys))
  expect_equal(nrow(polys), length(unique(cells)))
})

test_that("cell_to_sf works with HexData (no cell_id)", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  polys <- cell_to_sf(grid = result)

  expect_s3_class(polys, "sf")
  expect_equal(nrow(polys), length(unique(result@cell_id)))
})

test_that("cell_to_sf errors without cell_id for HexGridInfo", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(grid = grid), "cell_id required")
})

test_that("cell_to_sf errors on empty cell_id", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(cell_id = numeric(0), grid = grid), "No valid")
  expect_error(cell_to_sf(cell_id = c(NA, NA), grid = grid), "No valid")
})

test_that("extract_grid works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  g <- hexify:::extract_grid(result)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with legacy hexify_grid", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid errors on invalid input", {
  expect_error(hexify:::extract_grid(list(a = 1)), "Cannot extract grid")
  expect_error(hexify:::extract_grid(data.frame()), "Cannot extract grid")
})

test_that("cell_to_sf returns valid geometries for all cells", {
  skip_if_not_installed("sf")

  # Test that polar cells (at icosahedral vertices) have valid geometries
  # Cell 1 is always quad 0 (north pole), and the last cell is quad 11 (south pole)
  grid <- hex_grid(area_km2 = 100000)
  n_cells <- 10 * (as.integer(grid@aperture)^grid@resolution) + 2

  # Test north and south pole cells
  polar_cells <- c(1, n_cells)
  polys <- cell_to_sf(polar_cells, grid)

  # All geometries must be valid
  validity <- sf::st_is_valid(polys)
  expect_true(all(validity))

  # These pentagon cells should have 5 vertices (6 coords with closing)
  vertex_counts <- sapply(sf::st_geometry(polys), function(g) nrow(sf::st_coordinates(g)))
  expect_true(all(vertex_counts == 6))
})
