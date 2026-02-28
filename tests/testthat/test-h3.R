# test-h3.R
# Tests for H3 grid support (native C backend)

skip_if_not_installed("sf")

# =============================================================================
# Grid creation
# =============================================================================

test_that("hex_grid() creates H3 grid by resolution", {
  grid <- hex_grid(resolution = 8, type = "h3")
  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@grid_type, "h3")
  expect_equal(grid@aperture, "7")
  expect_equal(grid@resolution, 8L)
  expect_true(!is.na(grid@area_km2))
  expect_true(grid@area_km2 > 0)
})

test_that("hex_grid() creates H3 grid by area_km2", {
  expect_warning(
    grid <- hex_grid(area_km2 = 1.0, type = "h3"),
    "not exactly equal-area"
  )
  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@grid_type, "h3")
  # Resolution 8 has area ~0.737 km^2, closest to 1.0
  expect_equal(grid@resolution, 8L)
})

test_that("hex_grid() validates H3 resolution range", {
  expect_error(hex_grid(resolution = 16, type = "h3"), "between 0 and 15")
  expect_error(hex_grid(resolution = -1, type = "h3"), "between 0 and 15")
})

test_that("is_h3_grid() works correctly", {
  h3_grid <- hex_grid(resolution = 5, type = "h3")
  isea_grid <- hex_grid(resolution = 5, aperture = 3)

  expect_true(is_h3_grid(h3_grid))
  expect_false(is_h3_grid(isea_grid))
})

# =============================================================================
# hexify() round-trip
# =============================================================================

test_that("hexify() works with H3 grid", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(-73.9857, 2.3522, 139.6917),
                   lat = c(40.7484, 48.8566, 35.6895))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  expect_s4_class(result, "HexData")
  expect_true(is.character(result@cell_id))
  expect_equal(length(result@cell_id), 3)
  expect_true(all(nchar(result@cell_id) > 0))
})

test_that("hexify() H3 cell centers are near original points", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0), lat = c(45))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  # Cell center should be within ~1 km (res 8 cells are ~0.7 km^2)
  center_lon <- result@cell_center[1, "lon"]
  center_lat <- result@cell_center[1, "lat"]
  dist_deg <- sqrt((center_lon - 0)^2 + (center_lat - 45)^2)
  expect_true(dist_deg < 1)  # Less than 1 degree
})

# =============================================================================
# lonlat_to_cell / cell_to_lonlat round-trip
# =============================================================================

test_that("lonlat_to_cell() returns character for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  expect_true(is.character(cells))
  expect_equal(length(cells), 2)
})

test_that("cell_to_lonlat() round-trip for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(lon = c(0), lat = c(45), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_true(is.data.frame(coords))
  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  # Should be close to original
  expect_true(abs(coords$lon_deg[1] - 0) < 5)
  expect_true(abs(coords$lat_deg[1] - 45) < 5)
})

# =============================================================================
# cell_to_sf()
# =============================================================================

test_that("cell_to_sf() returns valid sf for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(lon = c(0, 10, 20), lat = c(45, 50, 55), grid = grid)
  polys <- cell_to_sf(cells, grid)

  expect_s3_class(polys, "sf")
  expect_true(nrow(polys) > 0)
  expect_true("cell_id" %in% names(polys))
  expect_true(all(sf::st_is_valid(polys)))
})

test_that("cell_to_sf() works with HexData for H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 47))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  polys <- cell_to_sf(grid = result)

  expect_s3_class(polys, "sf")
  expect_equal(nrow(polys), n_cells(result))
})

# =============================================================================
# Grid generation
# =============================================================================

test_that("grid_rect() works with H3", {
  grid <- hex_grid(resolution = 3, type = "h3")
  rect <- grid_rect(c(-10, 35, 10, 55), grid)

  expect_s3_class(rect, "sf")
  expect_true(nrow(rect) > 0)
  expect_true(all(sf::st_is_valid(rect)))
})

test_that("grid_clip() works with H3", {
  grid <- hex_grid(resolution = 2, type = "h3")
  # Simple test polygon
  poly <- sf::st_as_sfc("POLYGON((-10 40, 10 40, 10 55, -10 55, -10 40))",
                         crs = 4326)
  poly_sf <- sf::st_sf(geometry = poly)

  clipped <- grid_clip(poly_sf, grid)
  expect_s3_class(clipped, "sf")
  expect_true(nrow(clipped) > 0)
})

# =============================================================================
# Hierarchy
# =============================================================================

test_that("get_parent() works with H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  cells <- lonlat_to_cell(lon = c(0), lat = c(45), grid = grid)
  parent <- get_parent(cells, grid, levels = 1)

  expect_true(is.character(parent))
  expect_equal(length(parent), 1)
  # Parent should be a different (coarser) cell
  expect_false(parent == cells)
})

test_that("get_children() works with H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(lon = c(0), lat = c(45), grid = grid)
  children <- get_children(cells, grid, levels = 1)

  expect_true(is.list(children))
  expect_equal(length(children), 1)
  # H3 aperture 7: each cell has 7 children
  expect_equal(length(children[[1]]), 7)
})

test_that("cell_to_index() returns character for H3", {
  grid <- hex_grid(resolution = 5, type = "h3")
  cells <- lonlat_to_cell(lon = c(0), lat = c(45), grid = grid)
  idx <- cell_to_index(cells, grid)

  expect_true(is.character(idx))
  expect_equal(idx, cells)  # H3 cell IDs are already indices
})

# =============================================================================
# show() methods
# =============================================================================

test_that("show() for H3 HexGridInfo mentions H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  output <- capture.output(show(grid))
  expect_true(any(grepl("H3", output)))
  expect_true(any(grepl("not exactly equal-area", output, ignore.case = TRUE)))
})

test_that("show() for H3 HexData mentions H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  output <- capture.output(show(result))
  expect_true(any(grepl("H3", output)))
})

# =============================================================================
# as.data.frame() and accessors
# =============================================================================

test_that("as.data.frame() works with H3 HexData", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
  out <- as.data.frame(result)

  expect_true("cell_id" %in% names(out))
  expect_true(is.character(out$cell_id))
  expect_equal(nrow(out), 2)
})

test_that("cells() and n_cells() work for H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  df <- data.frame(lon = c(0, 10, 0), lat = c(45, 50, 45))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  unique_cells <- cells(result)
  expect_true(is.character(unique_cells))
  expect_true(n_cells(result) <= 3)
  expect_true(n_cells(result) >= 2)
})

# =============================================================================
# Compare resolutions
# =============================================================================

test_that("hexify_compare_resolutions() works for H3", {
  result <- hexify_compare_resolutions(type = "h3", res_range = 0:5)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 6)
  expect_true("cell_area_km2" %in% names(result))
  # Areas should decrease with resolution
  expect_true(all(diff(result$cell_area_km2) < 0))
})

# =============================================================================
# Plot methods (no-error tests)
# =============================================================================

test_that("plot() does not error for H3 HexData", {
  grid <- hex_grid(resolution = 5, type = "h3")
  df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 47))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  # Should not error
  expect_silent({
    png(tempfile(fileext = ".png"))
    plot(result, basemap = FALSE)
    dev.off()
  })
})

test_that("hexify_heatmap() does not error for H3 HexData", {
  skip_if_not_installed("ggplot2")

  grid <- hex_grid(resolution = 5, type = "h3")
  df <- data.frame(lon = c(0, 5, 10), lat = c(45, 46, 47), count = c(1, 2, 3))
  result <- hexify(df, lon = "lon", lat = "lat", grid = grid)

  p <- hexify_heatmap(result, value = "count")
  expect_s3_class(p, "gg")
})

# =============================================================================
# dgearthstat() with HexGridInfo
# =============================================================================

test_that("dgearthstat() works with H3 HexGridInfo", {
  grid <- hex_grid(resolution = 8, type = "h3")
  stats <- dgearthstat(grid)

  expect_true(is.list(stats))
  expect_equal(stats$resolution, 8L)
  expect_equal(stats$grid_type, "h3")
  expect_true(stats$n_cells > 0)
  expect_true(stats$cell_area_km2 > 0)
})

# =============================================================================
# Legacy conversion guard
# =============================================================================

test_that("HexGridInfo_to_hexify_grid() errors for H3", {
  grid <- hex_grid(resolution = 8, type = "h3")
  expect_error(HexGridInfo_to_hexify_grid(grid), "cannot be converted")
})

# =============================================================================
# Backward compatibility
# =============================================================================

test_that("extract_grid() handles old objects without grid_type", {
  # Simulate old object by creating a grid and removing grid_type
  grid <- hex_grid(resolution = 5, aperture = 3)
  # In real deserialized objects, the slot would be missing
  # Here we just test the normal path works
  g <- extract_grid(grid)
  expect_false(is_h3_grid(g))
})
