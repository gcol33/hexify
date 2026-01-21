
# tests/testthat/test-conversions.R
# Tests for coordinate system conversion functions
#
# Functions tested:
# - hexify_lonlat_to_quad_ij()
# - hexify_quad_ij_to_cell()
# - hexify_quad_ij_to_xy()
# - hexify_icosa_tri_to_quad_xy()
# - hexify_icosa_tri_to_quad_ij()
# - hexify_quad_xy_to_icosa_tri()
# - hexify_quad_ij_to_icosa_tri()

# =============================================================================
# LON/LAT TO QUAD IJ
# =============================================================================

test_that("lonlat_to_quad_ij returns valid structure", {
  result <- hexify_lonlat_to_quad_ij(lon = 2.35, lat = 48.86,
                                      resolution = 10, aperture = 3)

  expect_true("quad" %in% names(result))
  expect_true("i" %in% names(result))
  expect_true("j" %in% names(result))
  expect_true("icosa_triangle_face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
})

test_that("lonlat_to_quad_ij quad is in valid range", {
  result <- hexify_lonlat_to_quad_ij(lon = 0, lat = 0,
                                      resolution = 10, aperture = 3)

  expect_true(result$quad >= 0 && result$quad <= 11)
})

test_that("lonlat_to_quad_ij i,j are non-negative", {
  result <- hexify_lonlat_to_quad_ij(lon = 16.37, lat = 48.21,
                                      resolution = 10, aperture = 3)

  expect_true(result$i >= 0)
  expect_true(result$j >= 0)
})

test_that("lonlat_to_quad_ij works for all apertures", {
  for (ap in c(3, 4, 7)) {
    res <- if (ap == 7) 5 else 10

    result <- hexify_lonlat_to_quad_ij(lon = 0, lat = 45,
                                        resolution = res, aperture = ap)

    expect_true(result$quad >= 0 && result$quad <= 11,
                info = sprintf("aperture %d", ap))
  }
})

# =============================================================================
# QUAD IJ TO CELL
# =============================================================================

test_that("quad_ij_to_cell returns positive cell IDs", {
  cell_id <- hexify_quad_ij_to_cell(quad = 1, i = 100, j = 50,
                                     resolution = 10, aperture = 3)

  expect_true(cell_id > 0)
  expect_true(is.numeric(cell_id))
})

test_that("quad_ij_to_cell is consistent with lonlat_to_cell", {
  lon <- 16.37
  lat <- 48.21
  resolution <- 10

  for (ap in c(3, 4)) {
    # Get cell via lon/lat
    cell_direct <- hexify_lonlat_to_cell(lon, lat, resolution, aperture = ap)

    # Get cell via quad IJ
    quad_ij <- hexify_lonlat_to_quad_ij(lon, lat, resolution, aperture = ap)
    cell_indirect <- hexify_quad_ij_to_cell(quad_ij$quad, quad_ij$i, quad_ij$j,
                                             resolution, aperture = ap)

    expect_equal(cell_direct, cell_indirect,
                 info = sprintf("aperture %d", ap))
  }
})

# =============================================================================
# QUAD IJ TO XY
# =============================================================================

test_that("quad_ij_to_xy returns valid structure", {
  result <- hexify_quad_ij_to_xy(quad = 1, i = 100, j = 50,
                                  resolution = 10, aperture = 3)

  expect_true("quad_x" %in% names(result))
  expect_true("quad_y" %in% names(result))
  expect_true(is.numeric(result$quad_x))
  expect_true(is.numeric(result$quad_y))
})

# =============================================================================
# ICOSA TRI TO QUAD XY
# =============================================================================

test_that("icosa_tri_to_quad_xy returns valid structure", {
  hexify_build_icosa()

  fwd <- hexify_forward(lon = 2.35, lat = 48.86)

  result <- hexify_icosa_tri_to_quad_xy(
    icosa_triangle_face = fwd["face"],
    icosa_triangle_x = fwd["icosa_triangle_x"],
    icosa_triangle_y = fwd["icosa_triangle_y"]
  )

  expect_true("quad" %in% names(result))
  expect_true("quad_x" %in% names(result))
  expect_true("quad_y" %in% names(result))
  expect_true(result$quad >= 0 && result$quad <= 11)
})

# =============================================================================
# ICOSA TRI TO QUAD IJ
# =============================================================================

test_that("icosa_tri_to_quad_ij returns valid structure", {
  hexify_build_icosa()

  fwd <- hexify_forward(lon = 2.35, lat = 48.86)

  result <- hexify_icosa_tri_to_quad_ij(
    icosa_triangle_face = fwd["face"],
    icosa_triangle_x = fwd["icosa_triangle_x"],
    icosa_triangle_y = fwd["icosa_triangle_y"],
    resolution = 10,
    aperture = 3
  )

  expect_true("quad" %in% names(result))
  expect_true("i" %in% names(result))
  expect_true("j" %in% names(result))
  expect_true(result$quad >= 0 && result$quad <= 11)
})

# =============================================================================
# QUAD XY TO ICOSA TRI (INVERSE)
# =============================================================================

test_that("quad_xy_to_icosa_tri returns valid structure", {
  result <- hexify_quad_xy_to_icosa_tri(quad = 1, quad_x = 0.5, quad_y = 0.3)

  expect_true("icosa_triangle_face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
  expect_true(result$icosa_triangle_face >= 0 && result$icosa_triangle_face <= 19)
})

test_that("icosa_tri <-> quad_xy round-trip is consistent", {
  hexify_build_icosa()

  # Forward: lon/lat -> icosa tri -> quad xy
  fwd <- hexify_forward(lon = 10, lat = 45)

  quad_xy <- hexify_icosa_tri_to_quad_xy(
    icosa_triangle_face = fwd["face"],
    icosa_triangle_x = fwd["icosa_triangle_x"],
    icosa_triangle_y = fwd["icosa_triangle_y"]
  )

  # Inverse: quad xy -> icosa tri
  back <- hexify_quad_xy_to_icosa_tri(
    quad = quad_xy$quad,
    quad_x = quad_xy$quad_x,
    quad_y = quad_xy$quad_y
  )

  # Should recover original icosa tri coordinates
  expect_equal(as.numeric(fwd["face"]), back$icosa_triangle_face)
  expect_equal(as.numeric(fwd["icosa_triangle_x"]), back$icosa_triangle_x,
               tolerance = 1e-10)
  expect_equal(as.numeric(fwd["icosa_triangle_y"]), back$icosa_triangle_y,
               tolerance = 1e-10)
})

# =============================================================================
# QUAD IJ TO ICOSA TRI
# =============================================================================

test_that("quad_ij_to_icosa_tri returns valid structure", {
  result <- hexify_quad_ij_to_icosa_tri(quad = 1, i = 100, j = 50,
                                         resolution = 10, aperture = 3)

  expect_true("icosa_triangle_face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
  expect_true(result$icosa_triangle_face >= 0 && result$icosa_triangle_face <= 19)
})

# =============================================================================
# PIPELINE CONSISTENCY
# =============================================================================

test_that("coordinate pipeline is consistent end-to-end", {
  hexify_build_icosa()

  lon <- 16.37
  lat <- 48.21
  resolution <- 10
  aperture <- 3

  # Full pipeline: lon/lat -> quad IJ -> cell ID -> quad IJ (back)
  quad_ij <- hexify_lonlat_to_quad_ij(lon, lat, resolution, aperture)
  cell_id <- hexify_quad_ij_to_cell(quad_ij$quad, quad_ij$i, quad_ij$j,
                                     resolution, aperture)

  quad_ij_back <- hexify_cell_to_quad_ij(cell_id, resolution, aperture)

  expect_equal(quad_ij$quad, quad_ij_back$quad)
  expect_equal(quad_ij$i, quad_ij_back$i)
  expect_equal(quad_ij$j, quad_ij_back$j)
})

test_that("all apertures have consistent coordinate pipeline", {
  skip_on_cran()  # Detailed loop test
  hexify_build_icosa()

  lon <- 0
  lat <- 45

  for (ap in c(3, 4, 7)) {
    res <- if (ap == 7) 5 else 10

    # lon/lat -> cell -> lonlat
    cell_id <- hexify_lonlat_to_cell(lon, lat, res, aperture = ap)
    coords <- hexify_cell_to_lonlat(cell_id, res, aperture = ap)

    # Should be reasonably close
    lon_error <- abs(coords$lon_deg - lon)
    if (lon_error > 180) lon_error <- 360 - lon_error
    lat_error <- abs(coords$lat_deg - lat)

    max_error <- if (ap == 7) 5.0 else 2.0
    expect_true(lon_error < max_error,
                info = sprintf("aperture %d lon error", ap))
    expect_true(lat_error < max_error,
                info = sprintf("aperture %d lat error", ap))
  }
})

# =============================================================================
# INPUT VALIDATION
# =============================================================================

test_that("hexify_lonlat_to_quad_ij validates aperture", {
  expect_error(
    hexify_lonlat_to_quad_ij(lon = 0, lat = 0, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_lonlat_to_quad_ij validates resolution", {
  expect_error(
    hexify_lonlat_to_quad_ij(lon = 0, lat = 0, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
  expect_error(
    hexify_lonlat_to_quad_ij(lon = 0, lat = 0, resolution = 31, aperture = 3),
    "between 0 and 30"
  )
})

test_that("hexify_quad_ij_to_cell validates aperture", {
  expect_error(
    hexify_quad_ij_to_cell(quad = 0, i = 1, j = 1, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_quad_ij_to_cell validates resolution", {
  expect_error(
    hexify_quad_ij_to_cell(quad = 0, i = 1, j = 1, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
})

test_that("hexify_quad_ij_to_xy validates aperture", {
  expect_error(
    hexify_quad_ij_to_xy(quad = 1, i = 10, j = 5, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_icosa_tri_to_quad_ij validates aperture", {
  expect_error(
    hexify_icosa_tri_to_quad_ij(
      icosa_triangle_face = 0,
      icosa_triangle_x = 0.5,
      icosa_triangle_y = 0.5,
      resolution = 5,
      aperture = 5
    ),
    "3, 4, or 7"
  )
})

# =============================================================================
# CELL TO QUAD IJ
# =============================================================================

test_that("hexify_cell_to_quad_ij returns correct structure", {
  hexify_build_icosa()

  result <- hexify_cell_to_quad_ij(cell_id = 1000, resolution = 5, aperture = 3)

  expect_s3_class(result, "data.frame")
  expect_true("quad" %in% names(result))
  expect_true("i" %in% names(result))
  expect_true("j" %in% names(result))
})

test_that("hexify_cell_to_quad_ij validates aperture", {
  expect_error(
    hexify_cell_to_quad_ij(cell_id = 100, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_cell_to_quad_ij validates resolution", {
  expect_error(
    hexify_cell_to_quad_ij(cell_id = 100, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
})

# =============================================================================
# CELL TO ICOSA TRIANGLE
# =============================================================================

test_that("hexify_cell_to_icosa_tri returns correct structure", {
  hexify_build_icosa()

  result <- hexify_cell_to_icosa_tri(cell_id = 1000, resolution = 5, aperture = 3)

  expect_s3_class(result, "data.frame")
  expect_true("icosa_triangle_face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
})

test_that("hexify_cell_to_icosa_tri validates aperture", {
  expect_error(
    hexify_cell_to_icosa_tri(cell_id = 100, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_cell_to_icosa_tri validates resolution", {
  expect_error(
    hexify_cell_to_icosa_tri(cell_id = 100, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
})

test_that("hexify_quad_ij_to_icosa_tri validates aperture", {
  expect_error(
    hexify_quad_ij_to_icosa_tri(
      quad = 1, i = 10, j = 5, resolution = 5, aperture = 5
    ),
    "3, 4, or 7"
  )
})

test_that("hexify_quad_ij_to_icosa_tri validates resolution", {
  expect_error(
    hexify_quad_ij_to_icosa_tri(
      quad = 1, i = 10, j = 5, resolution = -1, aperture = 3
    ),
    "between 0 and 30"
  )
})

# =============================================================================
# CELL TO QUAD XY (CONTINUOUS COORDINATES)
# =============================================================================

test_that("hexify_cell_to_quad_xy returns correct structure", {
  hexify_build_icosa()

  result <- hexify_cell_to_quad_xy(cell_id = 1000, resolution = 5, aperture = 3)

  expect_s3_class(result, "data.frame")
  expect_true("quad" %in% names(result))
  expect_true("quad_x" %in% names(result))
  expect_true("quad_y" %in% names(result))
})

test_that("hexify_cell_to_quad_xy validates aperture", {
  expect_error(
    hexify_cell_to_quad_xy(cell_id = 100, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_cell_to_quad_xy validates resolution", {
  expect_error(
    hexify_cell_to_quad_xy(cell_id = 100, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
})

test_that("hexify_quad_xy_to_cell returns valid cell ID", {
  hexify_build_icosa()

  cell_id <- hexify_quad_xy_to_cell(
    quad = 1, quad_x = 0.5, quad_y = 0.3,
    resolution = 5, aperture = 3
  )

  expect_true(is.numeric(cell_id))
  expect_true(cell_id > 0)
})

test_that("hexify_quad_xy_to_cell validates aperture", {
  expect_error(
    hexify_quad_xy_to_cell(
      quad = 1, quad_x = 0.5, quad_y = 0.3, resolution = 5, aperture = 5
    ),
    "3, 4, or 7"
  )
})

test_that("hexify_quad_xy_to_cell validates resolution", {
  expect_error(
    hexify_quad_xy_to_cell(
      quad = 1, quad_x = 0.5, quad_y = 0.3, resolution = -1, aperture = 3
    ),
    "between 0 and 30"
  )
})

test_that("quad_xy round-trip is consistent", {
  hexify_build_icosa()

  cell_id <- 1000

  # cell -> quad_xy -> cell
  result1 <- hexify_cell_to_quad_xy(cell_id = cell_id, resolution = 5, aperture = 3)
  cell_id2 <- hexify_quad_xy_to_cell(
    quad = result1$quad, quad_x = result1$quad_x, quad_y = result1$quad_y,
    resolution = 5, aperture = 3
  )

  expect_equal(cell_id, cell_id2)
})

# =============================================================================
# PLANE COORDINATE CONVERSIONS
# =============================================================================

test_that("hexify_icosa_tri_to_plane returns correct structure", {
  hexify_build_icosa()

  fwd <- hexify_forward(lon = 10, lat = 50)

  result <- hexify_icosa_tri_to_plane(
    icosa_triangle_face = fwd["face"],
    icosa_triangle_x = fwd["icosa_triangle_x"],
    icosa_triangle_y = fwd["icosa_triangle_y"]
  )

  expect_s3_class(result, "data.frame")
  expect_true("plane_x" %in% names(result))
  expect_true("plane_y" %in% names(result))
})

test_that("hexify_cell_to_plane returns correct structure", {
  hexify_build_icosa()

  result <- hexify_cell_to_plane(cell_id = 1000, resolution = 5, aperture = 3)

  expect_s3_class(result, "data.frame")
  expect_true("plane_x" %in% names(result))
  expect_true("plane_y" %in% names(result))
})

test_that("hexify_cell_to_plane validates aperture", {
  expect_error(
    hexify_cell_to_plane(cell_id = 100, resolution = 5, aperture = 5),
    "3, 4, or 7"
  )
})

test_that("hexify_cell_to_plane validates resolution", {
  expect_error(
    hexify_cell_to_plane(cell_id = 100, resolution = -1, aperture = 3),
    "between 0 and 30"
  )
})

test_that("hexify_lonlat_to_plane returns correct structure", {
  hexify_build_icosa()

  result <- hexify_lonlat_to_plane(lon = 10, lat = 50)

  expect_s3_class(result, "data.frame")
  expect_true("plane_x" %in% names(result))
  expect_true("plane_y" %in% names(result))
})

test_that("hexify_lonlat_to_plane works for multiple points", {
  hexify_build_icosa()

  lons <- c(0, 45, -120, 180, -180)
  lats <- c(0, 30, -45, 60, -60)

  result <- hexify_lonlat_to_plane(lon = lons, lat = lats)

  expect_equal(nrow(result), 5)
  expect_true(all(is.finite(result$plane_x)))
  expect_true(all(is.finite(result$plane_y)))
})

test_that("hexify_cell_to_plane handles multiple cells", {
  hexify_build_icosa()

  result <- hexify_cell_to_plane(
    cell_id = c(100, 200, 300), resolution = 5, aperture = 3
  )

  expect_equal(nrow(result), 3)
  expect_true(all(is.finite(result$plane_x)))
  expect_true(all(is.finite(result$plane_y)))
})

# =============================================================================
# CROSS-APERTURE PLANE TESTS
# =============================================================================

test_that("plane conversions work for all apertures", {
  skip_on_cran()  # Loop test across apertures
  hexify_build_icosa()

  for (ap in c(3, 4, 7)) {
    # cell -> plane
    plane <- hexify_cell_to_plane(cell_id = 100, resolution = 5, aperture = ap)
    expect_true(is.finite(plane$plane_x))
    expect_true(is.finite(plane$plane_y))
  }
})

# =============================================================================
# HIERARCHICAL INDEX FUNCTIONS
# =============================================================================

test_that("hexify_lonlat_to_h_index returns correct structure", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_lonlat_to_h_index(grid, lon = 0, lat = 45)

  expect_s3_class(result, "data.frame")
  expect_true("h_index" %in% names(result))
  expect_true("face" %in% names(result))
  expect_type(result$h_index, "character")
  expect_type(result$face, "integer")
})

test_that("hexify_lonlat_to_h_index handles multiple points", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  lons <- c(0, 10, -5)
  lats <- c(45, 30, -20)

  result <- hexify_lonlat_to_h_index(grid, lon = lons, lat = lats)

  expect_equal(nrow(result), 3)
  expect_true(all(nchar(result$h_index) > 0))
})

test_that("hexify_lonlat_to_h_index handles NA values", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_lonlat_to_h_index(grid, lon = c(0, NA), lat = c(45, NA))

  expect_equal(nrow(result), 2)
  expect_true(is.na(result$h_index[2]))
  expect_true(is.na(result$face[2]))
})

test_that("hexify_lonlat_to_h_index validates grid object", {
  expect_error(
    hexify_lonlat_to_h_index(list(x = 1), lon = 0, lat = 0),
    "hexify_grid object"
  )
})

test_that("hexify_lonlat_to_h_index validates input lengths", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  expect_error(
    hexify_lonlat_to_h_index(grid, lon = c(0, 1), lat = 0),
    "same length"
  )
})

test_that("hexify_lonlat_to_h_index validates numeric input", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  expect_error(
    hexify_lonlat_to_h_index(grid, lon = "a", lat = "b"),
    "must be numeric"
  )
})

test_that("hexify_lonlat_to_h_index warns on out-of-range coordinates", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  expect_warning(
    hexify_lonlat_to_h_index(grid, lon = 200, lat = 45),
    "outside valid range"
  )

  expect_warning(
    hexify_lonlat_to_h_index(grid, lon = 0, lat = 100),
    "outside valid range"
  )
})

test_that("hexify_h_index_to_lonlat returns correct structure", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # First get an index
  h_result <- hexify_lonlat_to_h_index(grid, lon = 5, lat = 45)

  # Then convert back
  result <- hexify_h_index_to_lonlat(grid, h_result$h_index)

  expect_s3_class(result, "data.frame")
  expect_true("lon" %in% names(result))
  expect_true("lat" %in% names(result))
  expect_type(result$lon, "double")
  expect_type(result$lat, "double")
})

test_that("hexify_h_index_to_lonlat handles NA values", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # Get a valid h_index first
  valid_index <- hexify_lonlat_to_h_index(grid, lon = 5, lat = 45)$h_index

  result <- hexify_h_index_to_lonlat(grid, c(valid_index, NA_character_))

  expect_equal(nrow(result), 2)
  expect_true(!is.na(result$lon[1]))
  expect_true(is.na(result$lon[2]))
  expect_true(is.na(result$lat[2]))
})

test_that("hexify_h_index_to_lonlat validates grid object", {
  expect_error(
    hexify_h_index_to_lonlat(list(x = 1), "0100000"),
    "hexify_grid object"
  )
})

test_that("hexify_h_index_to_lonlat validates character input", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  expect_error(
    hexify_h_index_to_lonlat(grid, 12345),
    "must be a character"
  )
})

test_that("h_index round-trip is consistent", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  lon <- 10
  lat <- 45

  # Forward
  h_result <- hexify_lonlat_to_h_index(grid, lon = lon, lat = lat)

  # Inverse
  coords <- hexify_h_index_to_lonlat(grid, h_result$h_index)

  # Should be close (within cell distance)
  expect_true(abs(coords$lon - lon) < 5)
  expect_true(abs(coords$lat - lat) < 5)
})

# =============================================================================
# GRID-BASED WRAPPERS
# =============================================================================

test_that("hexify_grid_to_cell returns cell IDs", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(result > 0))
})

test_that("hexify_grid_to_cell validates grid object", {
  expect_error(
    hexify_grid_to_cell(list(x = 1), lon = 0, lat = 0),
    "hexify_grid object"
  )
})

test_that("hexify_grid_cell_to_lonlat returns coordinates", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  # Get some cell IDs first
  cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))

  # Convert back
  result <- hexify_grid_cell_to_lonlat(grid, cell_ids)

  expect_s3_class(result, "data.frame")
  expect_true("lon_deg" %in% names(result))
  expect_true("lat_deg" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify_grid_cell_to_lonlat validates grid object", {
  expect_error(
    hexify_grid_cell_to_lonlat(list(x = 1), cell_id = 1000),
    "hexify_grid object"
  )
})

test_that("grid wrapper round-trip is consistent", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  lon <- 5
  lat <- 45

  # Forward
  cell_id <- hexify_grid_to_cell(grid, lon = lon, lat = lat)

  # Inverse
  coords <- hexify_grid_cell_to_lonlat(grid, cell_id)

  # Should be close
  expect_true(abs(coords$lon_deg - lon) < 5)
  expect_true(abs(coords$lat_deg - lat) < 5)
})

# =============================================================================
# ROUNDTRIP TEST FUNCTION
# =============================================================================

test_that("hexify_roundtrip_test returns correct structure", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_roundtrip_test(grid, lon = 5, lat = 45)

  expect_type(result, "list")
  expect_true("original" %in% names(result))
  expect_true("h_index" %in% names(result))
  expect_true("reconstructed" %in% names(result))
  expect_true("error" %in% names(result))
  expect_true("units" %in% names(result))
})

test_that("hexify_roundtrip_test reports error in km", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_roundtrip_test(grid, lon = 5, lat = 45, units = "km")

  expect_equal(result$units, "km")
  expect_true(is.numeric(result$error))
  expect_true(result$error >= 0)
})

test_that("hexify_roundtrip_test reports error in degrees", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_roundtrip_test(grid, lon = 5, lat = 45, units = "degrees")

  expect_equal(result$units, "degrees")
  expect_true(is.numeric(result$error))
  expect_true(result$error >= 0)
})

test_that("hexify_roundtrip_test error is reasonable", {
  grid <- hexify_grid(area = 1000, aperture = 3)

  result <- hexify_roundtrip_test(grid, lon = 5, lat = 45, units = "km")

  # Error should be less than the cell diagonal (for ~1000 km2 area, diagonal ~ 35-40 km)
  expect_true(result$error < 50)
})

