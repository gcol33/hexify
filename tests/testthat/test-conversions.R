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
