
# tests/testthat/test-dggrid-compat.R
# Tests verifying hexify produces identical output to dggridR
#
# These tests compare hexify results against dggridR reference data to ensure
# compatibility. The reference data was generated using dggridR.

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

load_validation_data <- function(aperture, resolution) {
  file <- sprintf("data/dggrid_validation_ap%d_res%d.csv", aperture, resolution)
  read.csv(test_path(file))
}

# =============================================================================
# DGGRID VALIDATION TESTS - APERTURE 3
# =============================================================================

test_that("hexify matches dggridR for aperture 3, resolution 3", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 3)

  for (i in seq_len(nrow(ref))) {
    result <- hexify_lonlat_to_quad_ij(
      lon = ref$lon[i],
      lat = ref$lat[i],
      resolution = 3,
      aperture = 3
    )

    expect_equal(result$quad, ref$quad[i],
                 info = sprintf("Row %d: quad mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$i, ref$i[i],
                 info = sprintf("Row %d: i mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$j, ref$j[i],
                 info = sprintf("Row %d: j mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
  }
})

test_that("hexify matches dggridR for aperture 3, resolution 5", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 5)

  for (i in seq_len(nrow(ref))) {
    result <- hexify_lonlat_to_quad_ij(
      lon = ref$lon[i],
      lat = ref$lat[i],
      resolution = 5,
      aperture = 3
    )

    expect_equal(result$quad, ref$quad[i],
                 info = sprintf("Row %d: quad mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$i, ref$i[i],
                 info = sprintf("Row %d: i mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$j, ref$j[i],
                 info = sprintf("Row %d: j mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
  }
})

test_that("hexify matches dggridR for aperture 3, resolution 7", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 7)

  for (i in seq_len(nrow(ref))) {
    result <- hexify_lonlat_to_quad_ij(
      lon = ref$lon[i],
      lat = ref$lat[i],
      resolution = 7,
      aperture = 3
    )

    expect_equal(result$quad, ref$quad[i],
                 info = sprintf("Row %d: quad mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$i, ref$i[i],
                 info = sprintf("Row %d: i mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(result$j, ref$j[i],
                 info = sprintf("Row %d: j mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
  }
})

# =============================================================================
# VECTORIZED CELL ID TESTS
# =============================================================================

test_that("hexify batch cell_id matches dggridR for aperture 3, resolution 3", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 3)

  # Get hexify cell IDs (vectorized)
  hexify_cells <- hexify_lonlat_to_cell(
    lon = ref$lon,
    lat = ref$lat,
    resolution = 3,
    aperture = 3
  )

  # Convert dggridR quad/i/j to cell ID for comparison
  dggrid_cells <- hexify_quad_ij_to_cell(
    quad = ref$quad,
    i = ref$i,
    j = ref$j,
    resolution = 3,
    aperture = 3
  )

  expect_equal(hexify_cells, dggrid_cells)
})

test_that("hexify batch cell_id matches dggridR for aperture 3, resolution 5", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 5)

  hexify_cells <- hexify_lonlat_to_cell(
    lon = ref$lon,
    lat = ref$lat,
    resolution = 5,
    aperture = 3
  )

  dggrid_cells <- hexify_quad_ij_to_cell(
    quad = ref$quad,
    i = ref$i,
    j = ref$j,
    resolution = 5,
    aperture = 3
  )

  expect_equal(hexify_cells, dggrid_cells)
})

test_that("hexify batch cell_id matches dggridR for aperture 3, resolution 7", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 7)

  hexify_cells <- hexify_lonlat_to_cell(
    lon = ref$lon,
    lat = ref$lat,
    resolution = 7,
    aperture = 3
  )

  dggrid_cells <- hexify_quad_ij_to_cell(
    quad = ref$quad,
    i = ref$i,
    j = ref$j,
    resolution = 7,
    aperture = 3
  )

  expect_equal(hexify_cells, dggrid_cells)
})

# =============================================================================
# CELL ID VALIDATION TESTS
# =============================================================================

test_that("hexify cell_id matches dggridR cell IDs for aperture 3", {
  skip_on_cran()
  setup_icosa()

  for (res in c(3, 5, 7)) {
    ref <- load_validation_data(3, res)

    # Get hexify cell IDs
    hexify_cells <- hexify_lonlat_to_cell(
      lon = ref$lon,
      lat = ref$lat,
      resolution = res,
      aperture = 3
    )

    # Convert dggridR quad/i/j to cell ID for comparison
    dggrid_cells <- hexify_quad_ij_to_cell(
      quad = ref$quad,
      i = ref$i,
      j = ref$j,
      resolution = res,
      aperture = 3
    )

    expect_equal(hexify_cells, dggrid_cells,
                 info = sprintf("Cell ID mismatch for aperture 3, resolution %d", res))
  }
})

# =============================================================================
# ICOSAHEDRON FACE VALIDATION
# =============================================================================

test_that("hexify face assignment matches dggridR tnum", {
  skip_on_cran()
  setup_icosa()

  for (res in c(3, 5, 7)) {
    ref <- load_validation_data(3, res)

    for (i in seq_len(nrow(ref))) {
      face <- cpp_which_face(ref$lon[i], ref$lat[i])

      expect_equal(face, ref$tnum[i],
                   info = sprintf("Res %d, row %d: face mismatch for lon=%.4f, lat=%.4f (got %d, expected %d)",
                                  res, i, ref$lon[i], ref$lat[i], face, ref$tnum[i]))
    }
  }
})

# =============================================================================
# PROJECTION VALIDATION (tx, ty)
# =============================================================================

test_that("hexify projection matches dggridR tx/ty", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 5)  # Use res 5 data

  for (i in seq_len(nrow(ref))) {
    proj <- cpp_snyder_forward(ref$lon[i], ref$lat[i])
    face <- unname(proj["face"])
    tx <- unname(proj["icosa_triangle_x"])
    ty <- unname(proj["icosa_triangle_y"])

    # Face should match
    expect_equal(face, ref$tnum[i],
                 info = sprintf("Row %d: face mismatch", i))

    # tx/ty should be close (allow tolerance for floating point differences)
    expect_equal(tx, ref$tx[i], tolerance = 1e-6,
                 info = sprintf("Row %d: tx mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
    expect_equal(ty, ref$ty[i], tolerance = 1e-6,
                 info = sprintf("Row %d: ty mismatch for lon=%.4f, lat=%.4f",
                                i, ref$lon[i], ref$lat[i]))
  }
})

# =============================================================================
# ROUND-TRIP CONSISTENCY TESTS
# =============================================================================

test_that("lon/lat -> cell -> lon/lat round-trip lands in same cell", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 5)

  for (res in c(3, 5, 7)) {
    # Get cell IDs
    cells <- hexify_lonlat_to_cell(
      lon = ref$lon,
      lat = ref$lat,
      resolution = res,
      aperture = 3
    )

    # Get cell centers
    centers <- hexify_cell_to_lonlat(
      cell_id = cells,
      resolution = res,
      aperture = 3
    )

    # Convert centers back to cells
    cells2 <- hexify_lonlat_to_cell(
      lon = centers$lon_deg,
      lat = centers$lat_deg,
      resolution = res,
      aperture = 3
    )

    expect_equal(cells, cells2,
                 info = sprintf("Round-trip failed for resolution %d", res))
  }
})

test_that("quad/i/j -> cell -> quad/i/j round-trip is exact", {
  skip_on_cran()
  setup_icosa()

  ref <- load_validation_data(3, 5)

  for (res in c(3, 5, 7)) {
    ref_res <- load_validation_data(3, res)

    # quad/i/j -> cell
    cells <- hexify_quad_ij_to_cell(
      quad = ref_res$quad,
      i = ref_res$i,
      j = ref_res$j,
      resolution = res,
      aperture = 3
    )

    # cell -> quad/i/j
    result <- hexify_cell_to_quad_ij(
      cell_id = cells,
      resolution = res,
      aperture = 3
    )

    expect_equal(result$quad, ref_res$quad,
                 info = sprintf("Quad mismatch for resolution %d", res))
    expect_equal(result$i, ref_res$i,
                 info = sprintf("i mismatch for resolution %d", res))
    expect_equal(result$j, ref_res$j,
                 info = sprintf("j mismatch for resolution %d", res))
  }
})
