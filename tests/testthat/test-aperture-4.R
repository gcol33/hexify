
# tests/testthat/test-aperture-4.R
# Tests for aperture 4 (ISEA4H) hexagonal grid quantization
#
# Aperture 4 uses Class I (flat-top) at all resolutions.
# Each resolution has 4x the cells of the previous (area ratio 4:1).

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# ROUND-TRIP TESTS
# =============================================================================

test_that("aperture 4 round-trip works", {
  skip_on_cran()
  setup_icosa()

  test_points <- list(
    c(0.4, 0.35),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0),
    c(0.7, 0.5)
  )

  for (res in c(0, 1, 2, 3, 4, 5)) {
    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell <- cpp_hex_quantize_ap4(tx, ty, res)
      center <- cpp_hex_center_ap4(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantize_ap4(center["cx"], center["cy"], res)

      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("aperture 4 batch round-trip succeeds", {
  skip_on_cran()
  setup_icosa()

  set.seed(456)
  n <- 20
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)

  for (res in c(2, 3, 4, 5)) {
    result <- cpp_batch_test_roundtrip_ap4(tx, ty, res)
    n_pass <- sum(result$success)

    expect_equal(n_pass, n,
                 info = sprintf("res=%d: %d/%d passed", res, n_pass, n))
  }
})

# =============================================================================
# CLASS I CONSISTENCY
# =============================================================================

test_that("aperture 4 always uses Class I (no rotation)", {
  skip_on_cran()
  setup_icosa()

  for (res in c(0, 1, 2, 3, 4)) {
    # Cell (0,0) should always be at origin
    center <- cpp_hex_center_ap4(0, 0, res)
    expect_equal(as.numeric(center["cx"]), 0, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)

    # Cell (1,0) should follow Class I formula
    center <- cpp_hex_center_ap4(1, 0, res)
    expected_x <- 1.0 / (2^res)
    expect_equal(as.numeric(center["cx"]), expected_x, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)

    # Cell (0,1) should follow Class I formula: x = -0.5/scale
    center <- cpp_hex_center_ap4(0, 1, res)
    expected_x <- -0.5 / (2^res)
    expected_y <- sin(60 * pi / 180) / (2^res)
    expect_equal(as.numeric(center["cx"]), expected_x, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), expected_y, tolerance = 1e-10)
  }
})

# =============================================================================
# SCALING AND REFINEMENT
# =============================================================================

test_that("aperture 4 refines by factor of 4", {
  skip_on_cran()
  setup_icosa()

  # Scale should double each resolution (2^res), area quarters
  for (res in 0:4) {
    scale <- 2^res

    center <- cpp_hex_center_ap4(1, 0, res)
    expected_x <- 1.0 / scale

    expect_equal(as.numeric(center["cx"]), expected_x, tolerance = 1e-10,
                 info = sprintf("res=%d: cell (1,0) x-coord", res))
  }
})

test_that("aperture 4 cell area ratios are correct", {
  skip_on_cran()
  setup_icosa()

  # Adjacent resolutions should have area ratio of 4:1
  for (res in 1:4) {
    # Use cell spacing as proxy for area
    center_cur <- cpp_hex_center_ap4(1, 0, res)
    center_prev <- cpp_hex_center_ap4(1, 0, res - 1)

    spacing_cur <- as.numeric(center_cur["cx"])
    spacing_prev <- as.numeric(center_prev["cx"])

    ratio <- spacing_prev / spacing_cur
    expect_equal(ratio, 2, tolerance = 1e-10,
                 info = sprintf("res=%d spacing ratio", res))
  }
})

# =============================================================================
# HEXAGON CORNERS
# =============================================================================

test_that("aperture 4 corners form valid hexagons", {
  skip_on_cran()
  setup_icosa()

  for (res in c(0, 1, 2, 3)) {
    corners <- cpp_hex_corners_ap4(0, 0, res, 1.0)

    expect_equal(length(corners$x), 6)
    expect_equal(length(corners$y), 6)
    expect_true(all(is.finite(corners$x)))
    expect_true(all(is.finite(corners$y)))
  }
})

# =============================================================================
# LON/LAT WORKFLOW
# =============================================================================

test_that("aperture 4 lon/lat workflow works", {
  skip_on_cran()
  setup_icosa()

  lon <- 16.37  # Vienna
  lat <- 48.21
  res <- 6

  cell <- cpp_lonlat_to_cell_ap4(lon, lat, res)

  expect_true(cell["face"] >= 0 && cell["face"] < 20)
  expect_true(is.numeric(cell["i"]))
  expect_true(is.numeric(cell["j"]))

  ll <- cpp_cell_to_lonlat_ap4(cell["face"], cell["i"], cell["j"], res)

  # Should be reasonably close
  dist <- sqrt((ll["lon"] - lon)^2 + (ll["lat"] - lat)^2)
  expect_true(dist < 10.0)
})

# =============================================================================
# SINGLE-POINT ROUNDTRIP TEST HELPER
# =============================================================================

test_that("cpp_test_roundtrip_ap4 returns TRUE for valid points", {
  skip_on_cran()
  setup_icosa()

  test_points <- list(
    c(0.5, 0.3),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0)
  )

  for (res in c(2, 3, 4, 5)) {
    for (pt in test_points) {
      result <- cpp_test_roundtrip_ap4(pt[1], pt[2], res)
      expect_true(result, info = sprintf("res=%d, pt=(%.2f, %.2f)", res, pt[1], pt[2]))
    }
  }
})
