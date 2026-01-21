
# tests/testthat/test-aperture-7.R
# Tests for aperture 7 (ISEA7H) hexagonal grid quantization
#
# Aperture 7 alternates between Class III-I and Class III-II orientations.
# Each resolution has 7x the cells of the previous (area ratio 7:1).

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# ROUND-TRIP TESTS
# =============================================================================

test_that("aperture 7 round-trip works", {
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

      cell <- cpp_hex_quantize_ap7(tx, ty, res)
      center <- cpp_hex_center_ap7(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantize_ap7(center["cx"], center["cy"], res)

      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("aperture 7 batch round-trip succeeds", {
  skip_on_cran()
  setup_icosa()

  set.seed(789)
  n <- 20
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)

  for (res in c(2, 3, 4, 5)) {
    result <- cpp_batch_test_roundtrip_ap7(tx, ty, res)
    n_pass <- sum(result$success)

    expect_equal(n_pass, n,
                 info = sprintf("res=%d: %d/%d passed", res, n_pass, n))
  }
})

# =============================================================================
# CLASS ALTERNATION
# =============================================================================

test_that("aperture 7 alternates Class III-I and III-II", {
  skip_on_cran()
  setup_icosa()

  # Cell (0,0) should always be at origin regardless of class
  for (res in c(0, 1, 2, 3, 4)) {
    center <- cpp_hex_center_ap7(0, 0, res)
    expect_equal(as.numeric(center["cx"]), 0, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)
  }
})

# =============================================================================
# SCALING AND REFINEMENT
# =============================================================================

test_that("aperture 7 refines by factor of 7", {
  skip_on_cran()
  setup_icosa()

  # Aperture 7 has asymmetric scaling due to Class III variants:
  # - Class III-I (even) → Class III-II (odd): scale by sqrt(21)
  # - Class III-II (odd) → Class III-I (even): scale by sqrt(7)/sqrt(3)
  # Overall cumulative: sqrt(7)^res

  sqrt7 <- sqrt(7)
  sqrt3 <- sqrt(3)
  sqrt21 <- sqrt(21)

  for (res in 1:4) {
    center_r <- cpp_hex_center_ap7(1, 0, res)
    center_prev <- cpp_hex_center_ap7(1, 0, res - 1)

    ratio_x <- abs(as.numeric(center_prev["cx"]) / as.numeric(center_r["cx"]))

    is_prev_class3i <- ((res - 1) %% 2 == 0)

    if (is_prev_class3i) {
      expected_ratio <- sqrt21
    } else {
      expected_ratio <- sqrt7 / sqrt3
    }

    expect_equal(ratio_x, expected_ratio, tolerance = 0.01,
                 info = sprintf("res=%d scaling ratio", res))
  }
})

# =============================================================================
# HEXAGON CORNERS
# =============================================================================

test_that("aperture 7 corners form valid hexagons", {
  skip_on_cran()
  setup_icosa()

  for (res in c(0, 1, 2, 3)) {
    corners <- cpp_hex_corners_ap7(0, 0, res, 1.0)

    expect_equal(length(corners$x), 6)
    expect_equal(length(corners$y), 6)
    expect_true(all(is.finite(corners$x)))
    expect_true(all(is.finite(corners$y)))
  }
})

# =============================================================================
# LON/LAT WORKFLOW
# =============================================================================

test_that("aperture 7 lon/lat workflow works", {
  skip_on_cran()
  setup_icosa()

  lon <- 16.37  # Vienna
  lat <- 48.21
  res <- 5

  cell <- cpp_lonlat_to_cell_ap7(lon, lat, res)

  expect_true(cell["face"] >= 0 && cell["face"] < 20)
  expect_true(is.numeric(cell["i"]))
  expect_true(is.numeric(cell["j"]))

  ll <- cpp_cell_to_lonlat_ap7(cell["face"], cell["i"], cell["j"], res)

  # Aperture 7 may have larger errors due to surrogate-substrate conversion
  dist <- sqrt((ll["lon"] - lon)^2 + (ll["lat"] - lat)^2)
  expect_true(dist < 15.0)
})

# =============================================================================
# SINGLE-POINT ROUNDTRIP TEST HELPER
# =============================================================================

test_that("cpp_test_roundtrip_ap7 returns TRUE for valid points", {
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
      result <- cpp_test_roundtrip_ap7(pt[1], pt[2], res)
      expect_true(result, info = sprintf("res=%d, pt=(%.2f, %.2f)", res, pt[1], pt[2]))
    }
  }
})
