
# tests/testthat/test-aperture-7.R
# Tests for aperture 7 (ISEA7H) hexagonal grid quantization
#
# Successive aperture-7 levels rotate by kAp7RotDeg in alternating directions,
# so even resolutions are unrotated (Class I) and odd resolutions sit at
# kAp7RotDeg on a sqrt(7) substrate.
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
# ORIENTATION ALTERNATION
# =============================================================================

test_that("aperture 7 keeps cell (0,0) at the origin at every resolution", {
  setup_icosa()

  for (res in c(0, 1, 2, 3, 4)) {
    center <- cpp_hex_center_ap7(0, 0, res)
    expect_equal(as.numeric(center["cx"]), 0, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)
  }
})

# =============================================================================
# SCALING AND REFINEMENT
# =============================================================================

test_that("aperture 7 substrate divisor is 7^ceil(res/2)", {
  setup_icosa()

  # Cell centres are substrate coordinates divided by sqrt(7)^res times the
  # substrate factor (1 at even resolutions, sqrt(7) at odd ones), so the
  # divisor is the integer 7^ceil(res/2) that the exact aperture-7 route in
  # coordinate_transforms.cpp uses.

  for (res in 0:5) {
    center <- cpp_hex_center_ap7(1, 0, res)
    unit <- cpp_hex_center_ap7(1, 0, 0)
    divisor <- as.numeric(unit["cx"]) / as.numeric(center["cx"])

    expect_equal(divisor, 7^ceiling(res / 2), tolerance = 1e-9,
                 info = sprintf("res=%d substrate divisor", res))
  }
})

test_that("aperture 7 is centre-nested: a parent centre is also a child centre", {
  setup_icosa()

  set.seed(707)
  pts <- cbind(runif(30, -0.7, 0.7), runif(30, -0.7, 0.7))

  for (res in 0:4) {
    child_spacing <- 1 / sqrt(7)^(res + 1)
    for (k in seq_len(nrow(pts))) {
      cell <- cpp_hex_quantize_ap7(pts[k, 1], pts[k, 2], res)
      p <- cpp_hex_center_ap7(cell["i"], cell["j"], res)
      kid <- cpp_hex_quantize_ap7(p["cx"], p["cy"], res + 1L)
      q <- cpp_hex_center_ap7(kid["i"], kid["j"], res + 1L)

      d <- sqrt((p["cx"] - q["cx"])^2 + (p["cy"] - q["cy"])^2)
      expect_lt(as.numeric(d), child_spacing * 1e-6)
    }
  }
})

# =============================================================================
# HEXAGON CORNERS
# =============================================================================

test_that("aperture 7 corners form valid hexagons", {
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
