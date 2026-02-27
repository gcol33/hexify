# tests/testthat/test-aperture-3.R
# Tests for aperture 3 (ISEA3H) hexagonal grid quantization
#
# Aperture 3 uses Class I (flat-top) at even resolutions and
# Class II (pointy-top, rotated 30°) at odd resolutions.

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# ROUND-TRIP TESTS
# =============================================================================

test_that("aperture 3 round-trip works for Class I (even resolutions)", {
  skip_on_cran()
  setup_icosa()

  test_points <- list(
    c(0.5, 0.3),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0)
  )

  for (res in c(0, 2, 4, 6)) {
    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell <- cpp_hex_quantize_ap3(tx, ty, res)
      center <- cpp_hex_center_ap3(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantize_ap3(center["cx"], center["cy"], res)

      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("Class I res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("Class I res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("aperture 3 round-trip works for Class II (odd resolutions)", {
  skip_on_cran()
  setup_icosa()

  test_points <- list(
    c(0.5, 0.3),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0)
  )

  for (res in c(1, 3, 5)) {
    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell <- cpp_hex_quantize_ap3(tx, ty, res)
      center <- cpp_hex_center_ap3(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantize_ap3(center["cx"], center["cy"], res)

      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("Class II res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("Class II res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("aperture 3 batch round-trip succeeds", {
  skip_on_cran()
  setup_icosa()

  set.seed(123)
  n <- 20
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)

  for (res in c(2, 3, 4, 5)) {
    result <- cpp_batch_test_roundtrip_ap3(tx, ty, res)
    n_pass <- sum(result$success)

    expect_equal(n_pass, n,
                 info = sprintf("res=%d: %d/%d passed", res, n_pass, n))
  }
})

# =============================================================================
# SCALING AND REFINEMENT
# =============================================================================

test_that("aperture 3 cell spacing decreases with resolution", {
  skip_on_cran()
  setup_icosa()

  # Compare spacing between same-class resolutions
  # Class I: 0, 2, 4, 6 (even)
  # Class II: 1, 3, 5, 7 (odd)

  for (res in c(2, 4)) {
    center_lo_0 <- cpp_hex_center_ap3(0, 0, res)
    center_lo_1 <- cpp_hex_center_ap3(1, 0, res)
    spacing_lo <- abs(center_lo_1["cx"] - center_lo_0["cx"])

    center_hi_0 <- cpp_hex_center_ap3(0, 0, res + 2)
    center_hi_1 <- cpp_hex_center_ap3(1, 0, res + 2)
    spacing_hi <- abs(center_hi_1["cx"] - center_hi_0["cx"])

    expect_true(spacing_hi < spacing_lo,
                info = sprintf("res=%d vs res=%d spacing", res, res + 2))
  }
})

test_that("resolution refinement decreases distance to point", {
  skip_on_cran()
  setup_icosa()

  tx <- 0.4
  ty <- 0.35

  cell2 <- cpp_hex_quantize_ap3(tx, ty, 2)
  cell4 <- cpp_hex_quantize_ap3(tx, ty, 4)
  cell6 <- cpp_hex_quantize_ap3(tx, ty, 6)

  center2 <- cpp_hex_center_ap3(cell2["i"], cell2["j"], 2)
  center4 <- cpp_hex_center_ap3(cell4["i"], cell4["j"], 4)
  center6 <- cpp_hex_center_ap3(cell6["i"], cell6["j"], 6)

  dist2 <- sqrt((center2["cx"] - tx)^2 + (center2["cy"] - ty)^2)
  dist4 <- sqrt((center4["cx"] - tx)^2 + (center4["cy"] - ty)^2)
  dist6 <- sqrt((center6["cx"] - tx)^2 + (center6["cy"] - ty)^2)

  expect_true(dist4 < dist2)
  expect_true(dist6 < dist4)
})

# =============================================================================
# CLASS I TO CLASS II TRANSITION
# =============================================================================

test_that("Class I to Class II transition is consistent", {
  skip_on_cran()
  setup_icosa()

  tx <- 0.42
  ty <- 0.41

  cell2 <- cpp_hex_quantize_ap3(tx, ty, 2)  # Class I
  cell3 <- cpp_hex_quantize_ap3(tx, ty, 3)  # Class II
  cell4 <- cpp_hex_quantize_ap3(tx, ty, 4)  # Class I

  center2 <- cpp_hex_center_ap3(cell2["i"], cell2["j"], 2)
  center3 <- cpp_hex_center_ap3(cell3["i"], cell3["j"], 3)
  center4 <- cpp_hex_center_ap3(cell4["i"], cell4["j"], 4)

  dist2 <- sqrt((center2["cx"] - tx)^2 + (center2["cy"] - ty)^2)
  dist3 <- sqrt((center3["cx"] - tx)^2 + (center3["cy"] - ty)^2)
  dist4 <- sqrt((center4["cx"] - tx)^2 + (center4["cy"] - ty)^2)

  expect_true(dist3 < dist2)
  expect_true(dist4 < dist3)
})

# =============================================================================
# HEXAGON CORNERS
# =============================================================================

test_that("aperture 3 corners form valid hexagons", {
  skip_on_cran()
  setup_icosa()

  for (res in c(2, 3)) {
    corners <- cpp_hex_corners_ap3(0, 0, res, 1.0)

    expect_equal(length(corners$x), 6)
    expect_equal(length(corners$y), 6)
    expect_true(all(is.finite(corners$x)))
    expect_true(all(is.finite(corners$y)))
  }
})

test_that("hexagon corners centroid matches center", {
  skip_on_cran()
  setup_icosa()

  for (res in c(2, 4)) {
    for (i in 0:2) {
      for (j in 0:2) {
        center <- cpp_hex_center_ap3(i, j, res)
        corners <- cpp_hex_corners_ap3(i, j, res, 1.0)

        centroid_x <- mean(corners$x)
        centroid_y <- mean(corners$y)

        expect_equal(centroid_x, as.numeric(center["cx"]), tolerance = 1e-10)
        expect_equal(centroid_y, as.numeric(center["cy"]), tolerance = 1e-10)
      }
    }
  }
})

# =============================================================================
# LON/LAT WORKFLOW
# =============================================================================

test_that("aperture 3 lon/lat workflow works", {
  skip_on_cran()
  setup_icosa()

  lon <- 16.37  # Vienna
  lat <- 48.21
  res <- 6

  cell <- cpp_lonlat_to_cell_ap3(lon, lat, res)

  expect_true(cell["face"] >= 0 && cell["face"] < 20)
  expect_true(is.numeric(cell["i"]))
  expect_true(is.numeric(cell["j"]))

  ll <- cpp_cell_to_lonlat_ap3(cell["face"], cell["i"], cell["j"], res)

  # Should be reasonably close (within cell diameter)
  dist <- sqrt((ll["lon"] - lon)^2 + (ll["lat"] - lat)^2)
  expect_true(dist < 10.0)
})

# =============================================================================
# SINGLE-POINT ROUNDTRIP TEST HELPER
# =============================================================================

test_that("cpp_test_roundtrip_ap3 returns TRUE for valid points", {
  skip_on_cran()
  setup_icosa()

  # Test various points and resolutions
  test_points <- list(
    c(0.5, 0.3),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0)
  )

  for (res in c(2, 3, 4, 5)) {
    for (pt in test_points) {
      result <- cpp_test_roundtrip_ap3(pt[1], pt[2], res)
      expect_true(result, info = sprintf("res=%d, pt=(%.2f, %.2f)", res, pt[1], pt[2]))
    }
  }
})
