# tests/testthat/test-dggrid-class1.R
# Test Class I (even resolutions) aperture 3

test_that("Class I round-trip works", {
  # Initialize icosahedron
  cpp_build_icosa()

  # Test even resolutions (Class I)
  for (res in c(0, 2, 4, 6)) {
    # Test a few points
    test_points <- list(
      c(0.5, 0.3),
      c(-0.4, 0.2),
      c(0.1, -0.6),
      c(0.0, 0.0)
    )

    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      # Quantize to cell
      cell <- cpp_hex_quantify_ap3(tx, ty, res)

      # Get center
      center <- cpp_hex_center_ap3(cell["i"], cell["j"], res)

      # Re-quantize center
      cell2 <- cpp_hex_quantify_ap3(center["cx"], center["cy"], res)

      # Should match
      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("Class I batch test works", {
  cpp_build_icosa()

  set.seed(123)
  n <- 20  # Start small
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)

  for (res in c(2, 4)) {
    result <- cpp_batch_test_roundtrip_ap3(tx, ty, res)

    n_pass <- sum(result$success)

    if (n_pass < n) {
      # Show failures
      fail_idx <- which(!result$success)
      message(sprintf("\nResolution %d: %d/%d passed", res, n_pass, n))
      for (idx in head(fail_idx, 3)) {
        message(sprintf("  Failed: tx=%.6f, ty=%.6f", tx[idx], ty[idx]))
        message(sprintf("    orig: (i=%d, j=%d)",
                        result$i_orig[idx], result$j_orig[idx]))
        message(sprintf("    recomp: (i=%d, j=%d)",
                        result$i_recomp[idx], result$j_recomp[idx]))
      }
    }

    expect_equal(n_pass, n,
                 info = sprintf("res=%d: %d/%d passed", res, n_pass, n))
  }
})

test_that("Lon/lat workflow works", {
  cpp_build_icosa()

  # Test a known location
  lon <- 16.37  # Vienna
  lat <- 48.21
  res <- 6

  # Convert to cell
  cell <- cpp_lonlat_to_cell_ap3(lon, lat, res)

  expect_true(cell["face"] >= 0 && cell["face"] < 20)
  expect_true(is.numeric(cell["i"]))
  expect_true(is.numeric(cell["j"]))

  # Convert back
  ll <- cpp_cell_to_lonlat_ap3(cell["face"], cell["i"], cell["j"], res)

  # Should be reasonably close (within a few degrees at res 6)
  dist <- sqrt((ll["lon"] - lon)^2 + (ll["lat"] - lat)^2)
  expect_true(dist < 10.0,
              info = sprintf("Distance: %.3f degrees", dist))
})

test_that("Centers are near cell centroids", {
  cpp_build_icosa()

  set.seed(456)

  for (res in c(2, 4)) {
    # Test a few random cells
    for (k in 1:5) {
      tx <- runif(1, -0.7, 0.7)
      ty <- runif(1, -0.7, 0.7)

      cell <- cpp_hex_quantify_ap3(tx, ty, res)
      center <- cpp_hex_center_ap3(cell["i"], cell["j"], res)

      # Get corners
      corners <- cpp_hex_corners_ap3(cell["i"], cell["j"], res, 1.0)

      # Compute centroid of corners
      centroid_x <- mean(corners$x)
      centroid_y <- mean(corners$y)

      # Center should be close to centroid
      dist <- sqrt((center["cx"] - centroid_x)^2 +
                     (center["cy"] - centroid_y)^2)

      expect_true(dist < 0.01,
                  info = sprintf("res=%d: center-centroid dist=%.6f", res, dist))
    }
  }
})

test_that("Resolution refines correctly", {
  cpp_build_icosa()

  # Same point at different resolutions
  # Note: (0.5, 0.3) was a poor choice as it sits on aligned cell centers
  # Using (0.4, 0.35) instead for reliable refinement testing
  tx <- 0.4
  ty <- 0.35

  # Get cells
  cell2 <- cpp_hex_quantify_ap3(tx, ty, 2)
  cell4 <- cpp_hex_quantify_ap3(tx, ty, 4)
  cell6 <- cpp_hex_quantify_ap3(tx, ty, 6)

  # Get centers
  center2 <- cpp_hex_center_ap3(cell2["i"], cell2["j"], 2)
  center4 <- cpp_hex_center_ap3(cell4["i"], cell4["j"], 4)
  center6 <- cpp_hex_center_ap3(cell6["i"], cell6["j"], 6)

  # Distance to original point should decrease
  dist2 <- sqrt((center2["cx"] - tx)^2 + (center2["cy"] - ty)^2)
  dist4 <- sqrt((center4["cx"] - tx)^2 + (center4["cy"] - ty)^2)
  dist6 <- sqrt((center6["cx"] - tx)^2 + (center6["cy"] - ty)^2)

  expect_true(dist4 < dist2,
              info = sprintf("dist2=%.6f, dist4=%.6f", dist2, dist4))
  expect_true(dist6 < dist4,
              info = sprintf("dist4=%.6f, dist6=%.6f", dist4, dist6))
})
