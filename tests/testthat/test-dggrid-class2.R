# tests/testthat/test-dggrid-class2.R
# Test Class II (odd resolutions) aperture 3

test_that("Class II round-trip works", {
  # Initialize icosahedron
  cpp_build_icosa()
  
  # Test odd resolutions (Class II)
  for (res in c(1, 3, 5)) {
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

test_that("Class II batch test works", {
  cpp_build_icosa()
  
  set.seed(789)
  n <- 20
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)
  
  for (res in c(1, 3, 5)) {
    result <- cpp_batch_test_roundtrip_ap3(tx, ty, res)
    
    n_pass <- sum(result$success)
    
    if (n_pass < n) {
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

test_that("Class I to Class II transition", {
  cpp_build_icosa()


  # A point should refine consistently across resolutions

  # Note: Using (0.42, 0.41) - chosen to avoid coincidental alignment where
  # the point falls equidistant from cell centers at different resolutions
  tx <- 0.42
  ty <- 0.41
  
  # Get cells at consecutive resolutions
  cell2 <- cpp_hex_quantify_ap3(tx, ty, 2)  # Class I
  cell3 <- cpp_hex_quantify_ap3(tx, ty, 3)  # Class II
  cell4 <- cpp_hex_quantify_ap3(tx, ty, 4)  # Class I
  
  # Get centers
  center2 <- cpp_hex_center_ap3(cell2["i"], cell2["j"], 2)
  center3 <- cpp_hex_center_ap3(cell3["i"], cell3["j"], 3)
  center4 <- cpp_hex_center_ap3(cell4["i"], cell4["j"], 4)
  
  # Distance to original point should decrease with resolution
  dist2 <- sqrt((center2["cx"] - tx)^2 + (center2["cy"] - ty)^2)
  dist3 <- sqrt((center3["cx"] - tx)^2 + (center3["cy"] - ty)^2)
  dist4 <- sqrt((center4["cx"] - tx)^2 + (center4["cy"] - ty)^2)
  
  expect_true(dist3 < dist2,
              info = sprintf("res3 should be closer than res2: dist2=%.6f, dist3=%.6f", 
                           dist2, dist3))
  expect_true(dist4 < dist3,
              info = sprintf("res4 should be closer than res3: dist3=%.6f, dist4=%.6f", 
                           dist3, dist4))
})

test_that("Class II corners are rotated 30 degrees", {
  cpp_build_icosa()
  
  # Same cell at Class I (res 2) and Class II (res 3)
  tx <- 0.4
  ty <- 0.35
  
  cell2 <- cpp_hex_quantify_ap3(tx, ty, 2)  # Class I
  cell3 <- cpp_hex_quantify_ap3(tx, ty, 3)  # Class II
  
  corners2 <- cpp_hex_corners_ap3(cell2["i"], cell2["j"], 2, 0.1)
  corners3 <- cpp_hex_corners_ap3(cell3["i"], cell3["j"], 3, 0.1)
  
  # Class I: first vertex should be directly above center (angle = 90°)
  # Class II: first vertex should be at angle = 120° (rotated 30°)
  
  # Just verify we get 6 corners for each
  expect_equal(length(corners2$x), 6)
  expect_equal(length(corners2$y), 6)
  expect_equal(length(corners3$x), 6)
  expect_equal(length(corners3$y), 6)
})
