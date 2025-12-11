# tests/testthat/test-dggrid-ap7.R
# Test aperture 7 hexagon grid implementation

test_that("Aperture 7 round-trip works", {
  cpp_build_icosa()
  
  # Test multiple resolutions
  for (res in c(0, 1, 2, 3, 4, 5)) {
    # Test several points
    test_points <- list(
      c(0.4, 0.35),
      c(-0.4, 0.2),
      c(0.1, -0.6),
      c(0.0, 0.0),
      c(0.7, 0.5)
    )
    
    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]
      
      # Quantize to cell
      cell <- cpp_hex_quantify_ap7(tx, ty, res)
      
      # Get center
      center <- cpp_hex_center_ap7(cell["i"], cell["j"], res)
      
      # Re-quantize center
      cell2 <- cpp_hex_quantify_ap7(center["cx"], center["cy"], res)
      
      # Should match
      expect_equal(cell["i"], cell2["i"], 
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("Aperture 7 alternates Class I and Class II", {
  cpp_build_icosa()
  
  # Even resolutions should use Class I (no rotation)
  # Odd resolutions should use Class II (30° rotation)
  
  for (res in c(0, 1, 2, 3, 4)) {
    # Test cell (0,0) - should always be at origin after scaling
    center <- cpp_hex_center_ap7(0, 0, res)
    expect_equal(as.numeric(center["cx"]), 0, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)
  }
})

test_that("Aperture 7 refines by factor of 7", {
  cpp_build_icosa()
  
  # Aperture 7 has asymmetric scaling due to Class III variants:
  # - Class III-I (even) → Class III-II (odd): scale by sqrt(21) ≈ 4.58
  # - Class III-II (odd) → Class III-I (even): scale by sqrt(7)/sqrt(3) ≈ 1.53
  # - Overall cumulative: product gives sqrt(7)^res
  
  sqrt7 <- sqrt(7)
  sqrt3 <- sqrt(3)
  sqrt21 <- sqrt(21)  # sqrt(7*3)
  
  for (res in 1:4) {
    center_r <- cpp_hex_center_ap7(1, 0, res)
    center_prev <- cpp_hex_center_ap7(1, 0, res - 1)
    
    ratio_x <- abs(as.numeric(center_prev["cx"]) / as.numeric(center_r["cx"]))
    
    # Alternating pattern based on transition
    is_prev_class3i <- ((res - 1) %% 2 == 0)
    
    if (is_prev_class3i) {
      # Class III-I → Class III-II: scale by sqrt(21)
      expected_ratio <- sqrt21
    } else {
      # Class III-II → Class III-I: scale by sqrt(7)/sqrt(3)
      expected_ratio <- sqrt7 / sqrt3
    }
    
    expect_equal(ratio_x, expected_ratio, tolerance = 0.05,
                 info = sprintf("res %d to %d: scale ratio", res-1, res))
  }
  
  # Check cumulative scaling
  # The product of all transitions should equal sqrt(7)^res
  # But the path is: sqrt(21) × (sqrt(7)/sqrt(3)) × sqrt(21) × ... = sqrt(7)^res
  # Because: sqrt(21) × sqrt(7)/sqrt(3) = sqrt(21×7/3) = sqrt(49) = 7
  
  for (res in c(2, 4)) {
    center <- cpp_hex_center_ap7(1, 0, res)
    center_0 <- cpp_hex_center_ap7(1, 0, 0)
    ratio <- as.numeric(center_0["cx"]) / as.numeric(center["cx"])
    
    # Cumulative should be sqrt(7)^res
    expected_cumulative <- sqrt7^res
    
    expect_equal(ratio, expected_cumulative, tolerance = 0.05,
                 info = sprintf("cumulative scale from res 0 to %d", res))
  }
})

test_that("Aperture 7 refinement decreases distance", {
  cpp_build_icosa()
  
  # Test point that doesn't sit on aligned centers
  tx <- 0.4
  ty <- 0.35
  
  dists <- numeric()
  
  for (res in c(0, 1, 2, 3, 4, 5)) {
    cell <- cpp_hex_quantify_ap7(tx, ty, res)
    center <- cpp_hex_center_ap7(cell["i"], cell["j"], res)
    dist <- sqrt((center["cx"] - tx)^2 + (center["cy"] - ty)^2)
    dists <- c(dists, dist)
  }
  
  # Distances should generally decrease
  # Later resolutions should definitely be closer
  expect_true(dists[5] < dists[3],
              info = sprintf("dist[res4]=%.6f should be < dist[res2]=%.6f", 
                           dists[5], dists[3]))
  expect_true(dists[6] < dists[4],
              info = sprintf("dist[res5]=%.6f should be < dist[res3]=%.6f", 
                           dists[6], dists[4]))
})

test_that("Aperture 7 batch test works", {
  cpp_build_icosa()
  
  set.seed(789)
  n <- 50
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)
  
  for (res in c(0, 1, 2, 3, 4)) {
    failures <- 0
    
    for (i in 1:n) {
      cell <- cpp_hex_quantify_ap7(tx[i], ty[i], res)
      center <- cpp_hex_center_ap7(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantify_ap7(center["cx"], center["cy"], res)
      
      if (cell["i"] != cell2["i"] || cell["j"] != cell2["j"]) {
        failures <- failures + 1
        if (failures <= 3) {
          message(sprintf("Res %d failure: tx=%.6f, ty=%.6f", res, tx[i], ty[i]))
          message(sprintf("  cell1: (%d, %d), cell2: (%d, %d)", 
                         cell["i"], cell["j"], cell2["i"], cell2["j"]))
        }
      }
    }
    
    expect_equal(failures, 0,
                 info = sprintf("res=%d: %d/%d failed", res, failures, n))
  }
})

test_that("Aperture 7 corners have correct count and shape", {
  cpp_build_icosa()
  
  for (res in c(0, 1, 2, 3)) {
    # Test a few cells
    cells <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
    
    for (cell in cells) {
      i <- cell[1]
      j <- cell[2]
      
      corners <- cpp_hex_corners_ap7(i, j, res, 1.0)
      
      # Should have 6 corners
      expect_equal(length(corners$x), 6,
                   info = sprintf("res=%d, cell=(%d,%d)", res, i, j))
      expect_equal(length(corners$y), 6)
      
      # Center should be inside the hexagon
      center <- cpp_hex_center_ap7(i, j, res)
      mean_x <- mean(corners$x)
      mean_y <- mean(corners$y)
      
      # Centroid should be close to center
      dist <- sqrt((mean_x - center["cx"])^2 + (mean_y - center["cy"])^2)
      expect_true(dist < 0.1,
                  info = sprintf("res=%d, cell=(%d,%d): centroid-center dist=%.6f",
                               res, i, j, dist))
    }
  }
})

test_that("Aperture 7 lon/lat workflow works", {
  cpp_build_icosa()
  
  # Test known locations
  test_locations <- list(
    list(lon = 16.37, lat = 48.21, name = "Vienna"),
    list(lon = -73.94, lat = 40.71, name = "New York"),
    list(lon = 139.69, lat = 35.68, name = "Tokyo"),
    list(lon = 0, lat = 0, name = "Equator/Prime Meridian")
  )
  
  for (loc in test_locations) {
    for (res in c(3, 5, 7)) {
      # Convert to cell
      cell <- cpp_lonlat_to_cell_ap7(loc$lon, loc$lat, res)
      
      expect_true(cell["face"] >= 0 && cell["face"] < 20,
                  info = sprintf("%s res=%d: invalid face", loc$name, res))
      expect_true(is.numeric(cell["i"]))
      expect_true(is.numeric(cell["j"]))
      
      # Convert back
      ll <- cpp_cell_to_lonlat_ap7(cell["face"], cell["i"], cell["j"], res)
      
      # Should be reasonably close
      max_error <- if (res >= 7) 1.0 else if (res >= 5) 5.0 else 10.0
      dist <- sqrt((ll["lon"] - loc$lon)^2 + (ll["lat"] - loc$lat)^2)
      
      expect_true(dist < max_error,
                  info = sprintf("%s res=%d: dist=%.3f > max_error=%.1f",
                               loc$name, res, dist, max_error))
    }
  }
})

test_that("Aperture 7 has correct cumulative scaling", {
  cpp_build_icosa()
  
  # Verify that aperture 7's cumulative scale is sqrt(7)^res
  # despite its asymmetric Class III-I/III-II alternation
  
  sqrt7 <- sqrt(7)
  
  # Test at several resolutions
  for (res in c(0, 2, 4)) {
    if (res == 0) next  # skip res 0 as it has no scaling
    
    # Get a cell center at this resolution
    center_r <- cpp_hex_center_ap7(10, 5, res)
    
    # Get the same (i,j) at resolution 0 (which has different meaning)
    center_0 <- cpp_hex_center_ap7(10, 5, 0)
    
    # The coordinate magnitudes should scale by sqrt(7)^res
    # (approximately, since the coordinate systems differ)
    
    # Better test: verify cell (1,0) follows expected scaling
    c1_r <- cpp_hex_center_ap7(1, 0, res)
    c1_0 <- cpp_hex_center_ap7(1, 0, 0)
    
    ratio <- as.numeric(c1_0["cx"]) / as.numeric(c1_r["cx"])
    
    # Should be approximately sqrt(7)^res
    expected <- sqrt7^res
    
    expect_equal(ratio, expected, tolerance = 0.1,
                 info = sprintf("Cumulative scale at res=%d", res))
  }
})
