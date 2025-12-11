# tests/testthat/test-dggrid-ap4.R
# Test aperture 4 hexagon grid implementation

test_that("Aperture 4 round-trip works", {
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
      cell <- cpp_hex_quantify_ap4(tx, ty, res)
      
      # Get center
      center <- cpp_hex_center_ap4(cell["i"], cell["j"], res)
      
      # Re-quantize center
      cell2 <- cpp_hex_quantify_ap4(center["cx"], center["cy"], res)
      
      # Should match
      expect_equal(cell["i"], cell2["i"], 
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("res=%d, tx=%.3f, ty=%.3f", res, tx, ty))
    }
  }
})

test_that("Aperture 4 always uses Class I (no rotation)", {
  cpp_build_icosa()
  
  # Unlike aperture 3, aperture 4 should have the same orientation at all resolutions
  # All cells should use the Class I formula: x = i - 0.5*j
  
  for (res in c(0, 1, 2, 3, 4)) {
    # Test cell (0,0) - should always be at origin after scaling
    center <- cpp_hex_center_ap4(0, 0, res)
    expect_equal(as.numeric(center["cx"]), 0, tolerance = 1e-10)
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)
    
    # Test cell (1,0) - should be at (1/scale, 0)
    center <- cpp_hex_center_ap4(1, 0, res)
    expected_x <- 1.0 / (2^res)
    expect_equal(as.numeric(center["cx"]), expected_x, tolerance = 1e-10,
                 info = sprintf("res=%d: cell (1,0) x-coord", res))
    expect_equal(as.numeric(center["cy"]), 0, tolerance = 1e-10)
    
    # Test cell (0,1) - should use Class I formula: x = -0.5/scale
    center <- cpp_hex_center_ap4(0, 1, res)
    expected_x <- -0.5 / (2^res)
    expected_y <- sin(60 * pi/180) / (2^res)
    expect_equal(as.numeric(center["cx"]), expected_x, tolerance = 1e-10,
                 info = sprintf("res=%d: cell (0,1) x-coord", res))
    expect_equal(as.numeric(center["cy"]), expected_y, tolerance = 1e-10,
                 info = sprintf("res=%d: cell (0,1) y-coord", res))
  }
})

test_that("Aperture 4 refines by factor of 4", {
  cpp_build_icosa()
  
  # At each resolution, cells should be 1/4 the size
  # Scale should double (2^res)
  
  for (res in 0:4) {
    scale <- 2^res
    
    # Cell (1,0) x-coordinate should be 1/scale
    center <- cpp_hex_center_ap4(1, 0, res)
    expect_equal(as.numeric(center["cx"]), 1.0 / scale, tolerance = 1e-10,
                 info = sprintf("res=%d scale check", res))
  }
  
  # Cell count should increase by 4× at each resolution
  # This is implicit in the quantify logic
})

test_that("Aperture 4 refinement decreases distance", {
  cpp_build_icosa()
  
  # Test point that doesn't sit on aligned centers
  tx <- 0.4
  ty <- 0.35
  
  dists <- numeric()
  
  for (res in c(0, 1, 2, 3, 4, 5)) {
    cell <- cpp_hex_quantify_ap4(tx, ty, res)
    center <- cpp_hex_center_ap4(cell["i"], cell["j"], res)
    dist <- sqrt((center["cx"] - tx)^2 + (center["cy"] - ty)^2)
    dists <- c(dists, dist)
  }
  
  # Distances should generally decrease (may not be strictly monotonic due to geometry)
  # But later resolutions should definitely be closer
  expect_true(dists[5] < dists[3],
              info = sprintf("dist[res4]=%.6f should be < dist[res2]=%.6f", 
                           dists[5], dists[3]))
  expect_true(dists[6] < dists[4],
              info = sprintf("dist[res5]=%.6f should be < dist[res3]=%.6f", 
                           dists[6], dists[4]))
})

test_that("Aperture 4 batch test works", {
  cpp_build_icosa()
  
  set.seed(456)
  n <- 50
  tx <- runif(n, -0.8, 0.8)
  ty <- runif(n, -0.8, 0.8)
  
  for (res in c(0, 1, 2, 3, 4)) {
    failures <- 0
    
    for (i in 1:n) {
      cell <- cpp_hex_quantify_ap4(tx[i], ty[i], res)
      center <- cpp_hex_center_ap4(cell["i"], cell["j"], res)
      cell2 <- cpp_hex_quantify_ap4(center["cx"], center["cy"], res)
      
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

test_that("Aperture 4 corners have correct count and shape", {
  cpp_build_icosa()
  
  for (res in c(0, 1, 2, 3)) {
    # Test a few cells
    cells <- list(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
    
    for (cell in cells) {
      i <- cell[1]
      j <- cell[2]
      
      corners <- cpp_hex_corners_ap4(i, j, res, 1.0)
      
      # Should have 6 corners
      expect_equal(length(corners$x), 6,
                   info = sprintf("res=%d, cell=(%d,%d)", res, i, j))
      expect_equal(length(corners$y), 6)
      
      # Center should be inside the hexagon (approximate check)
      center <- cpp_hex_center_ap4(i, j, res)
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

test_that("Aperture 4 lon/lat workflow works", {
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
      cell <- cpp_lonlat_to_cell_ap4(loc$lon, loc$lat, res)
      
      expect_true(cell["face"] >= 0 && cell["face"] < 20,
                  info = sprintf("%s res=%d: invalid face", loc$name, res))
      expect_true(is.numeric(cell["i"]))
      expect_true(is.numeric(cell["j"]))
      
      # Convert back
      ll <- cpp_cell_to_lonlat_ap4(cell["face"], cell["i"], cell["j"], res)
      
      # Should be reasonably close
      # Higher resolutions should be more accurate
      max_error <- if (res >= 7) 1.0 else if (res >= 5) 5.0 else 10.0
      dist <- sqrt((ll["lon"] - loc$lon)^2 + (ll["lat"] - loc$lat)^2)
      
      expect_true(dist < max_error,
                  info = sprintf("%s res=%d: dist=%.3f > max_error=%.1f",
                               loc$name, res, dist, max_error))
    }
  }
})

test_that("Aperture 4 is simpler than aperture 3", {
  cpp_build_icosa()
  
  # Aperture 4 should not have the Class I/II alternation
  # All resolutions should behave the same way (modulo scaling)
  
  # Test that the same (i,j) at different resolutions scales linearly
  i <- 5
  j <- 3
  
  centers <- list()
  for (res in 0:4) {
    center <- cpp_hex_center_ap4(i, j, res)
    centers[[res + 1]] <- center
  }
  
  # Each resolution should be half the previous (scale = 2)
  for (res in 1:4) {
    ratio_x <- as.numeric(centers[[res]]["cx"]) / as.numeric(centers[[res + 1]]["cx"])
    ratio_y <- as.numeric(centers[[res]]["cy"]) / as.numeric(centers[[res + 1]]["cy"])
    
    expect_equal(ratio_x, 2.0, tolerance = 1e-10,
                 info = sprintf("res %d to %d: x ratio", res-1, res))
    expect_equal(ratio_y, 2.0, tolerance = 1e-10,
                 info = sprintf("res %d to %d: y ratio", res-1, res))
  }
})
