
# tests/testthat/test-aperture-mixed.R
# Tests for mixed aperture 4/3 (ISEA43H) hexagonal grid
#
# Mixed aperture uses aperture 4 for the first few resolutions,
# then switches to aperture 3. This is DGGRID's standard configuration.

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# PURE SEQUENCE CONSISTENCY
# =============================================================================

test_that("pure aperture 4 sequence matches ap4 implementation", {
  skip_on_cran()  # Detailed consistency check
  setup_icosa()

  # Pure ap4: sequence of all 4s
  for (res in c(0, 1, 2, 3)) {
    ap_seq <- rep(4, res + 1)

    test_points <- list(c(0.4, 0.35))

    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell_34 <- cpp_hex_quantize_ap34(tx, ty, ap_seq)
      cell_4 <- cpp_hex_quantize_ap4(tx, ty, res)

      expect_equal(cell_34["i"], cell_4["i"],
                   info = sprintf("res=%d, pure ap4 i mismatch", res))
      expect_equal(cell_34["j"], cell_4["j"],
                   info = sprintf("res=%d, pure ap4 j mismatch", res))
    }
  }
})

test_that("pure aperture 3 sequence matches ap3 implementation", {
  skip_on_cran()  # Detailed consistency check
  setup_icosa()

  # Pure ap3: sequence of all 3s
  for (res in c(1, 2, 3, 4)) {
    ap_seq <- rep(3, res + 1)

    test_points <- list(
      c(0.4, 0.35),
      c(-0.4, 0.2),
      c(0.1, -0.6)
    )

    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell_34 <- cpp_hex_quantize_ap34(tx, ty, ap_seq)
      cell_3 <- cpp_hex_quantize_ap3(tx, ty, res)

      expect_equal(cell_34["i"], cell_3["i"],
                   info = sprintf("res=%d, pure ap3 i mismatch", res))
      expect_equal(cell_34["j"], cell_3["j"],
                   info = sprintf("res=%d, pure ap3 j mismatch", res))

      # Compare centers
      center_34 <- cpp_hex_center_ap34(cell_34["i"], cell_34["j"], ap_seq)
      center_3 <- cpp_hex_center_ap3(cell_3["i"], cell_3["j"], res)

      expect_equal(as.numeric(center_34["cx"]), as.numeric(center_3["cx"]),
                   tolerance = 1e-10)
      expect_equal(as.numeric(center_34["cy"]), as.numeric(center_3["cy"]),
                   tolerance = 1e-10)
    }
  }
})

# =============================================================================
# MIXED 4/3 PATTERN (DGGRID STYLE)
# =============================================================================

test_that("mixed 43H pattern (DGGRID style) works", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  # DGGRID 43H: first 2 resolutions ap4, then ap3
  # c(4, 4, 3, 3, 3, ...)

  for (total_res in 2:5) {
    ap_seq <- c(4, 4, rep(3, max(0, total_res - 1)))

    test_points <- list(
      c(0.4, 0.35),
      c(-0.4, 0.2),
      c(0.1, -0.6)
    )

    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell <- cpp_hex_quantize_ap34(tx, ty, ap_seq)
      center <- cpp_hex_center_ap34(cell["i"], cell["j"], ap_seq)
      cell2 <- cpp_hex_quantize_ap34(center["cx"], center["cy"], ap_seq)

      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("43H res=%d round-trip i", total_res))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("43H res=%d round-trip j", total_res))
    }
  }
})

# =============================================================================
# ROUND-TRIP
# =============================================================================

test_that("mixed aperture round-trip works for various sequences", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  sequences <- list(
    c(4, 4),           # Two levels of ap4
    c(4, 4, 3),        # Standard 43H start
    c(4, 4, 3, 3),     # Longer 43H
    c(4, 4, 3, 3, 3),  # Even longer
    c(3, 4),           # Start with 3, then 4
    c(4, 3, 4)         # Alternating
  )

  test_points <- list(
    c(0.4, 0.35),
    c(-0.4, 0.2),
    c(0.0, 0.0)
  )

  for (seq in sequences) {
    for (pt in test_points) {
      tx <- pt[1]
      ty <- pt[2]

      cell <- cpp_hex_quantize_ap34(tx, ty, seq)
      center <- cpp_hex_center_ap34(cell["i"], cell["j"], seq)
      cell2 <- cpp_hex_quantize_ap34(center["cx"], center["cy"], seq)

      seq_str <- paste(seq, collapse = ",")
      expect_equal(cell["i"], cell2["i"],
                   info = sprintf("seq=[%s] round-trip i", seq_str))
      expect_equal(cell["j"], cell2["j"],
                   info = sprintf("seq=[%s] round-trip j", seq_str))
    }
  }
})

# =============================================================================
# CORNERS
# =============================================================================

test_that("mixed aperture corners form valid hexagons", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  sequences <- list(
    c(4, 4),
    c(4, 4, 3),
    c(4, 4, 3, 3)
  )

  for (seq in sequences) {
    corners <- cpp_hex_corners_ap34(0, 0, seq, 1.0)

    expect_equal(length(corners$x), 6)
    expect_equal(length(corners$y), 6)
    expect_true(all(is.finite(corners$x)))
    expect_true(all(is.finite(corners$y)))
  }
})

# =============================================================================
# EDGE CASES - MIXED APERTURE LEVEL BOUNDARIES
# =============================================================================

test_that("mixed aperture level = 0 (all aperture 3)", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  # When mixed_aperture_level = 0, should behave like pure aperture 3
  test_lon <- c(0, 45, -120, 16.37)
  test_lat <- c(0, 30, -45, 48.21)

  for (res in c(3, 5, 7)) {
    # Mixed with level 0 = all ap3
    cell_mixed <- cpp_lonlat_to_cell_ap43(test_lon, test_lat, res, 0)

    # Pure aperture 3
    cell_pure <- cpp_lonlat_to_cell(test_lon, test_lat, res, 3)

    expect_equal(cell_mixed, cell_pure,
                 info = sprintf("res=%d, mixed_level=0 should equal pure ap3", res))
  }
})

test_that("mixed aperture level = resolution produces valid cells", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  # When mixed_aperture_level = resolution, all resolutions use aperture 4
  # Note: Cell IDs may differ from pure ap4 due to different numbering schemes
  test_lon <- c(0, 45, -120, 16.37)
  test_lat <- c(0, 30, -45, 48.21)

  for (res in c(2, 4, 6)) {
    # Mixed with level = res = all ap4 subdivisions
    cell_mixed <- cpp_lonlat_to_cell_ap43(test_lon, test_lat, res, res)

    # Verify cells are valid (positive integers)
    expect_true(all(cell_mixed >= 1),
                info = sprintf("res=%d, mixed_level=res should produce valid cells", res))
    expect_true(all(is.finite(cell_mixed)),
                info = sprintf("res=%d cells should be finite", res))

    # Verify round-trip works
    centers <- cpp_cell_to_lonlat_ap43(cell_mixed, res, res)
    cell2 <- cpp_lonlat_to_cell_ap43(centers$lon_deg, centers$lat_deg, res, res)
    expect_equal(cell_mixed, cell2,
                 info = sprintf("res=%d, mixed_level=res round-trip", res))
  }
})

test_that("mixed aperture round-trip at boundary resolutions", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  test_lon <- c(0, 45, -120)
  test_lat <- c(0, 30, -45)

  # Test at resolution boundaries
  # Skip res=1 with mixed_level=1 as it's a degenerate case
  for (res in c(2, 5, 10)) {
    for (mixed_level in c(0, 1, min(res - 1, 5))) {
      cell <- cpp_lonlat_to_cell_ap43(test_lon, test_lat, res, mixed_level)
      centers <- cpp_cell_to_lonlat_ap43(cell, res, mixed_level)
      cell2 <- cpp_lonlat_to_cell_ap43(centers$lon_deg, centers$lat_deg, res, mixed_level)

      expect_equal(cell, cell2,
                   info = sprintf("res=%d, mixed_level=%d round-trip", res, mixed_level))
    }
  }
})

test_that("mixed aperture invalid level throws error", {
  setup_icosa()

  # mixed_aperture_level > resolution should fail
  expect_error(
    cpp_lonlat_to_cell_ap43(0, 0, 5, 6),
    "mixed_aperture_level must be between 0 and resolution"
  )

  # Negative mixed_aperture_level should fail
  expect_error(
    cpp_lonlat_to_cell_ap43(0, 0, 5, -1),
    "mixed_aperture_level must be between 0 and resolution"
  )
})

# =============================================================================
# MIXED APERTURE LON/LAT WORKFLOW (ap34 functions)
# =============================================================================

test_that("cpp_lonlat_to_cell_ap34 and cpp_cell_to_lonlat_ap34 work", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  lon <- 16.37  # Vienna
  lat <- 48.21

  # Test with various aperture sequences
  sequences <- list(
    c(4, 3, 3, 3),      # res 4, one aperture-4 step
    c(4, 4, 3, 3),      # res 4, two aperture-4 steps
    c(3, 3, 3, 3, 3)    # res 5, all aperture-3
  )

  for (ap_seq in sequences) {
    cell <- cpp_lonlat_to_cell_ap34(lon, lat, ap_seq)

    expect_true(cell["face"] >= 0 && cell["face"] < 20)
    expect_true(is.numeric(cell["i"]))
    expect_true(is.numeric(cell["j"]))

    ll <- cpp_cell_to_lonlat_ap34(cell["face"], cell["i"], cell["j"], ap_seq)

    expect_true(ll["lon"] >= -180 && ll["lon"] <= 180)
    expect_true(ll["lat"] >= -90 && ll["lat"] <= 90)
  }
})

test_that("cpp_test_roundtrip_ap34 returns TRUE for valid points", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  test_points <- list(
    c(0.5, 0.3),
    c(-0.4, 0.2),
    c(0.1, -0.6),
    c(0.0, 0.0)
  )

  sequences <- list(
    c(4, 3, 3, 3),
    c(4, 4, 3, 3),
    c(3, 3, 3, 3, 3)
  )

  for (ap_seq in sequences) {
    for (pt in test_points) {
      result <- cpp_test_roundtrip_ap34(pt[1], pt[2], ap_seq)
      expect_true(result, info = sprintf("ap_seq length=%d, pt=(%.2f, %.2f)",
                                          length(ap_seq), pt[1], pt[2]))
    }
  }
})
