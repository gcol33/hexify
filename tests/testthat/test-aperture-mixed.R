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
