
# tests/testthat/test-quantizer-nearest.R
# Quantization returns the cell whose centre is nearest the point.
#
# The candidate cells are collected by quantizing a scan of offsets covering
# one cell spacing around the point: the nearest cell's centre lies inside that
# disc, and a point at a centre quantizes to its own cell, so every cell that
# could beat the returned one is in the set. This holds for every grid form,
# rotated or not, and needs no separate model of which (i, j) are grid cells.

scan_offsets <- function(spacing, n = 9) {
  s <- seq(-1.2 * spacing, 1.2 * spacing, length.out = n)
  expand.grid(dx = s, dy = s)
}

nearest_cell_holds <- function(quantize, center, px, py, spacing) {
  cell <- quantize(px, py)
  c0 <- center(cell[["i"]], cell[["j"]])
  best <- (c0[["cx"]] - px)^2 + (c0[["cy"]] - py)^2

  off <- scan_offsets(spacing)
  for (k in seq_len(nrow(off))) {
    q <- quantize(px + off$dx[k], py + off$dy[k])
    cc <- center(q[["i"]], q[["j"]])
    d <- (cc[["cx"]] - px)^2 + (cc[["cy"]] - py)^2
    if (d < best - 1e-12) return(FALSE)
  }
  TRUE
}

# Cell spacing of a sequence: one over the product of sqrt(aperture) per step
seq_spacing <- function(ap_seq) 1 / prod(sqrt(ap_seq[-1]))

test_that("pure aperture quantization returns the nearest cell", {
  skip_on_cran()
  cpp_build_icosa()

  set.seed(101)
  pts <- cbind(runif(12, -0.7, 0.7), runif(12, -0.7, 0.7))

  quantizers <- list(
    "3" = list(q = cpp_hex_quantize_ap3, c = cpp_hex_center_ap3),
    "4" = list(q = cpp_hex_quantize_ap4, c = cpp_hex_center_ap4),
    "7" = list(q = cpp_hex_quantize_ap7, c = cpp_hex_center_ap7)
  )

  for (ap in names(quantizers)) {
    fns <- quantizers[[ap]]
    for (res in 0:3) {
      spacing <- 1 / (sqrt(as.numeric(ap))^res)
      for (k in seq_len(nrow(pts))) {
        ok <- nearest_cell_holds(
          function(x, y) fns$q(x, y, res),
          function(i, j) fns$c(i, j, res),
          pts[k, 1], pts[k, 2], spacing
        )
        expect_true(ok, info = sprintf("aperture %s res %d, point (%.3f, %.3f)",
                                       ap, res, pts[k, 1], pts[k, 2]))
      }
    }
  }
})

test_that("mixed sequence quantization returns the nearest cell", {
  skip_on_cran()
  cpp_build_icosa()

  set.seed(202)
  pts <- cbind(runif(8, -0.7, 0.7), runif(8, -0.7, 0.7))

  sequences <- list(
    c(4, 4, 4, 3, 3),
    c(4, 4, 3, 7),
    c(3, 7, 4, 3),
    c(7, 7, 7, 3),
    c(3, 3, 4, 4, 7)
  )

  for (ap_seq in sequences) {
    spacing <- seq_spacing(ap_seq)
    for (k in seq_len(nrow(pts))) {
      ok <- nearest_cell_holds(
        function(x, y) cpp_hex_quantize_mixed(x, y, ap_seq),
        function(i, j) cpp_hex_center_mixed(i, j, ap_seq),
        pts[k, 1], pts[k, 2], spacing
      )
      expect_true(ok, info = sprintf("sequence %s, point (%.3f, %.3f)",
                                     paste(ap_seq, collapse = "/"),
                                     pts[k, 1], pts[k, 2]))
    }
  }
})
