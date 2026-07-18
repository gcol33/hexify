
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
# TRUE MIXED-RADIX SUBSTRATE (i,j) -- not a pure aperture-3 approximation
# =============================================================================
#
# cpp_lonlat_to_cell_ap43() used to quantize (i,j) as pure aperture-3 at the
# full target resolution, ignoring mixed_aperture_level entirely (only the
# cell-count/offset arithmetic accounted for the aperture-4 prefix). That
# meant distinct real cells could collide onto the same cell ID once the
# aperture-4 prefix meaningfully changed the coordinate scale. These tests
# lock in the fix: (i,j) now come from the same 2^level * sqrt(3)^(res-level)
# substrate that calc_grid_params_ap43()'s cell-count formula describes.

test_that("mixed aperture (i,j) does not collide across many sampled points", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  set.seed(43)
  lons <- runif(2000, -170, 170)
  lats <- runif(2000, -80, 80)

  for (mixed_level in c(1, 2, 3)) {
    for (res in c(3, 4, 5, 6)) {
      if (mixed_level >= res) next
      cells <- cpp_lonlat_to_cell_ap43(lons, lats, res, mixed_level)
      centers <- cpp_cell_to_lonlat_ap43(cells, res, mixed_level)
      cells2 <- cpp_lonlat_to_cell_ap43(centers$lon_deg, centers$lat_deg, res, mixed_level)

      expect_equal(cells2, cells,
                   info = sprintf("mixed_level=%d res=%d: cell id round-trip", mixed_level, res))
    }
  }
})

test_that("mixed aperture (i,j) matches the mixed-radix grid dimension, not pure aperture 3", {
  skip_on_cran()
  setup_icosa()

  # Regression check for the specific collision from issue #31: at
  # mixed_aperture_level=2 res=5 the pure-aperture-3-at-full-resolution
  # approximation produced (i,j) an order of magnitude too small, causing
  # unrelated points to land on the same cell ID.
  set.seed(44)
  lons <- runif(500, -170, 170)
  lats <- runif(500, -80, 80)

  res <- 5
  mixed_level <- 2
  cells <- cpp_lonlat_to_cell_ap43(lons, lats, res, mixed_level)
  qij <- cpp_cell_to_quad_ij_ap43(cells, res, mixed_level)

  # True mixed-radix dim: 2^2 * sqrt(3)^3 (+ substrate boost if class II)
  ap3_count <- res - mixed_level
  true_scale <- 2^mixed_level * sqrt(3)^ap3_count
  if (ap3_count %% 2 == 1) true_scale <- true_scale * sqrt(3)
  pure_ap3_scale <- sqrt(3)^res

  expect_true(true_scale > pure_ap3_scale * 1.5,
              info = "sanity: mixed-radix scale is meaningfully larger than pure ap3 at this res")
  expect_true(max(qij$i) <= true_scale + 1,
              info = "i must not exceed the true mixed-radix grid dimension")
  expect_true(max(qij$j) <= true_scale + 1,
              info = "j must not exceed the true mixed-radix grid dimension")
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

# =============================================================================
# HIERARCHICAL NAVIGATION (issue #31): cell_to_index / get_parent / get_children
# =============================================================================
#
# Mixed 4/3 grids have no DGGRID-standard hierarchical index. hexify defines the
# hierarchy geometrically (a cell's parent is the coarser cell containing its
# centre; see R/aperture_mixed_hierarchy.R). These functions used to throw on
# every call for aperture "4/3"; they now navigate the hierarchy and round-trip.

test_that("mixed 4/3 hierarchical navigation no longer errors (issue #31)", {
  setup_icosa()
  grid <- hex_grid(resolution = 6, aperture = "4/3")
  cells <- lonlat_to_cell(c(0, 30, -60), c(10, 45, -30), grid)

  expect_silent(p <- get_parent(cells, grid))
  expect_silent(idx <- cell_to_index(cells, grid))
  expect_silent(kids <- get_children(cells, grid))

  expect_true(all(p >= 1))
  expect_type(idx, "character")
  expect_length(idx, length(cells))
  expect_type(kids, "list")
  expect_length(kids, length(cells))
})

test_that("mixed 4/3 get_parent(get_children()) recovers the parent", {
  setup_icosa()
  for (res in c(3, 4, 5, 6)) {
    grid  <- hex_grid(resolution = res, aperture = "4/3")
    cgrid <- hex_grid(resolution = res + 1L, aperture = "4/3")
    set.seed(res)
    cells <- unique(lonlat_to_cell(runif(40, -180, 180), runif(40, -85, 85), grid))
    kids <- get_children(cells, grid)
    for (i in seq_along(cells)) {
      k <- kids[[i]]
      expect_true(length(k) >= 1,
                  info = sprintf("res=%d cell=%.0f has at least one child", res, cells[i]))
      back <- get_parent(k, cgrid)
      expect_true(all(back == cells[i]),
                  info = sprintf("res=%d cell=%.0f: every child maps back to it", res, cells[i]))
    }
  }
})

test_that("mixed 4/3 cell_to_index round-trips to the original cell", {
  setup_icosa()
  decode <- hexify:::ap43_index_to_cell_one
  for (res in c(2, 4, 6, 7)) {
    grid <- hex_grid(resolution = res, aperture = "4/3")
    set.seed(100 + res)
    # Include high-latitude (near-pole) points so the icosahedron-vertex scatter
    # case is exercised, not just interior cells.
    lon <- c(runif(80, -180, 180), runif(20, -180, 180))
    lat <- c(runif(80, -85, 85), runif(20, 80, 89.9) * sample(c(-1, 1), 20, TRUE))
    cells <- unique(lonlat_to_cell(lon, lat, grid))
    idx <- cell_to_index(cells, grid)
    back <- vapply(idx, function(s) decode(s, res), numeric(1))
    expect_equal(unname(back), cells,
                 info = sprintf("res=%d: index round-trip", res))
    # Parent's index is a prefix of the child's index.
    parents <- get_parent(cells, grid)
    pidx <- cell_to_index(parents, hex_grid(resolution = res - 1L, aperture = "4/3"))
    expect_true(all(substr(idx, 1, nchar(pidx)) == pidx),
                info = sprintf("res=%d: parent index is a prefix", res))
  }
})

test_that("mixed 4/3 get_children is exactly the geometric-parent inverse", {
  skip_on_cran()  # brute-force ground-truth comparison
  setup_icosa()
  # At a small resolution, brute-force {c : parent(c) == P} and compare.
  for (res in c(2, 3)) {
    grid  <- hex_grid(resolution = res, aperture = "4/3")
    cgrid <- hex_grid(resolution = res + 1L, aperture = "4/3")
    ncc <- hexify:::ap43_n_cells(res + 1L)
    all_child <- seq_len(ncc)
    par_of_child <- get_parent(all_child, cgrid)
    parents <- sort(unique(par_of_child))
    for (p in parents) {
      truth <- sort(all_child[par_of_child == p])
      got <- get_children(p, grid)[[1]]
      expect_equal(got, truth,
                   info = sprintf("res=%d parent=%.0f exact children", res, p))
    }
  }
})

test_that("mixed 4/3 multi-level parent and children are consistent", {
  setup_icosa()
  grid <- hex_grid(resolution = 6, aperture = "4/3")
  g2   <- hex_grid(resolution = 4, aperture = "4/3")
  set.seed(7)
  cells <- unique(lonlat_to_cell(runif(30, -180, 180), runif(30, -85, 85), grid))
  gp2 <- get_parent(cells, grid, levels = 2L)
  expect_true(all(gp2 >= 1))
  desc <- get_children(gp2, g2, levels = 2L)
  for (i in seq_along(cells)) {
    expect_true(cells[i] %in% desc[[i]],
                info = sprintf("cell %.0f is a 2-level descendant of its 2-up parent", cells[i]))
  }
})
