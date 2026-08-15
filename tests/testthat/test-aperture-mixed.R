
# tests/testthat/test-aperture-mixed.R
# Tests for mixed aperture sequences
#
# An aperture sequence gives the aperture of every resolution level. ISEA43H is
# aperture 4 for the first few levels then aperture 3, which is DGGRID's
# standard configuration; apertures 3, 4 and 7 may appear in any order.

# =============================================================================
# SETUP
# =============================================================================

setup_icosa <- function() {
  cpp_build_icosa()
}

# The ISEA43H sequence: `level` aperture-4 refinements then aperture 3.
seq43 <- function(res, level) {
  c(if (level > 0) 4 else 3, rep(4, level), rep(3, res - level))
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

      cell_34 <- cpp_hex_quantize_mixed(tx, ty, ap_seq)
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

      cell_34 <- cpp_hex_quantize_mixed(tx, ty, ap_seq)
      cell_3 <- cpp_hex_quantize_ap3(tx, ty, res)

      expect_equal(cell_34["i"], cell_3["i"],
                   info = sprintf("res=%d, pure ap3 i mismatch", res))
      expect_equal(cell_34["j"], cell_3["j"],
                   info = sprintf("res=%d, pure ap3 j mismatch", res))

      # Compare centers
      center_34 <- cpp_hex_center_mixed(cell_34["i"], cell_34["j"], ap_seq)
      center_3 <- cpp_hex_center_ap3(cell_3["i"], cell_3["j"], res)

      expect_equal(as.numeric(center_34["cx"]), as.numeric(center_3["cx"]),
                   tolerance = 1e-10)
      expect_equal(as.numeric(center_34["cy"]), as.numeric(center_3["cy"]),
                   tolerance = 1e-10)
    }
  }
})

test_that("pure aperture 7 sequence matches ap7 implementation", {
  skip_on_cran()  # Detailed consistency check
  setup_icosa()

  set.seed(77)
  pts <- cbind(runif(40, -0.8, 0.8), runif(40, -0.8, 0.8))

  for (res in 0:5) {
    ap_seq <- rep(7, res + 1)

    for (k in seq_len(nrow(pts))) {
      tx <- pts[k, 1]
      ty <- pts[k, 2]

      cell_mixed <- cpp_hex_quantize_mixed(tx, ty, ap_seq)
      cell_7 <- cpp_hex_quantize_ap7(tx, ty, res)

      expect_equal(cell_mixed["i"], cell_7["i"],
                   info = sprintf("res=%d, pure ap7 i mismatch", res))
      expect_equal(cell_mixed["j"], cell_7["j"],
                   info = sprintf("res=%d, pure ap7 j mismatch", res))

      center_mixed <- cpp_hex_center_mixed(cell_mixed["i"], cell_mixed["j"], ap_seq)
      center_7 <- cpp_hex_center_ap7(cell_7["i"], cell_7["j"], res)

      expect_equal(as.numeric(center_mixed["cx"]), as.numeric(center_7["cx"]),
                   tolerance = 1e-10)
      expect_equal(as.numeric(center_mixed["cy"]), as.numeric(center_7["cy"]),
                   tolerance = 1e-10)
    }
  }
})

test_that("aperture sequences reject unsupported apertures", {
  setup_icosa()

  expect_error(cpp_hex_quantize_mixed(0.1, 0.2, c(3, 5)), "aperture must be 3, 4, or 7")
  expect_error(cpp_hex_quantize_mixed(0.1, 0.2, c(5, 3)), "aperture must be 3, 4, or 7")
  expect_error(cpp_hex_quantize_mixed(0.1, 0.2, integer(0)), "cannot be empty")
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

      cell <- cpp_hex_quantize_mixed(tx, ty, ap_seq)
      center <- cpp_hex_center_mixed(cell["i"], cell["j"], ap_seq)
      cell2 <- cpp_hex_quantize_mixed(center["cx"], center["cy"], ap_seq)

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

      cell <- cpp_hex_quantize_mixed(tx, ty, seq)
      center <- cpp_hex_center_mixed(cell["i"], cell["j"], seq)
      cell2 <- cpp_hex_quantize_mixed(center["cx"], center["cy"], seq)

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
    c(4, 4, 3, 3),
    c(4, 7),
    c(7, 4),
    c(4, 7, 3, 7)
  )

  for (seq in sequences) {
    corners <- cpp_hex_corners_mixed(0, 0, seq, 1.0)

    expect_equal(length(corners$x), 6)
    expect_equal(length(corners$y), 6)
    expect_true(all(is.finite(corners$x)))
    expect_true(all(is.finite(corners$y)))
  }
})

# =============================================================================
# APERTURE-7 SEQUENCES (issue #55)
# =============================================================================
#
# Aperture 7 rotates the lattice by kAp7RotDeg per level rather than by a
# multiple of 30 degrees, so it composes with apertures 3 and 4 through the
# Eisenstein-generator model in grid_math.h rather than a Class I/II flag.

test_that("sequences mixing aperture 7 with 3 and 4 round-trip", {
  setup_icosa()

  sequences <- list(
    c(4, 7),
    c(7, 4),
    c(4, 4, 7),
    c(7, 7, 4),
    c(3, 7),
    c(7, 3),
    c(4, 7, 3),
    c(7, 4, 3, 7),
    c(4, 4, 7, 7, 3)
  )

  set.seed(755)
  pts <- cbind(runif(25, -0.8, 0.8), runif(25, -0.8, 0.8))

  for (seq in sequences) {
    seq_str <- paste(seq, collapse = ",")
    for (k in seq_len(nrow(pts))) {
      cell <- cpp_hex_quantize_mixed(pts[k, 1], pts[k, 2], seq)
      center <- cpp_hex_center_mixed(cell["i"], cell["j"], seq)
      cell2 <- cpp_hex_quantize_mixed(center["cx"], center["cy"], seq)

      expect_equal(cell["i"], cell2["i"], info = sprintf("seq=[%s] i", seq_str))
      expect_equal(cell["j"], cell2["j"], info = sprintf("seq=[%s] j", seq_str))
    }
  }
})

test_that("every aperture sequence is centre-nested in its refinements", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  # A cell's centre must also be a centre in the next finer grid, whatever
  # aperture that level uses. This is what breaks when the orientation of a
  # level is computed independently of the levels above it.
  parents <- list(
    c(3), c(4), c(7),
    c(3, 3), c(4, 4), c(7, 7),
    c(4, 3), c(3, 4), c(4, 7), c(7, 4), c(3, 7), c(7, 3),
    c(4, 4, 3), c(4, 4, 7), c(7, 7, 4), c(4, 7, 3),
    c(4, 4, 3, 3), c(7, 4, 3, 7)
  )

  set.seed(756)
  pts <- cbind(runif(25, -0.7, 0.7), runif(25, -0.7, 0.7))

  for (parent in parents) {
    for (nxt in c(3, 4, 7)) {
      child <- c(parent, nxt)
      child_spacing <- 1 / prod(sqrt(child[-1]))
      seq_str <- paste(child, collapse = ",")

      for (k in seq_len(nrow(pts))) {
        cell <- cpp_hex_quantize_mixed(pts[k, 1], pts[k, 2], parent)
        p <- cpp_hex_center_mixed(cell["i"], cell["j"], parent)
        kid <- cpp_hex_quantize_mixed(p["cx"], p["cy"], child)
        q <- cpp_hex_center_mixed(kid["i"], kid["j"], child)

        d <- sqrt((p["cx"] - q["cx"])^2 + (p["cy"] - q["cy"])^2)
        expect_lt(as.numeric(d), child_spacing * 1e-6,
                  label = sprintf("seq=[%s] point %d centre offset", seq_str, k))
      }
    }
  }
})

test_that("aperture order does not change the grid", {
  setup_icosa()

  # Scale is a product and orientation composes commutatively, so sequences
  # that use the same multiset of apertures describe the same grid. Each pair
  # swaps two adjacent steps, which leaves the accumulated scale bit-identical.
  set.seed(757)
  pts <- cbind(runif(25, -0.8, 0.8), runif(25, -0.8, 0.8))

  pairs <- list(
    list(c(4, 4, 3, 7), c(4, 3, 4, 7)),
    list(c(3, 3, 4, 7), c(3, 4, 3, 7)),
    list(c(7, 7, 4, 3), c(7, 4, 7, 3))
  )

  for (pr in pairs) {
    for (k in seq_len(nrow(pts))) {
      a <- cpp_hex_quantize_mixed(pts[k, 1], pts[k, 2], pr[[1]])
      b <- cpp_hex_quantize_mixed(pts[k, 1], pts[k, 2], pr[[2]])
      expect_equal(a["i"], b["i"])
      expect_equal(a["j"], b["j"])
    }
  }
})

# =============================================================================
# TRUE MIXED-RADIX SUBSTRATE (i,j) -- not a pure aperture-3 approximation
# =============================================================================
#
# The mixed cell-ID path used to quantize (i,j) as pure aperture-3 at the
# full target resolution, ignoring mixed_aperture_level entirely (only the
# cell-count/offset arithmetic accounted for the aperture-4 prefix). That
# meant distinct real cells could collide onto the same cell ID once the
# aperture-4 prefix meaningfully changed the coordinate scale. These tests
# lock in the fix: (i,j) now come from the same 2^level * sqrt(3)^(res-level)
# substrate that calc_grid_params_mixed()'s cell-count formula describes.

test_that("mixed aperture (i,j) does not collide across many sampled points", {
  skip_on_cran()  # Detailed loop test
  setup_icosa()

  set.seed(43)
  lons <- runif(2000, -170, 170)
  lats <- runif(2000, -80, 80)

  for (mixed_level in c(1, 2, 3)) {
    for (res in c(3, 4, 5, 6)) {
      if (mixed_level >= res) next
      cells <- cpp_lonlat_to_cell_seq(lons, lats, seq43(res, mixed_level))
      centers <- cpp_cell_to_lonlat_seq(cells, seq43(res, mixed_level))
      cells2 <- cpp_lonlat_to_cell_seq(centers$lon_deg, centers$lat_deg, seq43(res, mixed_level))

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
  cells <- cpp_lonlat_to_cell_seq(lons, lats, seq43(res, mixed_level))
  qij <- cpp_cell_to_quad_ij_seq(cells, seq43(res, mixed_level))

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
    cell_mixed <- cpp_lonlat_to_cell_seq(test_lon, test_lat, seq43(res, 0))

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
    cell_mixed <- cpp_lonlat_to_cell_seq(test_lon, test_lat, seq43(res, res))

    # Verify cells are valid (positive integers)
    expect_true(all(cell_mixed >= 1),
                info = sprintf("res=%d, mixed_level=res should produce valid cells", res))
    expect_true(all(is.finite(cell_mixed)),
                info = sprintf("res=%d cells should be finite", res))

    # Verify round-trip works
    centers <- cpp_cell_to_lonlat_seq(cell_mixed, seq43(res, res))
    cell2 <- cpp_lonlat_to_cell_seq(centers$lon_deg, centers$lat_deg, seq43(res, res))
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
      cell <- cpp_lonlat_to_cell_seq(test_lon, test_lat, seq43(res, mixed_level))
      centers <- cpp_cell_to_lonlat_seq(cell, seq43(res, mixed_level))
      cell2 <- cpp_lonlat_to_cell_seq(centers$lon_deg, centers$lat_deg, seq43(res, mixed_level))

      expect_equal(cell, cell2,
                   info = sprintf("res=%d, mixed_level=%d round-trip", res, mixed_level))
    }
  }
})

test_that("mixed aperture rejects an unsupported aperture", {
  setup_icosa()

  expect_error(
    cpp_lonlat_to_cell_seq(0, 0, c(4, 4, 5)),
    "aperture must be 3, 4, or 7"
  )

  expect_error(
    cpp_lonlat_to_cell_seq(0, 0, integer(0)),
    "ap_seq must name at least the base grid"
  )
})

# =============================================================================
# MIXED APERTURE LON/LAT WORKFLOW
# =============================================================================

test_that("cpp_lonlat_to_cell_mixed and cpp_cell_to_lonlat_mixed work", {
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
    cell <- cpp_lonlat_to_cell_mixed(lon, lat, ap_seq)

    expect_true(cell["face"] >= 0 && cell["face"] < 20)
    expect_true(is.numeric(cell["i"]))
    expect_true(is.numeric(cell["j"]))

    ll <- cpp_cell_to_lonlat_mixed(cell["face"], cell["i"], cell["j"], ap_seq)

    expect_true(ll["lon"] >= -180 && ll["lon"] <= 180)
    expect_true(ll["lat"] >= -90 && ll["lat"] <= 90)
  }
})

test_that("cpp_test_roundtrip_mixed returns TRUE for valid points", {
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
      result <- cpp_test_roundtrip_mixed(pt[1], pt[2], ap_seq)
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
  decode <- hexify:::mixed_index_to_cell_one
  for (res in c(2, 4, 6, 7)) {
    grid <- hex_grid(resolution = res, aperture = "4/3")
    set.seed(100 + res)
    # Include high-latitude (near-pole) points so the icosahedron-vertex scatter
    # case is exercised, not just interior cells.
    lon <- c(runif(80, -180, 180), runif(20, -180, 180))
    lat <- c(runif(80, -85, 85), runif(20, 80, 89.9) * sample(c(-1, 1), 20, TRUE))
    cells <- unique(lonlat_to_cell(lon, lat, grid))
    idx <- cell_to_index(cells, grid)
    back <- vapply(idx, function(s) decode(s, res, "4/3"), numeric(1))
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
    ncc <- hexify:::aperture_n_cells("4/3", res + 1L)
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

# =============================================================================
# R-LEVEL SPELLINGS (issue #57)
# =============================================================================
#
# hex_grid() takes a family name ("4/3", "4/7", "7/4") or one aperture per
# resolution level (c(4, 4, 7, 3)). Sequences containing an odd number of
# aperture-7 levels sit on a norm-7 or norm-21 substrate, which packs cell IDs
# on the same sublattice numbering as the Class II lattice of aperture 3.

test_that("hex_grid accepts family spellings and per-level sequences", {
  setup_icosa()

  g <- hex_grid(resolution = 4, aperture = "4/7")
  expect_equal(g@aperture, "4/7")
  expect_equal(hexify:::grid_ap_seq(g), c(4L, 4L, 4L, 7L, 7L))

  g2 <- hex_grid(resolution = 4, aperture = c(4, 7, 3, 7))
  expect_equal(g2@aperture, "4/7/3/7")
  expect_equal(hexify:::grid_ap_seq(g2), c(4L, 4L, 7L, 3L, 7L))

  # Cell count is the product over the sequence, plus the two poles
  expect_equal(hexify:::aperture_n_cells("4/7", 4), 10 * 4 * 4 * 7 * 7 + 2)
  expect_equal(hexify:::aperture_n_cells("4/7/3/7", 4), 10 * 4 * 7 * 3 * 7 + 2)

  expect_error(hex_grid(resolution = 4, aperture = c(4, 5, 3, 7)),
               "must be one of")
  expect_error(hex_grid(resolution = 4, aperture = c(4, 7, 3)),
               "one aperture per resolution level")
})

test_that("mixed sequence cell IDs are a bijection with (quad, i, j)", {
  skip_on_cran()  # Enumerates every cell of the grid
  setup_icosa()

  for (ap in c("4/3", "4/7", "7/4", "3/7", "7/3")) {
    for (res in 1:3) {
      g <- hex_grid(resolution = res, aperture = ap)
      ap_seq <- hexify:::grid_ap_seq(g)
      ids <- seq_len(hexify:::aperture_n_cells(ap, res))

      qij <- cpp_cell_to_quad_ij_seq(as.numeric(ids), ap_seq)
      back <- cpp_quad_ij_to_cell_seq(qij$quad, qij$i, qij$j, ap_seq)

      expect_equal(back, as.numeric(ids),
                   info = sprintf("%s res %d: cell ID round-trip", ap, res))
      expect_equal(length(unique(paste(qij$quad, qij$i, qij$j))), length(ids),
                   info = sprintf("%s res %d: distinct (quad, i, j) per cell", ap, res))

      dim <- cpp_ap_seq_edge_dim(ap_seq)
      expect_true(all(qij$i >= 0 & qij$i < dim & qij$j >= 0 & qij$j < dim),
                  info = sprintf("%s res %d: (i, j) within the quad", ap, res))
    }
  }
})

test_that("mixed sequence points land in a cell centred near them", {
  skip_on_cran()  # Sampled over the sphere
  setup_icosa()

  gc_km <- function(lon1, lat1, lon2, lat2) {
    r <- pi / 180
    6371 * acos(pmin(1, pmax(-1, sin(lat1 * r) * sin(lat2 * r) +
                               cos(lat1 * r) * cos(lat2 * r) * cos((lon2 - lon1) * r))))
  }

  set.seed(57)
  n <- 200
  lon <- 360 * runif(n) - 180
  lat <- asin(2 * runif(n) - 1) * 180 / pi

  for (ap in c("4/7", "7/4", "3/7")) {
    for (res in c(3, 4)) {
      g <- hex_grid(resolution = res, aperture = ap)
      cells <- lonlat_to_cell(lon, lat, g)
      ll <- cell_to_lonlat(cells, g)

      expect_true(all(cells >= 1 & cells <= hexify:::aperture_n_cells(ap, res)),
                  info = sprintf("%s res %d: cell IDs in range", ap, res))

      # Quantizing a cell centre returns that cell
      expect_equal(lonlat_to_cell(ll$lon_deg, ll$lat_deg, g), cells,
                   info = sprintf("%s res %d: centres quantize to their own cell", ap, res))

      # Cell centres sit within a cell spacing of the point
      spacing_km <- sqrt(2 * (EARTH_SURFACE_KM2 / hexify:::aperture_n_cells(ap, res)) / sqrt(3))
      d <- gc_km(lon, lat, ll$lon_deg, ll$lat_deg)
      expect_lt(median(d), spacing_km)
    }
  }
})

test_that("mixed sequence hierarchy navigates for aperture-7 spellings", {
  skip_on_cran()  # Walks children of sampled cells
  setup_icosa()

  for (ap in c("4/7", "7/4")) {
    grid <- hex_grid(resolution = 3, aperture = ap)
    cgrid <- hex_grid(resolution = 4, aperture = ap)
    set.seed(58)
    cells <- unique(lonlat_to_cell(runif(6, -180, 180), runif(6, -80, 80), grid))

    kids <- get_children(cells, grid)
    for (i in seq_along(cells)) {
      expect_true(length(kids[[i]]) >= 1,
                  info = sprintf("%s: cell %.0f has children", ap, cells[i]))
      expect_true(all(get_parent(kids[[i]], cgrid) == cells[i]),
                  info = sprintf("%s: every child maps back to its parent", ap))
    }
  }
})
