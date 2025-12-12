# tests/testthat/test-projection-inverse.R
# Tests for Snyder ISEA inverse projection
#
# Functions tested:
# - hexify_inverse()
# - hexify_set_precision()
# - hexify_get_precision()
# - hexify_projection_stats()

# =============================================================================
# BASIC INVERSE PROJECTION
# =============================================================================

test_that("inverse projection returns valid lon/lat", {
  hexify_build_icosa()

  result <- hexify_inverse(0.5, 0.3, 0)

  expect_true("lon" %in% names(result))
  expect_true("lat" %in% names(result))
  expect_true(is.finite(result["lon"]))
  expect_true(is.finite(result["lat"]))
})

test_that("inverse projection returns coordinates in valid range", {
  hexify_build_icosa()

  for (face in 0:19) {
    result <- hexify_inverse(0.5, 0.3, face)

    expect_true(result["lon"] >= -180 && result["lon"] <= 180)
    expect_true(result["lat"] >= -90 && result["lat"] <= 90)
  }
})

# =============================================================================
# ROUND-TRIP CONSISTENCY
# =============================================================================

test_that("forward-inverse round-trip works near face centers", {
  hexify_build_icosa()
  centers <- hexify_face_centers()

  for (face in 0:19) {
    # Forward projection
    fwd <- hexify_forward_to_face(face, centers$lon[face + 1], centers$lat[face + 1])

    # Inverse projection
    inv <- hexify_inverse(fwd["icosa_triangle_x"], fwd["icosa_triangle_y"], face)

    expect_true(abs(inv["lat"] - centers$lat[face + 1]) < 1e-6,
                info = sprintf("Face %d lat mismatch", face))
  }
})

test_that("forward-inverse round-trip works for random points", {
  hexify_build_icosa()
  hexify_set_precision("high")

  set.seed(123)

  for (i in 1:50) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -85, 85)  # Avoid extreme poles

    fwd <- hexify_forward(lon, lat)
    face <- as.integer(fwd["face"])

    inv <- hexify_inverse(fwd["icosa_triangle_x"], fwd["icosa_triangle_y"], face)

    lon_diff <- abs(inv["lon"] - lon)
    # Handle longitude wrap-around
    if (lon_diff > 180) lon_diff <- 360 - lon_diff

    expect_true(lon_diff < 1e-5,
                info = sprintf("lon error %.8f at (%.2f, %.2f)", lon_diff, lon, lat))
    expect_true(abs(inv["lat"] - lat) < 1e-5,
                info = sprintf("lat error at (%.2f, %.2f)", lon, lat))
  }
})

# =============================================================================
# PRECISION SETTINGS
# =============================================================================

test_that("precision presets are accepted", {
  expect_no_error(hexify_set_precision("fast"))
  expect_no_error(hexify_set_precision("default"))
  expect_no_error(hexify_set_precision("high"))
  expect_no_error(hexify_set_precision("ultra"))
})

test_that("get_precision returns valid values", {
  hexify_set_precision("fast")
  p <- hexify_get_precision()

  expect_true("tol" %in% names(p))
  expect_true("max_iters" %in% names(p))
  expect_true(is.numeric(p["tol"]))
  expect_true(is.numeric(p["max_iters"]))
})

test_that("precision settings affect iteration count", {
  hexify_build_icosa()

  hexify_set_precision("fast")
  fast_precision <- hexify_get_precision()

  hexify_set_precision("ultra")
  ultra_precision <- hexify_get_precision()

  # Ultra should have tighter tolerance or more iterations
  expect_true(ultra_precision["tol"] <= fast_precision["tol"] ||
                ultra_precision["max_iters"] >= fast_precision["max_iters"])
})

# =============================================================================
# PROJECTION STATS
# =============================================================================

test_that("projection_stats returns valid structure", {
  hexify_build_icosa()

  # Perform some projections
  for (i in 1:10) {
    hexify_inverse(0.5, 0.3, 0)
  }

  stats <- hexify_projection_stats()

  expect_true("calls" %in% names(stats))
  expect_true("iters_total" %in% names(stats))
  expect_true("iters_max" %in% names(stats))
})

test_that("projection_stats tracks calls", {
  hexify_build_icosa()

  # Reset stats
  hexify_projection_stats()

  # Perform known number of calls
  n_calls <- 5
  for (i in 1:n_calls) {
    hexify_inverse(0.5, 0.3, 0)
  }

  stats <- hexify_projection_stats()
  expect_equal(as.integer(stats["calls"]), n_calls)
})
