
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
  skip_on_cran()
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
  skip_on_cran()
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
  skip_on_cran()
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

# =============================================================================
# HEXIFY_SET_VERBOSE
# =============================================================================

test_that("hexify_set_verbose accepts TRUE and FALSE", {
  expect_no_error(hexify_set_verbose(TRUE))
  expect_no_error(hexify_set_verbose(FALSE))
})

# =============================================================================
# CUSTOM PRECISION SETTINGS
# =============================================================================

test_that("hexify_set_precision accepts custom tol", {
  expect_no_error(hexify_set_precision(tol = 1e-10))
})

test_that("hexify_set_precision accepts custom max_iters", {
  expect_no_error(hexify_set_precision(max_iters = 50))
})

test_that("hexify_set_precision accepts both custom parameters", {
  expect_no_error(hexify_set_precision(tol = 1e-12, max_iters = 100))
})

# =============================================================================
# HEXIFY_BUILD_ICOSA
# =============================================================================

test_that("hexify_build_icosa with custom parameters", {
  # Custom vertex position
  expect_no_error(hexify_build_icosa(vert0_lon = 0, vert0_lat = 90, azimuth = 0))

  # Reset to standard orientation
  hexify_build_icosa()
})

# =============================================================================
# HEXIFY_FACE_CENTERS
# =============================================================================

test_that("hexify_face_centers returns 20 faces", {
  hexify_build_icosa()

  centers <- hexify_face_centers()

  expect_s3_class(centers, "data.frame")
  expect_equal(nrow(centers), 20)
  expect_true(all(c("lon", "lat") %in% names(centers)))
})

test_that("hexify_face_centers returns valid coordinates", {
  hexify_build_icosa()

  centers <- hexify_face_centers()

  expect_true(all(centers$lon >= -180 & centers$lon <= 180))
  expect_true(all(centers$lat >= -90 & centers$lat <= 90))
})

# =============================================================================
# HEXIFY_WHICH_FACE
# =============================================================================

test_that("hexify_which_face returns valid face indices", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(42)
  for (i in 1:50) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -89, 89)

    face <- hexify_which_face(lon, lat)
    expect_true(face >= 0 && face <= 19)
  }
})

test_that("hexify_which_face is consistent with hexify_forward", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(123)
  for (i in 1:30) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -85, 85)

    face <- hexify_which_face(lon, lat)
    forward_result <- hexify_forward(lon, lat)

    expect_equal(face, as.integer(forward_result["face"]))
  }
})

# =============================================================================
# HEXIFY_INVERSE WITH CUSTOM PARAMETERS
# =============================================================================

test_that("hexify_inverse with custom tol parameter", {
  hexify_build_icosa()

  result <- hexify_inverse(0.5, 0.3, face = 0, tol = 1e-10)

  expect_true(is.finite(result["lon"]))
  expect_true(is.finite(result["lat"]))
})

test_that("hexify_inverse with custom max_iters parameter", {
  hexify_build_icosa()

  result <- hexify_inverse(0.5, 0.3, face = 0, max_iters = 50)

  expect_true(is.finite(result["lon"]))
  expect_true(is.finite(result["lat"]))
})

test_that("hexify_inverse validates input lengths", {
  hexify_build_icosa()

  expect_error(hexify_inverse(c(0.5, 0.6), 0.3, face = 0))
  expect_error(hexify_inverse(0.5, c(0.3, 0.4), face = 0))
  expect_error(hexify_inverse(0.5, 0.3, face = c(0, 1)))
})

# =============================================================================
# INTERNAL CPP FUNCTION TESTS
# =============================================================================

test_that("cpp_icosa_face_params returns valid face parameters", {
  skip_on_cran()
  hexify_build_icosa()

  for (face in 0:19) {
    params <- cpp_icosa_face_params(face)

    expect_true("cen_lat" %in% names(params))
    expect_true("cen_lon" %in% names(params))
    expect_true("face_azimuth_offset" %in% names(params))

    expect_true(params["cen_lat"] >= -90 && params["cen_lat"] <= 90)
    expect_true(params["cen_lon"] >= -180 && params["cen_lon"] <= 180)
  }
})

test_that("cpp_icosa_face_params errors on invalid face", {
  hexify_build_icosa()

  expect_error(cpp_icosa_face_params(-1), "face must be 0..19")
  expect_error(cpp_icosa_face_params(20), "face must be 0..19")
})

test_that("cpp_hex_index_face_to_lonlat works with degrees=TRUE", {
  hexify_build_icosa()

  # Get face 0 parameters
  params <- cpp_icosa_face_params(0)

  result <- cpp_hex_index_face_to_lonlat(
    x = 0.5,
    y = 0.3,
    cen_lat = params["cen_lat"],
    cen_lon = params["cen_lon"],
    face_azimuth_offset = params["face_azimuth_offset"],
    degrees = TRUE
  )

  expect_length(result, 2)
  expect_true(result[1] >= -180 && result[1] <= 180)  # lon
  expect_true(result[2] >= -90 && result[2] <= 90)    # lat
})

test_that("cpp_hex_index_face_to_lonlat works with degrees=FALSE", {
  hexify_build_icosa()

  # Get face 0 parameters
  params <- cpp_icosa_face_params(0)

  result <- cpp_hex_index_face_to_lonlat(
    x = 0.5,
    y = 0.3,
    cen_lat = params["cen_lat"],
    cen_lon = params["cen_lon"],
    face_azimuth_offset = params["face_azimuth_offset"],
    degrees = FALSE
  )

  expect_length(result, 2)
  # Radians: lon in [-pi, pi], lat in [-pi/2, pi/2]
  expect_true(result[1] >= -pi && result[1] <= pi)
  expect_true(result[2] >= -pi / 2 && result[2] <= pi / 2)
})

test_that("cpp_hex_index_face_to_lonlat works with custom tolerance", {
  hexify_build_icosa()

  params <- cpp_icosa_face_params(5)

  result <- cpp_hex_index_face_to_lonlat(
    x = 0.5,
    y = 0.3,
    cen_lat = params["cen_lat"],
    cen_lon = params["cen_lon"],
    face_azimuth_offset = params["face_azimuth_offset"],
    degrees = TRUE,
    tol = 1e-10,
    max_iters = 50
  )

  expect_length(result, 2)
  expect_true(is.finite(result[1]))
  expect_true(is.finite(result[2]))
})
