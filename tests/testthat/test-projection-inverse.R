
# tests/testthat/test-projection-inverse.R
# Tests for Snyder ISEA inverse projection
#
# Functions tested:
# - hexify_inverse()

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
# CLOSED FORM
# =============================================================================

# Great-circle angle between points given in degrees, in radians
arc_between <- function(lon1, lat1, lon2, lat2) {
  d2r <- pi / 180
  u <- cbind(cos(lat1 * d2r) * cos(lon1 * d2r), cos(lat1 * d2r) * sin(lon1 * d2r),
             sin(lat1 * d2r))
  v <- cbind(cos(lat2 * d2r) * cos(lon2 * d2r), cos(lat2 * d2r) * sin(lon2 * d2r),
             sin(lat2 * d2r))
  cr <- cbind(u[, 2] * v[, 3] - u[, 3] * v[, 2], u[, 3] * v[, 1] - u[, 1] * v[, 3],
              u[, 1] * v[, 2] - u[, 2] * v[, 1])
  atan2(sqrt(rowSums(cr^2)), rowSums(u * v))
}

# Face-plane points: a barycentric grid over the triangle, edges included,
# and points on the three radii from the centre to the vertices, where
# Snyder's projection has its cusps.
face_plane_samples <- function(n = 30) {
  vx <- c(0, 1, 0.5)
  vy <- c(0, 0, sqrt(3) / 2)
  g <- expand.grid(i = 0:n, j = 0:n)
  g <- g[g$i + g$j <= n, ]
  b1 <- g$i / n
  b2 <- g$j / n
  b0 <- 1 - b1 - b2
  t <- seq(0.01, 1, length.out = 40)
  cx <- mean(vx)
  cy <- mean(vy)
  data.frame(
    x = c(b0 * vx[1] + b1 * vx[2] + b2 * vx[3], cx + outer(t, vx - cx)),
    y = c(b0 * vy[1] + b1 * vy[2] + b2 * vy[3], cy + outer(t, vy - cy))
  )
}

test_that("closed-form inverse agrees with Newton's method on every solid", {
  skip_on_cran()
  pts <- face_plane_samples()
  for (solid in c("icosahedron", "octahedron", "tetrahedron")) {
    icosa <- projection_icosa("isea", solid)
    n_faces <- nrow(hexify_face_centers(solid))
    for (face in unique(c(0L, n_faces %/% 2L, n_faces - 1L))) {
      closed <- t(mapply(function(x, y) cpp_face_xy_to_ll(icosa, x, y, face),
                         pts$x, pts$y))
      newton <- t(mapply(function(x, y) cpp_face_xy_to_ll(icosa, x, y, face, newton = TRUE),
                         pts$x, pts$y))
      err <- arc_between(closed[, "lon"], closed[, "lat"],
                         newton[, "lon"], newton[, "lat"])
      expect_lt(max(err), 1e-13, label = sprintf("%s face %d", solid, face))
    }
  }
})

test_that("forward then closed-form inverse returns the point at machine precision", {
  skip_on_cran()
  set.seed(84)
  lon <- runif(400, -180, 180)
  lat <- asin(runif(400, -1, 1)) * 180 / pi
  for (solid in c("icosahedron", "octahedron", "tetrahedron")) {
    back <- t(mapply(function(lo, la) {
      f <- hexify_forward(lo, la, polyhedron = solid)
      hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]],
                     as.integer(f[["face"]]), polyhedron = solid)
    }, lon, lat))
    err <- arc_between(lon, lat, back[, "lon"], back[, "lat"])
    expect_lt(max(err), 1e-13, label = solid)
  }
})

unit_rows <- function(m) m / sqrt(rowSums(m^2))

test_that("Snyder's map on a hemisphere split in four is Collignon's", {
  # Recht (2021): with v0 the pole and v1, v2 on the equator a quarter turn
  # apart, x = b1 + b2 = sqrt(2) sin(pi/4 - phi/2) and y = b2 - b1 = (4/pi) x lambda.
  g <- expand.grid(lam = seq(-pi / 4, pi / 4, length.out = 41),
                   phi = seq(0.001, pi / 2, length.out = 41))
  v <- cbind(cos(g$phi) * cos(g$lam), cos(g$phi) * sin(g$lam), sin(g$phi))
  b <- cpp_snyder_triangle_forward(c(0, 0, 1), c(1, -1, 0) / sqrt(2),
                                   c(1, 1, 0) / sqrt(2), v)
  x <- sqrt(2) * sin(pi / 4 - g$phi / 2)
  expect_lt(max(abs(b[, 1] + b[, 2] - x)), 1e-14)
  expect_lt(max(abs(b[, 2] - b[, 1] - 4 / pi * x * g$lam)), 1e-14)
})

test_that("Snyder's map on a cube face is the COBE sky-cube formula", {
  # Recht (2021): v0 = (0, 0, 1), v1 = (1, -1, 1)/sqrt(3), v2 = (1, 1, 1)/sqrt(3),
  # b = sqrt(2 vx^2 + vy^2), x = sqrt(b (b + vx) / (1 + vz)),
  # y = x (12/pi) atan(vy / (b + 2 vx)).
  g <- expand.grid(X = seq(0.01, 1, length.out = 40), s = seq(-1, 1, length.out = 41))
  v <- unit_rows(cbind(g$X, g$s * g$X, 1))
  b <- cpp_snyder_triangle_forward(c(0, 0, 1), c(1, -1, 1) / sqrt(3),
                                   c(1, 1, 1) / sqrt(3), v)
  bb <- sqrt(2 * v[, 1]^2 + v[, 2]^2)
  x <- sqrt(bb * (bb + v[, 1]) / (1 + v[, 3]))
  y <- x * 12 / pi * atan(v[, 2] / (bb + 2 * v[, 1]))
  expect_lt(max(abs(b[, 1] + b[, 2] - x)), 1e-14)
  expect_lt(max(abs(b[, 2] - b[, 1] - y)), 1e-14)
})

test_that("Snyder's map inverts on irregular triangles of either orientation", {
  set.seed(96)
  for (k in 1:50) {
    centre <- rnorm(3)
    centre <- centre / sqrt(sum(centre^2))
    tri <- unit_rows(matrix(centre, 3, 3, byrow = TRUE) + matrix(runif(9, -0.6, 0.6), 3))
    b1 <- runif(60)
    b2 <- runif(60) * (1 - b1)
    v <- cpp_snyder_triangle_inverse(tri[1, ], tri[2, ], tri[3, ], cbind(b1, b2))
    b <- cpp_snyder_triangle_forward(tri[1, ], tri[2, ], tri[3, ], v)
    expect_lt(max(abs(b - cbind(b1, b2))), 1e-11)
  }
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
    params <- cpp_icosa_face_params(numeric(0), face)

    expect_true("cen_lat" %in% names(params))
    expect_true("cen_lon" %in% names(params))
    expect_true("face_azimuth_offset" %in% names(params))

    expect_true(params["cen_lat"] >= -90 && params["cen_lat"] <= 90)
    expect_true(params["cen_lon"] >= -180 && params["cen_lon"] <= 180)
  }
})

test_that("cpp_icosa_face_params errors on invalid face", {
  hexify_build_icosa()

  expect_error(cpp_icosa_face_params(numeric(0), -1), "face out of range")
  expect_error(cpp_icosa_face_params(numeric(0), 20), "face out of range")
})

test_that("cpp_hex_index_face_to_lonlat works with degrees=TRUE", {
  hexify_build_icosa()

  # Get face 0 parameters
  params <- cpp_icosa_face_params(numeric(0), 0)

  result <- cpp_hex_index_face_to_lonlat(numeric(0), 
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
  params <- cpp_icosa_face_params(numeric(0), 0)

  result <- cpp_hex_index_face_to_lonlat(numeric(0), 
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
