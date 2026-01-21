
# tests/testthat/test-projection-forward.R
# Tests for Snyder ISEA forward projection
#
# Functions tested:
# - hexify_forward()
# - hexify_forward_to_face()

# =============================================================================
# BASIC FORWARD PROJECTION
# =============================================================================

test_that("forward projection returns valid structure", {
  hexify_build_icosa()

  result <- hexify_forward(16.37, 48.21)

  expect_true("face" %in% names(result))
  expect_true("icosa_triangle_x" %in% names(result))
  expect_true("icosa_triangle_y" %in% names(result))
})

test_that("forward projection returns finite values", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(123)

  n <- 100
  lon <- runif(n, -180, 180)
  lat <- runif(n, -89.99, 89.99)

  for (i in seq_len(n)) {
    result <- hexify_forward(lon[i], lat[i])

    expect_true(is.finite(result["face"]))
    expect_true(is.finite(result["icosa_triangle_x"]))
    expect_true(is.finite(result["icosa_triangle_y"]))
  }
})

test_that("forward projection face values are in valid range", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(456)
  n <- 100

  for (i in 1:n) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -89.99, 89.99)

    result <- hexify_forward(lon, lat)
    expect_true(result["face"] >= 0 && result["face"] <= 19)
  }
})

test_that("forward projection is deterministic", {
  hexify_build_icosa()

  lon <- 16.37
  lat <- 48.21

  result1 <- hexify_forward(lon, lat)
  result2 <- hexify_forward(lon, lat)

  expect_identical(as.numeric(result1["face"]), as.numeric(result2["face"]))
  expect_identical(result1["icosa_triangle_x"], result2["icosa_triangle_x"])
  expect_identical(result1["icosa_triangle_y"], result2["icosa_triangle_y"])
})

# =============================================================================
# FORWARD TO FACE
# =============================================================================

test_that("forward_to_face matches forward for detected face", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(789)
  n <- 100

  for (i in 1:n) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -89.99, 89.99)

    full <- hexify_forward(lon, lat)
    face <- as.integer(full["face"])

    fixed <- hexify_forward_to_face(face, lon, lat)

    expect_equal(as.numeric(full["icosa_triangle_x"]),
                 as.numeric(fixed["icosa_triangle_x"]),
                 tolerance = 1e-15)
    expect_equal(as.numeric(full["icosa_triangle_y"]),
                 as.numeric(fixed["icosa_triangle_y"]),
                 tolerance = 1e-15)
  }
})

test_that("forward_to_face returns values in [0, 1] range", {
  skip_on_cran()
  hexify_build_icosa()

  set.seed(111)
  slack <- 1e-10

  for (i in 1:100) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -89.99, 89.99)

    result <- hexify_forward(lon, lat)
    tx <- result["icosa_triangle_x"]
    ty <- result["icosa_triangle_y"]

    expect_true(tx >= -slack && tx <= 1 + slack,
                info = sprintf("tx=%.10f out of range", tx))
    expect_true(ty >= -slack && ty <= 1 + slack,
                info = sprintf("ty=%.10f out of range", ty))
  }
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("forward projection handles near-pole coordinates", {
  skip_on_cran()
  hexify_build_icosa()

  eps <- 1e-9
  lons <- c(-180, -90, 0, 90, 180)

  for (L in lons) {
    north <- hexify_forward(L, 90 - eps)
    south <- hexify_forward(L, -90 + eps)

    expect_true(is.finite(north["icosa_triangle_x"]))
    expect_true(is.finite(north["icosa_triangle_y"]))
    expect_true(is.finite(south["icosa_triangle_x"]))
    expect_true(is.finite(south["icosa_triangle_y"]))
  }
})

test_that("forward projection handles antimeridian", {
  hexify_build_icosa()

  # Points near +180/-180
  eps <- 1e-7

  east <- hexify_forward(180 - eps, 0)
  west <- hexify_forward(-180 + eps, 0)

  expect_true(is.finite(east["icosa_triangle_x"]))
  expect_true(is.finite(west["icosa_triangle_x"]))
})

# =============================================================================
# REFERENCE DATA VALIDATION
# =============================================================================

test_that("forward projection matches reference data", {
  skip_on_cran()
  hexify_build_icosa(11.25, 58.28252559, 0)
  hexify_set_precision("ultra")

  truth <- read.csv(testthat::test_path("data/truth_txty.csv"),
                    stringsAsFactors = FALSE)

  expect_true(all(c("lon", "lat", "tnum", "tx", "ty") %in% names(truth)))

  atol <- 2e-8
  mismatches <- 0L
  msgs <- character()

  for (i in seq_len(nrow(truth))) {
    result <- hexify_forward(truth$lon[i], truth$lat[i])
    face <- as.integer(result["face"])
    tx <- as.numeric(result["icosa_triangle_x"])
    ty <- as.numeric(result["icosa_triangle_y"])

    ok <- (face == as.integer(truth$tnum[i])) &&
      (abs(tx - truth$tx[i]) <= atol) &&
      (abs(ty - truth$ty[i]) <= atol)

    if (!ok) {
      mismatches <- mismatches + 1L
      if (mismatches <= 5L) {
        msgs <- c(msgs, sprintf(
          "lon=%.10f lat=%.10f: face %d vs %d, tx %.12g vs %.12g, ty %.12g vs %.12g",
          truth$lon[i], truth$lat[i],
          face, truth$tnum[i], tx, truth$tx[i], ty, truth$ty[i]
        ))
      }
    }
  }

  if (mismatches > 0L) {
    message(paste(msgs, collapse = "\n"))
  }
  expect_equal(mismatches, 0L)
})
