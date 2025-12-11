# tests/testthat/test_snyder_forward_properties.R

test_that("Snyder forward: finite, in-range, and deterministic", {
  hexify_build_icosa()

  set.seed(123)
  n <- 200
  lon <- runif(n, -180, 180)
  lat <- runif(n, -89.999999, 89.999999)

  # 1) Finite + face in bounds + tx,ty roughly in [0,1]
  #    (projection normalizes by icosaEdge; allow small numerical slack)
  slack <- 1e-12
  for (i in seq_len(n)) {
    out <- hexify_forward(lon[i], lat[i])
    expect_true(is.finite(out["tx"]))
    expect_true(is.finite(out["ty"]))
    expect_true(is.finite(out["face"]))
    expect_true(0 - slack <= out["tx"] && out["tx"] <= 1 + slack)
    expect_true(0 - slack <= out["ty"] && out["ty"] <= 1 + slack)
    expect_true(out["face"] >= 0 && out["face"] <= 19)
  }

  # 2) Determinism: same inputs -> identical outputs
  idx <- sample.int(n, 50)
  for (i in idx) {
    a <- hexify_forward(lon[i], lat[i])
    b <- hexify_forward(lon[i], lat[i])
    expect_identical(as.numeric(a["face"]), as.numeric(b["face"]))
    expect_identical(as.numeric(a["tx"]),   as.numeric(b["tx"]))
    expect_identical(as.numeric(a["ty"]),   as.numeric(b["ty"]))
  }
})

test_that("Snyder forward: full vs fixed-face agree", {
  hexify_build_icosa()

  set.seed(456)
  n <- 100
  lon <- runif(n, -180, 180)
  lat <- runif(n, -89.999999, 89.999999)

  for (i in seq_len(n)) {
    f <- hexify_forward(lon[i], lat[i])
    face <- as.integer(f["face"])
    # use the C++ helper that projects with a fixed face
    xy_face <- hexify_forward_to_face(face, lon[i], lat[i])
    # They should be identical (same branch and arithmetic)
    expect_equal(as.numeric(f["tx"]), as.numeric(xy_face["tx"]))
    expect_equal(as.numeric(f["ty"]), as.numeric(xy_face["ty"]))
  }
})

test_that("Snyder forward: stable near poles and antimeridian", {
  hexify_build_icosa()

  eps <- 1e-9
  # Poles (lon shouldn't matter)
  lons <- c(-180, -90, 0, 90, 180)
  for (L in lons) {
    north <- hexify_forward(L,  90 - eps)
    south <- hexify_forward(L, -90 + eps)
    # just assert finiteness and in-range (same as earlier)
    expect_true(is.finite(north["tx"]) && is.finite(north["ty"]))
    expect_true(is.finite(south["tx"]) && is.finite(south["ty"]))
    expect_true(0 - 1e-12 <= north["tx"] && north["tx"] <= 1 + 1e-12)
    expect_true(0 - 1e-12 <= north["ty"] && north["ty"] <= 1 + 1e-12)
    expect_true(0 - 1e-12 <= south["tx"] && south["tx"] <= 1 + 1e-12)
    expect_true(0 - 1e-12 <= south["ty"] && south["ty"] <= 1 + 1e-12)
  }

  # Antimeridian wrap-around consistency
  # Points just left/right of +180/-180 should project close if on same face
  pts <- data.frame(
    lon = c(180 - 1e-7, -180 + 1e-7, 180 - 1e-6, -180 + 1e-6),
    lat = c(0, 0, 45, -45)
  )
  for (i in seq_len(nrow(pts))) {
    a <- hexify_forward(pts$lon[i], pts$lat[i])
    # Wrap counterpart
    b <- hexify_forward(if (pts$lon[i] > 0) pts$lon[i] - 360 else pts$lon[i] + 360, pts$lat[i])
    # If faces match, outputs should be nearly identical
    if (as.integer(a["face"]) == as.integer(b["face"])) {
      expect_lt(abs(as.numeric(a["tx"]) - as.numeric(b["tx"])), 1e-10)
      expect_lt(abs(as.numeric(a["ty"]) - as.numeric(b["ty"])), 1e-10)
    }
  }
})
