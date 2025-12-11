test_that("hexify_forward_to_face matches hexify_forward face path", {
  hexify_build_icosa(11.25, 58.28252559, 0)

  # Vienna sanity
  v <- hexify_forward(16.3738, 48.2082)
  xy <- hexify_forward_to_face(as.integer(v[["face"]]), 16.3738, 48.2082)
  expect_equal(as.numeric(v[["tx"]]), as.numeric(xy[["tx"]]))
  expect_equal(as.numeric(v[["ty"]]), as.numeric(xy[["ty"]]))

  # Random points agreement
  set.seed(99)
  lon <- runif(200, -180, 180)
  lat <- runif(200, -89.999999, 89.999999)
  for (i in seq_along(lon)) {
    v  <- hexify_forward(lon[i], lat[i])
    xy <- hexify_forward_to_face(as.integer(v[["face"]]), lon[i], lat[i])
    expect_equal(as.numeric(v[["tx"]]), as.numeric(xy[["tx"]]))
    expect_equal(as.numeric(v[["ty"]]), as.numeric(xy[["ty"]]))
  }
})

test_that("hexify_forward_to_face returns finite tx,ty and face-local coords in [0,1]", {
  hexify_build_icosa(11.25, 58.28252559, 0)
  set.seed(123)
  for (i in 1:200) {
    lon <- runif(1, -180, 180)
    lat <- runif(1, -89.999999, 89.999999)
    v   <- hexify_forward(lon, lat)
    xy  <- hexify_forward_to_face(as.integer(v[["face"]]), lon, lat)
    expect_true(is.finite(xy[["tx"]]))
    expect_true(is.finite(xy[["ty"]]))
    expect_true(xy[["tx"]] >= -1e-12 && xy[["tx"]] <= 1 + 1e-12)
    expect_true(xy[["ty"]] >= -1e-12 && xy[["ty"]] <= 1 + 1e-12)
  }
})
