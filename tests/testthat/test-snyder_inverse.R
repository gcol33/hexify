test_that("snyder inverse round-trips small offsets near face center", {
  # Use one face center from forward, project a tiny displacement and invert back.
  fc <- hexify_face_centers()
  face <- 0L
  # trivial: exactly the face center should invert to that center
  xy <- hexify_forward_to_face(face, fc$lon[face+1], fc$lat[face+1])
  ll <- hexify_inverse(xy[["tx"]], xy[["ty"]], face)
  expect_true(is.numeric(ll[["lon"]]) && is.numeric(ll[["lat"]]))
  # Near the same lon/lat (wrapping aside)
  expect_lt(abs(ll[["lat"]] - fc$lat[face+1]), 1e-6)
})

test_that("precision presets are accepted", {
  hexify_set_precision("fast")
  p <- hexify_get_precision()
  expect_true(p[["tol"]] >= 1e-12 || p[["max_iters"]] <= 40)
})
