# Fuller's face projection (Gray 1995, Crider 2008). Agreement with DGGRID's
# FULLER grids is in test-orientation.R.

test_that("hex_grid() carries the projection", {
  g <- hex_grid(resolution = 4, aperture = 4, projection = "fuller")
  expect_identical(g@projection, "fuller")
  expect_identical(hex_grid(resolution = 4)@projection, "isea")
  expect_identical(icosa_arg(g)[4], 1)
  expect_identical(summary(g)$projection, "fuller")
  expect_output(print(g), "Fuller")
  expect_error(hex_grid(resolution = 4, type = "h3", projection = "fuller"),
               "H3 fixes its own projection")
  expect_error(hex_grid(resolution = 4, projection = "gnomonic"))
})

test_that("forward then closed-form inverse returns the point at machine precision", {
  set.seed(1)
  lon <- runif(2000, -180, 180)
  lat <- asin(runif(2000, -1, 1)) * 180 / pi
  back <- t(mapply(function(lo, la) {
    f <- hexify_forward(lo, la, projection = "fuller")
    hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]],
                   as.integer(f[["face"]]), projection = "fuller")
  }, lon, lat))
  expect_lt(max(arc_between(lon, lat, back[, "lon"], back[, "lat"])), 1e-13)
})

test_that("the closed-form inverse holds at the vertices, the centre and outside the face", {
  # At the centre two roots of the cubic meet, away from the face's own root.
  # Points up to a quarter of the face outside it are in the closed form's
  # range.
  icosa <- projection_icosa("fuller", "icosahedron")
  pts <- face_plane_samples(12)
  cx <- 0.5
  cy <- 0.5 / sqrt(3)
  xy <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2), c(cx, cy),
              c(cx, cy + 1e-9), c(cx + 1e-9, cy),
              cbind(cx + 1.1 * (pts$x - cx), cy + 1.1 * (pts$y - cy)),
              cbind(cx + 1.25 * (pts$x - cx), cy + 1.25 * (pts$y - cy)))
  for (face in c(0L, 7L, 19L)) {
    closed <- t(apply(xy, 1, function(p) cpp_face_xy_to_ll(icosa, p[1], p[2], face)))
    newton <- t(apply(xy, 1, function(p) cpp_face_xy_to_ll(icosa, p[1], p[2], face,
                                                            newton = TRUE)))
    err <- arc_between(closed[, "lon"], closed[, "lat"], newton[, "lon"], newton[, "lat"])
    expect_lt(max(err), 1e-13, label = sprintf("face %d", face))
  }
})

test_that("Fuller keeps lengths along face edges", {
  # A point a fraction s along a spherical face edge lands the same fraction
  # along the plane triangle's edge.
  fc <- cpp_face_centers(numeric(0))
  for (s in c(0.1, 0.37, 0.5, 0.81)) {
    f0 <- hexify_forward(fc$lon[1] * 180 / pi, fc$lat[1] * 180 / pi, projection = "fuller")
    v1 <- hexify_inverse(0, 0, 0, projection = "fuller")
    v2 <- hexify_inverse(1, 0, 0, projection = "fuller")
    a <- hexify:::unit_vec(v1[[1]], v1[[2]]); b <- hexify:::unit_vec(v2[[1]], v2[[2]])
    arc <- acos(sum(a * b))
    p <- (sin((1 - s) * arc) * a + sin(s * arc) * b) / sin(arc)
    ll <- hexify:::vec_lonlat(p)
    xy <- hexify_forward_to_face(0, ll[1], ll[2], projection = "fuller")
    expect_equal(unname(xy), c(s, 0), tolerance = 1e-9)
    expect_equal(unname(f0[2:3]), c(0.5, 0.5 / sqrt(3)), tolerance = 1e-12)
  }
})

test_that("Fuller cells tile the sphere and differ in area", {
  for (ap in list(3, 4, 7, "4/3")) {
    g <- hex_grid(resolution = 4, aperture = ap, projection = "fuller")
    ids <- grid_global(g)$cell_id
    a <- cell_area(ids, g)
    expect_equal(sum(a), body_surface_km2(EARTH_RADIUS_KM), tolerance = 1e-9)
    expect_gt(max(a) / min(a), 1.1)
    expect_equal(mean(a), g@area_km2, tolerance = 1e-9)
  }
  gi <- hex_grid(resolution = 4, aperture = 4)
  ids <- grid_global(gi)$cell_id
  expect_equal(max(cell_area(ids, gi)),
               body_surface_km2(EARTH_RADIUS_KM) / (length(ids) - 2))
})

test_that("a Fuller cell's centre lies in its cell and shares its ID", {
  g <- hex_grid(resolution = 7, aperture = 3, projection = "fuller")
  set.seed(2)
  lon <- runif(500, -180, 180)
  lat <- asin(runif(500, -1, 1)) * 180 / pi
  ids <- lonlat_to_cell(lon, lat, g)
  ctr <- cell_to_lonlat(ids, g)
  expect_identical(lonlat_to_cell(ctr[[1]], ctr[[2]], g), ids)
  gi <- hex_grid(resolution = 7, aperture = 3)
  expect_false(isTRUE(all.equal(cell_to_lonlat(ids, gi), ctr)))
})

test_that("Fuller travels through the DGGRID conversions", {
  g <- hex_grid(resolution = 5, aperture = 4, projection = "fuller")
  d <- as_dggrid(g)
  expect_identical(d$projection, "FULLER")
  expect_identical(extract_grid(from_dggrid(d))@projection, "fuller")
})
