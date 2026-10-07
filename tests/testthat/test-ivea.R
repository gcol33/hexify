# van Leeuwen and Strebe's (2006) vertex-oriented equal-area projection
# (IVEA). Agreement with DGGAL's IVEA grids on the full set of zones is
# paper/bench/bench_dggal_agreement.R; the fixture holds a sample.

dggal_orientation <- c(11.20, 58.282525588538995, 0)

test_that("hex_grid() carries the projection on every grid solid", {
  g <- hex_grid(resolution = 4, aperture = 3, projection = "ivea")
  expect_identical(g@projection, "ivea")
  expect_identical(icosa_arg(g)[4], 2)
  expect_identical(summary(g)$projection, "ivea")
  expect_output(print(g), "IVEA")
  expect_output(print(g), "Area: ")
  go <- hex_grid(resolution = 3, aperture = 4, projection = "ivea", polyhedron = "octahedron")
  expect_identical(go@projection, "ivea")
  expect_error(hex_grid(resolution = 4, type = "h3", projection = "ivea"),
               "H3 fixes its own projection")
})

test_that("forward then closed-form inverse returns the point on every solid", {
  set.seed(92)
  lon <- runif(1500, -180, 180)
  lat <- asin(runif(1500, -1, 1)) * 180 / pi
  for (solid in c("icosahedron", "octahedron", "tetrahedron")) {
    back <- t(mapply(function(lo, la) {
      f <- hexify_forward(lo, la, projection = "ivea", polyhedron = solid)
      hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]],
                     as.integer(f[["face"]]), projection = "ivea", polyhedron = solid)
    }, lon, lat))
    expect_lt(max(arc_between(lon, lat, back[, "lon"], back[, "lat"])), 1e-13,
              label = solid)
  }
})

test_that("the closed-form inverse matches Newton's method and undoes the forward", {
  icosa <- projection_icosa("ivea", "icosahedron")
  pts <- face_plane_samples(12)
  cx <- 0.5
  cy <- 0.5 / sqrt(3)
  # corners, edge midpoints, points beside the centre and on the lines to
  # corners and midpoints, and the face sampled with its edges and radii
  xy <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2),
              c(0.5, 0), c(0.25, sqrt(3) / 4), c(0.75, sqrt(3) / 4),
              c(cx, cy + 1e-6), c(cx + 1e-6, cy), c(cx, 0.6), c(0.4, cy - 0.1 / sqrt(3)),
              cbind(pts$x, pts$y))
  for (face in c(0L, 7L, 19L)) {
    closed <- t(apply(xy, 1, function(p) cpp_face_xy_to_ll(icosa, p[1], p[2], face)))
    newton <- t(apply(xy, 1, function(p) cpp_face_xy_to_ll(icosa, p[1], p[2], face,
                                                            newton = TRUE)))
    err <- arc_between(closed[, "lon"], closed[, "lat"], newton[, "lon"], newton[, "lat"])
    expect_lt(max(err), 1e-13, label = sprintf("face %d", face))
    back <- t(apply(closed, 1, function(p)
      hexify_forward_to_face(face, p[1], p[2], projection = "ivea")))
    expect_lt(max(abs(back - xy)), 1e-12, label = sprintf("face %d", face))
  }
})

test_that("a face's vertices, edge midpoints and centre keep their places", {
  icosa <- projection_icosa("ivea", "icosahedron")
  for (xy in list(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2), c(0.5, 0), c(0.5, 0.5 / sqrt(3)))) {
    isea <- cpp_face_xy_to_ll(projection_icosa("isea", "icosahedron"), xy[1], xy[2], 3L)
    ivea <- cpp_face_xy_to_ll(icosa, xy[1], xy[2], 3L)
    expect_lt(arc_between(isea[["lon"]], isea[["lat"]], ivea[["lon"]], ivea[["lat"]]),
              1e-13)
  }
})

test_that("IVEA cells are equal-area on the icosahedron and the octahedron", {
  for (solid in c("icosahedron", "octahedron")) for (ap in list(3, 4, 7, "4/3")) {
    res <- if (identical(ap, 7)) 2 else 3
    g <- hex_grid(resolution = res, aperture = ap, projection = "ivea", polyhedron = solid)
    ids <- as_cell_id(seq_len(n_cells(g)))
    lv <- isea_levels(g@aperture, g@resolution)
    sa <- cpp_cell_solid_angle(icosa_arg(g), ids, lv$resolution, lv$aperture, lv$ap_seq, 1e-4)
    sides <- isea_cell_sides(ids, g)
    hex <- 4 * pi / (length(ids) - 2)
    expect_equal(sa, sides / 6 * hex, tolerance = 1e-6,
                 label = sprintf("%s aperture %s", solid, format_aperture(ap, res)))
    a <- cell_area(ids, g)
    expect_equal(sum(a), body_surface_km2(EARTH_RADIUS_KM), tolerance = 1e-12)
  }
})

test_that("an IVEA cell's centre lies in its cell and moves off ISEA's", {
  g <- hex_grid(resolution = 7, aperture = 3, projection = "ivea")
  set.seed(3)
  lon <- runif(500, -180, 180)
  lat <- asin(runif(500, -1, 1)) * 180 / pi
  ids <- lonlat_to_cell(lon, lat, g)
  ctr <- cell_to_lonlat(ids, g)
  expect_identical(lonlat_to_cell(ctr[[1]], ctr[[2]], g), ids)
  gi <- hex_grid(resolution = 7, aperture = 3)
  expect_false(isTRUE(all.equal(cell_to_lonlat(ids, gi), ctr)))
})

test_that("IVEA cells and corners match DGGAL's", {
  d <- read.csv(test_path("data", "dggal_ivea.csv"))
  for (name in unique(d$dggrs)) {
    di <- d[d$dggrs == name, ]
    g <- hex_grid(resolution = di$level[1], aperture = if (name == "IVEA3H") 3 else 7,
                  projection = "ivea", orientation = dggal_orientation)
    ce <- di[di$kind == "centroid", ]
    id <- lonlat_to_cell(ce$lon, ce$lat, g)
    expect_false(anyDuplicated(id) > 0, label = name)
    hc <- cell_to_lonlat(id, g)
    expect_lt(max(arc_between(ce$lon, ce$lat, hc[[1]], hc[[2]])), 1e-12, label = name)

    rings <- isea_cell_rings(id, g@resolution, g@aperture, icosa_arg(g), 0)
    gap <- vapply(seq_len(nrow(ce)), function(i) {
      v <- di[di$kind == "vertex" & di$zone == ce$zone[i], ]
      h <- unique(rings[[i]][, 1:2, drop = FALSE])
      if (nrow(v) != nrow(h)) return(Inf)
      max(vapply(seq_len(nrow(v)), function(j)
        min(arc_between(rep(v$lon[j], nrow(h)), rep(v$lat[j], nrow(h)), h[, 1], h[, 2])),
        numeric(1)))
    }, numeric(1))
    expect_lt(max(gap), 1e-12, label = name)
  }
})

test_that("IVEA has no DGGRID form", {
  g <- hex_grid(resolution = 5, aperture = 4, projection = "ivea")
  expect_error(as_dggrid(g), "DGGRID has no IVEA")
})
