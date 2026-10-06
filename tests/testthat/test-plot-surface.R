# tests/testthat/test-plot-surface.R
# Grid drawn on the sphere and the icosahedron: globe_centers, resolve_center(),
# the face-plane cell boundaries behind both surfaces, and the plot method.

# =============================================================================
# globe_centers and resolve_center()
# =============================================================================

test_that("globe_centers is a named list of lon/lat presets", {
  expect_type(globe_centers, "list")
  expect_true(all(c("europe", "north_america", "south_america", "africa",
                     "asia", "oceania", "middle_east", "south_asia",
                     "pacific", "caribbean", "arctic", "antarctic")
                   %in% names(globe_centers)))

  for (nm in names(globe_centers)) {
    preset <- globe_centers[[nm]]
    expect_length(preset, 2)
    expect_equal(names(preset), c("lon", "lat"), info = nm)
    expect_true(preset["lon"] >= -180 && preset["lon"] <= 180, info = nm)
    expect_true(preset["lat"] >= -90 && preset["lat"] <= 90, info = nm)
  }
})

test_that("globe_centers$europe matches documented value", {
  expect_equal(unname(globe_centers$europe), c(10, 50))
})

test_that("resolve_center accepts a valid preset name", {
  expect_equal(hexify:::resolve_center("europe"), globe_centers$europe)
})

test_that("resolve_center errors on an unknown preset name", {
  expect_error(hexify:::resolve_center("nowhere"), "Unknown center preset")
})

test_that("resolve_center accepts numeric c(lon, lat), named or not", {
  expect_equal(hexify:::resolve_center(c(lon = 5, lat = 40)), c(lon = 5, lat = 40))
  expect_equal(hexify:::resolve_center(c(5, 40)), c(lon = 5, lat = 40))
})

test_that("resolve_center errors on invalid input", {
  expect_error(hexify:::resolve_center(c(1, 2, 3)), "preset name or numeric")
  expect_error(hexify:::resolve_center(TRUE), "preset name or numeric")
})

# =============================================================================
# Camera: orthographic and perspective views
# =============================================================================

test_that("the orthographic view puts the centre at the origin and north up", {
  view <- hexify:::surface_view(c(lon = 40, lat = 20))
  P <- hexify:::unit_vec(c(40, 40), c(20, 30))
  xy <- hexify:::project(P, view)
  expect_equal(xy[1, ], c(0, 0), tolerance = 1e-12)
  expect_equal(xy[2, 1], 0, tolerance = 1e-12)
  expect_gt(xy[2, 2], 0)
  rim <- hexify:::project(hexify:::horizon_ring(view), view)
  expect_equal(rowSums(rim^2), rep(1, nrow(rim)), tolerance = 1e-12)
})

test_that("the default frame holds the unit disc in both projections", {
  ortho <- hexify:::surface_view(c(lon = 0, lat = 0))
  expect_equal(hexify:::view_frame(ortho, NA_real_), c(0, 0, 1.02))
  for (d in c(1.2, 3)) {
    view <- hexify:::surface_view(c(lon = 0, lat = 0), distance = d)
    fov <- hexify:::resolve_camera("perspective", d, 0, 0)$fov
    expect_equal(hexify:::view_frame(view, fov), c(0, 0, 1.02), tolerance = 1e-12,
                 info = paste(d))
  }
})

test_that("the default frame holds the whole visible sphere of a tilted camera", {
  for (cfg in list(c(3, 40), c(1.6, 25), c(2, -20))) {
    info <- paste(cfg, collapse = " ")
    view <- hexify:::surface_view(c(lon = 30, lat = 10), distance = cfg[1], tilt = cfg[2])
    rim <- hexify:::project(hexify:::horizon_ring(view, 721L), view)
    frame <- hexify:::view_frame(view, NA_real_)
    expect_false(anyNA(rim), info = info)
    off <- abs(sweep(rim, 2, frame[1:2]))
    expect_lte(max(off), frame[3] / 1.02 + 1e-12)
    # the rim touches the frame on two opposite sides
    expect_equal(max(apply(rim, 2, function(x) diff(range(x)))), 2 * frame[3] / 1.02,
                 tolerance = 1e-12, info = info)
  }
  # close and steeply tilted, the nearby surface spreads wider than 120
  # degrees: the frame holds the rim cut to that square and touches it on
  # its wider side
  for (cfg in list(c(1.6, 60), c(1.2, 80))) {
    info <- paste(cfg, collapse = " ")
    view <- hexify:::surface_view(c(lon = 30, lat = 10), distance = cfg[1], tilt = cfg[2])
    widest <- view$scale * tan(60 * pi / 180)
    rim <- hexify:::project(hexify:::horizon_ring(view, 721L), view)
    rim <- pmin(pmax(rim[stats::complete.cases(rim), ], -widest), widest)
    frame <- hexify:::view_frame(view, NA_real_)
    expect_lte(frame[3], 1.02 * widest + 1e-12)
    expect_lte(max(abs(sweep(rim, 2, frame[1:2]))), frame[3] / 1.02 + 1e-12)
    expect_equal(max(apply(rim, 2, function(x) diff(range(x)))), 2 * frame[3] / 1.02,
                 tolerance = 1e-12, info = info)
    expect_equal(min(rim), -widest, info = info)
  }
})

test_that("an icosahedron is framed by the corners of the faces the camera sees", {
  view <- hexify:::surface_view(c(lon = 10, lat = 45), distance = 1.8, tilt = 35)
  ico <- hexify:::view_frame(view, NA_real_, hexify:::surface_outline("icosahedron", view, numeric(0)))
  sph <- hexify:::view_frame(view, NA_real_, hexify:::surface_outline("sphere", view, numeric(0)))
  expect_lt(ico[3], sph[3])
  s <- hexify:::icosa_solid(numeric(0))
  xy <- hexify:::project(s$vertices[unique(as.vector(s$faces[hexify:::icosa_front(s, view), ])), ], view)
  expect_lte(max(abs(sweep(xy, 2, ico[1:2]))), ico[3] / 1.02 + 1e-12)
})

test_that("an untilted perspective view shows the sphere as the unit disc", {
  for (d in c(1.2, 3, 6.6)) {
    view <- hexify:::surface_view(c(lon = -70, lat = 45), distance = d)
    rim <- hexify:::project(hexify:::horizon_ring(view), view)
    expect_equal(rowSums(rim^2), rep(1, nrow(rim)), tolerance = 1e-12, info = d)
    expect_equal(drop(hexify:::horizon_ring(view) %*% view$dir),
                 rep(1 / d, nrow(rim)), tolerance = 1e-12, info = d)
  }
})

test_that("faces_camera keeps the sphere points the camera has a line of sight to", {
  set.seed(1)
  P <- matrix(rnorm(3000), ncol = 3)
  P <- P / sqrt(rowSums(P^2))
  view <- hexify:::surface_view(c(lon = 100, lat = -30), distance = 2.2)
  # the ray from the camera meets the sphere first at P when P lies on the
  # near side of the tangent plane through P
  seen <- rowSums((matrix(view$eye, nrow(P), 3, byrow = TRUE) - P) * P) > 0
  expect_identical(hexify:::faces_camera(P, view), seen)
})

test_that("an icosahedron face is drawn exactly when the camera is outside its plane", {
  s <- hexify:::icosa_solid(numeric(0))
  for (d in c(1.05, 1.5, 4)) {
    view <- hexify:::surface_view(c(lon = 10, lat = 60), distance = d)
    h <- rowSums(s$normals * s$vertices[s$faces[, 1], ])
    drawn <- drop(s$normals %*% view$dir) > h * view$horizon
    centroid <- t(apply(s$faces, 1, function(f) colMeans(s$vertices[f, ])))
    outside <- rowSums((matrix(view$eye, 20, 3, byrow = TRUE) - centroid) * s$normals) > 0
    expect_identical(drawn, outside, info = d)
  }
})

test_that("project_ring keeps the part of a ring in front of the camera", {
  view <- hexify:::surface_view(c(lon = 0, lat = 0), distance = 1.2, tilt = 80)
  rim <- hexify:::horizon_ring(view)
  q <- hexify:::camera_coords(rim, view)
  expect_true(any(q[, 3] < view$near) && any(q[, 3] >= view$near))
  xy <- hexify:::project_ring(rim, view)
  expect_false(anyNA(xy))
  # every kept vertex is in front of the near plane, and the two cut points
  # lie on it
  kept <- sum(q[, 3] >= view$near)
  expect_equal(nrow(xy), kept + 2L)
  # a ring wholly behind the camera leaves nothing
  behind <- sweep(diag(3) * 0.1, 2, view$eye + 0.5 * view$cam[, 3], "+")
  expect_null(hexify:::project_ring(behind, view))
})

test_that("orient_ring turns a ring to the asked direction", {
  sq <- cbind(c(0, 1, 1, 0), c(0, 0, 1, 1))
  area2 <- function(p) {
    n <- nrow(p)
    sum(p[, 1] * p[c(2:n, 1), 2] - p[c(2:n, 1), 1] * p[, 2])
  }
  for (p in list(sq, sq[4:1, ])) {
    expect_gt(area2(hexify:::orient_ring(p, TRUE)), 0)
    expect_lt(area2(hexify:::orient_ring(p, FALSE)), 0)
  }
  expect_identical(hexify:::orient_ring(sq, TRUE), sq)
})

test_that("rotation turns north clockwise", {
  north <- hexify:::unit_vec(0, 90)
  view <- hexify:::surface_view(c(lon = 0, lat = 0), rotation = 90)
  expect_equal(drop(hexify:::project(north, view)), c(1, 0), tolerance = 1e-12)
})

test_that("tilt swings the camera about the centre point, which stays in the middle", {
  target <- hexify:::unit_vec(-40, 25)
  for (tilt in c(-60, 0, 35, 80)) {
    info <- paste("tilt", tilt)
    view <- hexify:::surface_view(c(lon = -40, lat = 25), distance = 1.8, tilt = tilt)
    expect_equal(drop(hexify:::project(target, view)), c(0, 0), tolerance = 1e-12, info = info)
    offset <- view$eye - drop(target)
    expect_equal(sqrt(sum(offset^2)), 0.8, tolerance = 1e-12, info = info)
    # the line to the camera leans from the local vertical by the tilt
    expect_equal(sum(offset * drop(target)) / 0.8, cos(tilt * pi / 180),
                 tolerance = 1e-12, info = info)
  }
  # a positive tilt moves the camera south of the point, so it looks north
  view <- hexify:::surface_view(c(lon = 0, lat = 0), distance = 2, tilt = 30)
  expect_lt(view$eye[3], 0)
})

test_that("projection arguments are checked", {
  cam <- hexify:::resolve_camera
  expect_identical(cam("orthographic", NULL, 0, 0), list(distance = Inf, fov = NA_real_))
  expect_identical(cam("perspective", NULL, 0, 0), list(distance = 3, fov = NA_real_))
  expect_identical(cam("perspective", 2, 0, 0, fov = 20)$fov, 20)
  expect_error(cam("orthographic", 2, 0, 0), "perspective camera")
  expect_error(cam("orthographic", NULL, 10, 0), "perspective camera")
  expect_error(cam("orthographic", NULL, 0, 0, fov = 30), "perspective camera")
  expect_error(cam("perspective", 2, 0, 0, fov = 170), "between 0 and 170")
  expect_error(cam("perspective", 1, 0, 0), "greater than 1")
  expect_error(cam("perspective", 2, 90, 0), "between -90 and 90")
  expect_identical(cam("perspective", 1.2, 80, 0)$distance, 1.2)
  expect_error(cam("perspective", 2, NA_real_, 0), "single number")
})

# =============================================================================
# Face-plane cell boundaries
# =============================================================================

surface_paths_for <- function(g, step = 0.02) {
  hexify:::grid_surface_paths(g, NULL, step)
}

to_lonlat <- function(S) {
  cbind(atan2(S[, 2], S[, 1]) * 180 / pi, asin(pmax(-1, pmin(1, S[, 3]))) * 180 / pi)
}

test_that("icosahedron solid has 12 unit vertices, 20 faces and 30 edges", {
  s <- hexify:::icosa_solid(numeric(0))
  expect_equal(dim(s$vertices), c(12L, 3L))
  expect_equal(rowSums(s$vertices^2), rep(1, 12), tolerance = 1e-12)
  expect_equal(dim(s$faces), c(20L, 3L))
  expect_equal(nrow(s$edges), 30L)
  # every edge has the same length on the inscribed icosahedron
  len <- sqrt(rowSums((s$vertices[s$edges[, "v1"], ] - s$vertices[s$edges[, "v2"], ])^2))
  expect_equal(len, rep(len[1], 30), tolerance = 1e-12)
  # neighbouring vertices subtend arccos(1/sqrt(5))
  expect_equal(len[1], sqrt(2 - 2 / sqrt(5)), tolerance = 1e-14)
})

test_that("cell boundary points sit on their face plane and on the sphere", {
  s <- hexify:::icosa_solid(numeric(0))
  inradius <- sqrt(sum(colMeans(s$vertices[s$faces[1, ], ])^2))
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 2, aperture = ap)
    P <- surface_paths_for(g)
    Q <- P[, c("solid_x", "solid_y", "solid_z")]
    S <- P[, c("sphere_x", "sphere_y", "sphere_z")]
    on_plane <- rowSums(Q * s$normals[P[, "face"] + 1, ])
    expect_equal(on_plane, rep(inradius, nrow(P)), tolerance = 1e-7, info = ap)
    expect_equal(rowSums(S^2), rep(1, nrow(P)), tolerance = 1e-12, info = ap)
    expect_setequal(unique(P[, "cell"]), seq_len(n_cells(g)))
  }
})

test_that("each cell's path is closed and continuous", {
  g <- hex_grid(resolution = 3, aperture = 3)
  P <- surface_paths_for(g, step = 0.02)
  S <- P[, c("sphere_x", "sphere_y", "sphere_z")]
  first <- !duplicated(P[, "cell"])
  last <- !duplicated(P[, "cell"], fromLast = TRUE)
  expect_equal(S[first, ], S[last, ], tolerance = 1e-12)
  same <- P[-1, "cell"] == P[-nrow(P), "cell"]
  step_deg <- acos(pmin(1, rowSums(S[-1, ][same, ] * S[-nrow(S), ][same, ]))) * 180 / pi
  # a step of 0.02 face edges is at most about 1.5 degrees of arc
  expect_lt(max(step_deg), 2)
})

test_that("cell boundaries separate the cell from its neighbours", {
  for (cfg in list(list(3, 3), list(4, 2), list(7, 2), list("4/7", 2), list("4/7", 1))) {
    g <- hex_grid(resolution = cfg[[2]], aperture = cfg[[1]])
    P <- surface_paths_for(g, step = 0.05)
    S <- P[, c("sphere_x", "sphere_y", "sphere_z")]
    ids <- seq_len(n_cells(g))
    ctr <- cell_to_lonlat(ids, g)
    C <- hexify:::unit_vec(ctr$lon_deg, ctr$lat_deg)[P[, "cell"], ]
    inward <- S + 0.005 * (C - S)
    q <- to_lonlat(inward / sqrt(rowSums(inward^2)))
    expect_equal(lonlat_to_cell(q[, 1], q[, 2], g), ids[P[, "cell"]],
                 info = paste(cfg[[1]], cfg[[2]]))
  }
})

test_that("a vertex cell of a 4,4,7,3 sequence drops the corner no face reads", {
  g <- hex_grid(resolution = 4, aperture = c(4, 4, 7, 3))
  seq <- hexify:::isea_levels(g@aperture, g@resolution)$ap_seq
  qij <- hexify:::cpp_cell_to_quad_ij(as.numeric(seq_len(n_cells(g))), 4L, 0L, seq)
  vertex <- which(qij$i == 0 & qij$j == 0)
  expect_length(vertex, 12L)
  P <- hexify:::cpp_cell_surface_paths(numeric(0), as.numeric(vertex), 4L, 0L, seq, 0.05)
  expect_setequal(unique(P[, "cell"]), seq_along(vertex))
})

# =============================================================================
# plot(<HexGridInfo>)
# =============================================================================

draws <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(expr)
}

test_that("plot draws a grid on the sphere and on the icosahedron", {
  g <- hex_grid(resolution = 2, aperture = 3)
  expect_identical(draws(plot(g)), g)
  expect_identical(draws(plot(g, surface = "icosahedron")), g)
  expect_identical(draws(plot(g, surface = "icosahedron", land = FALSE,
                              center = "pacific")), g)
  expect_identical(draws(plot(g, cells = 1:5, land_border = NA,
                              center = c(-60, -15))), g)
})

test_that("plot draws perspective, tilted and rotated views", {
  g <- hex_grid(resolution = 2, aperture = 3)
  expect_identical(draws(plot(g, projection = "perspective", distance = 1.5)), g)
  expect_identical(draws(plot(g, surface = "icosahedron", projection = "perspective",
                              distance = 2, tilt = 30, rotation = -45)), g)
  expect_identical(draws(plot(g, rotation = 120, center = "arctic")), g)
  expect_identical(draws(plot(g, projection = "perspective", distance = 1.3,
                              tilt = 30, fov = 40)), g)
  expect_identical(draws(plot(g, surface = "icosahedron", projection = "perspective",
                              distance = 1.15, tilt = 80, fov = 70)), g)
  expect_error(draws(plot(g, distance = 2)), "perspective camera")
})

test_that("plot takes land as an sf object", {
  g <- hex_grid(resolution = 1, aperture = 4)
  europe <- hexify_world[hexify_world$continent == "Europe", ]
  expect_identical(draws(plot(g, surface = "icosahedron", land = europe)), g)
  expect_error(draws(plot(g, land = "europe")), "land must be")
})

test_that("an H3 grid draws on the sphere only", {
  h <- hex_grid(resolution = 0, type = "h3")
  expect_identical(suppressMessages(draws(plot(h, land = FALSE))), h)
  expect_error(draws(plot(h, surface = "icosahedron")), "needs an ISEA grid")
  expect_error(draws(plot(h, face_edges = TRUE)), "not built on")
})

# =============================================================================
# The unfolded net
# =============================================================================

test_that("boundary points of one face run unbroken in the plane", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 2, aperture = ap)
    P <- surface_paths_for(g, step = 0.02)
    n <- nrow(P)
    same <- P[-1, "cell"] == P[-n, "cell"] & P[-1, "face"] == P[-n, "face"]
    gap <- sqrt(rowSums((P[-1, c("plane_x", "plane_y")] -
                         P[-n, c("plane_x", "plane_y")])^2))[same]
    expect_lte(max(gap), 0.02 + 1e-9)
  }
})

test_that("boundary points lie inside their face's triangle of the net", {
  tris <- hexify:::net_triangles(numeric(0))
  P <- surface_paths_for(hex_grid(resolution = 3, aperture = 3))
  for (f in unique(P[, "face"])) {
    T <- tris[[f + 1]]$plane
    p <- P[P[, "face"] == f, c("plane_x", "plane_y"), drop = FALSE]
    M <- cbind(T[1, ] - T[3, ], T[2, ] - T[3, ])
    w <- solve(M, t(p) - T[3, ])
    expect_gte(min(w, 1 - colSums(w)), -1e-9)
  }
})

test_that("a cell on one face is centred on its PLANE centre", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 3, aperture = ap)
    P <- surface_paths_for(g, step = 0.01)
    one_face <- tapply(P[, "face"], P[, "cell"], function(f) length(unique(f)) == 1)
    cells <- as.integer(names(one_face)[one_face])[1:20]
    ctr <- hexify_cell_to_plane(cells, 3, ap)
    for (i in seq_along(cells)) {
      ring <- P[P[, "cell"] == cells[i], c("plane_x", "plane_y")]
      xy <- sf::st_coordinates(sf::st_centroid(sf::st_polygon(list(ring))))
      expect_equal(unname(xy[1, 1:2]), c(ctr$plane_x[i], ctr$plane_y[i]),
                   tolerance = 1e-6, info = paste(ap, cells[i]))
    }
  }
})

test_that("plot draws a grid on the net", {
  g <- hex_grid(resolution = 2, aperture = 3)
  expect_identical(draws(plot(g, surface = "net")), g)
  expect_identical(draws(plot(g, surface = "net", land = FALSE, face_edges = FALSE,
                              cells = 1:10)), g)
  expect_identical(draws(plot(hex_grid(resolution = 3, aperture = "4/3"),
                              surface = "net", land = FALSE)), hex_grid(resolution = 3, aperture = "4/3"))
  expect_error(draws(plot(g, surface = "net", center = "europe")), "drawn flat")
  expect_error(draws(plot(g, surface = "net", rotation = 10)), "drawn flat")
  h <- hex_grid(resolution = 0, type = "h3")
  expect_error(suppressMessages(draws(plot(h, surface = "net"))), "needs an ISEA grid")
})
