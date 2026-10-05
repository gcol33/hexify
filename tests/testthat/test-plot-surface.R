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
# Face-plane cell boundaries
# =============================================================================

surface_paths_for <- function(g, step = 0.02) {
  hexify:::grid_surface_paths(g, NULL, step)
}

to_lonlat <- function(S) {
  cbind(atan2(S[, 2], S[, 1]) * 180 / pi, asin(pmax(-1, pmin(1, S[, 3]))) * 180 / pi)
}

test_that("icosahedron solid has 12 unit vertices, 20 faces and 30 edges", {
  s <- hexify:::icosa_solid()
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
  s <- hexify:::icosa_solid()
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
  seq <- hexify:::grid_ap_seq(g)
  qij <- hexify:::cpp_cell_to_quad_ij_seq(as.numeric(seq_len(n_cells(g))), seq)
  vertex <- which(qij$i == 0 & qij$j == 0)
  expect_length(vertex, 12L)
  P <- hexify:::cpp_cell_surface_paths(as.numeric(vertex), 4L, 0L, seq, 0.05)
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
