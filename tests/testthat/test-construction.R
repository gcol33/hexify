# Lambert's construction of Snyder's projection, step by step (#98)

unit_rows <- function(m) m / sqrt(rowSums(m^2))

# Points along the great-circle arc from a to b (unit vectors), n of them
arc_points <- function(a, b, n) {
  t <- seq(0, 1, length.out = n)
  unit_rows(outer(1 - t, a) + outer(t, b))
}

shoelace <- function(u, v) abs(sum(u * c(v[-1], v[1]) - c(u[-1], u[1]) * v)) / 2

construction_grids <- function() {
  list(
    hex_grid(resolution = 3, aperture = 3),
    hex_grid(resolution = 3, aperture = 4, orientation = c(-40, 20, 33)),
    hex_grid(resolution = 2, aperture = 4, polyhedron = "octahedron"),
    hex_grid(resolution = 3, aperture = 7, ellipsoid = "WGS84", orientation = "ogc")
  )
}

random_lonlat <- function(n, seed) {
  set.seed(seed)
  cbind(runif(n, -180, 180), asin(runif(n, -1, 1)) * 180 / pi)
}

test_that("the last step is the forward projection, bit for bit", {
  p <- random_lonlat(3000, 981)
  for (g in construction_grids()) {
    s <- projection_stages(g, p[, 1], p[, 2])
    ic <- hexify:::icosa_arg(g)
    fw <- t(vapply(seq_len(nrow(p)), function(k) {
      hexify:::cpp_icosa_forward(ic, p[k, 1], p[k, 2])
    }, numeric(3)))
    expect_identical(s$face, as.integer(fw[, "face"]))
    expect_identical(s$tx, unname(fw[, "icosa_triangle_x"]))
    expect_identical(s$ty, unname(fw[, "icosa_triangle_y"]))
  }
})

test_that("the Lambert point keeps its distance |TL| = |TS| = 2 sin(z / 2)", {
  p <- random_lonlat(3000, 982)
  for (g in construction_grids()) {
    s <- projection_stages(g, p[, 1], p[, 2])
    z <- s$arc * pi / 180
    chord <- 2 * sin(z / 2)
    ts <- sqrt(s$sphere_u^2 + s$sphere_v^2 + (cos(z) - 1)^2)
    tl <- sqrt(s$lambert_u^2 + s$lambert_v^2)
    expect_equal(ts, chord, tolerance = 1e-12)
    expect_equal(tl, chord, tolerance = 1e-12)
    # S, L, T and the sphere's centre lie in one plane: L is S swung down
    # about T, at the azimuth of S
    expect_lt(max(abs(s$sphere_u * s$lambert_v - s$sphere_v * s$lambert_u)), 1e-12)
    expect_gt(min(s$sphere_u * s$lambert_u + s$sphere_v * s$lambert_v), -1e-15)
    # Snyder's point lies at radius 2 R' f sin(z / 2), the nudged one at
    # 2 f sin(z / 2), both at azimuth Az'
    r1 <- attr(s, "r1")
    expect_equal(sqrt(s$snyder_u^2 + s$snyder_v^2), r1 * s$f * chord, tolerance = 1e-12)
    expect_equal(s$nudge_u * r1, s$snyder_u, tolerance = 1e-14)
    a <- s$az_prime * pi / 180
    expect_equal(s$snyder_u, r1 * s$f * chord * sin(a), tolerance = 1e-12)
    expect_equal(s$snyder_v, r1 * s$f * chord * cos(a), tolerance = 1e-12)
  }
})

test_that("each step's plane coordinates are the points the globe draws", {
  # The globe shader builds the steps from the sphere point S, the face
  # centre T and the flat face: L = T + |TS| (S - (S . T) T) / |S - (S . T) T|,
  # Snyder's point at R' / cos(g) times the inscribed solid's point, and the
  # nudged point that over R'.
  p <- random_lonlat(2000, 983)
  for (g in construction_grids()) {
    ic <- hexify:::icosa_arg(g)
    s <- projection_stages(g, p[, 1], p[, 2])
    r1 <- attr(s, "r1")
    k <- r1 / cos(hexify:::solid_info(hexify:::grid_polyhedron(g))$g_deg * pi / 180)
    normals <- hexify:::icosa_solid(ic)$normals
    for (f in unique(s$face)) {
      at <- which(s$face == f)
      corner <- hexify:::cpp_face_tri_to_solid(ic, f, c(0, 1, 0), c(0, 0, 1))
      corner <- unname(corner)
      ex <- (corner[2, ] - corner[1, ]) / sqrt(sum((corner[2, ] - corner[1, ])^2))
      ey <- (corner[3, ] - corner[1, ]) / sqrt(sum((corner[3, ] - corner[1, ])^2))
      tc <- unname(normals[f + 1, ])
      expect_equal(sum(ex * ey), 0, tolerance = 1e-12)
      expect_equal(sum(ex * tc), 0, tolerance = 1e-12)
      # x, y and the outward normal turn as the azimuth does, clockwise seen
      # from outside
      expect_equal(sum(hexify:::cross3(ex, ey) * tc), 1, tolerance = 1e-12)
      place <- function(u, v, h) {
        outer(rep(h, length.out = length(u)), tc) + outer(u, ex) + outer(v, ey)
      }
      S <- unname(hexify:::unit_vec(p[at, 1], p[at, 2], ic))
      z <- s$arc[at] * pi / 180
      expect_equal(place(s$sphere_u[at], s$sphere_v[at], cos(z)), S, tolerance = 1e-12)
      off <- S - outer(drop(S %*% tc), tc)
      lambert <- outer(rep(1, length(at)), tc) +
        sqrt(rowSums((S - outer(rep(1, length(at)), tc))^2)) * unit_rows(off)
      expect_equal(place(s$lambert_u[at], s$lambert_v[at], 1), lambert, tolerance = 1e-12)
      solid <- hexify:::cpp_face_tri_to_solid(ic, f, s$tx[at], s$ty[at])
      expect_equal(place(s$snyder_u[at], s$snyder_v[at], r1), k * unname(solid),
                   tolerance = 1e-12)
      expect_equal(place(s$nudge_u[at], s$nudge_v[at], 1), k / r1 * unname(solid),
                   tolerance = 1e-12)
    }
  }
})

test_that("the Lambert step is equal-area and Snyder's steps keep area in proportion", {
  for (poly in c("icosahedron", "octahedron")) {
    g <- hex_grid(resolution = 2, aperture = 4, polyhedron = poly)
    ic <- hexify:::icosa_arg(g)
    solid <- hexify:::icosa_solid(ic)
    V <- unit_rows(solid$vertices)
    n_faces <- nrow(solid$faces)
    r1 <- attr(projection_stages(g, 0, 0), "r1")
    for (f in c(0L, n_faces - 1L)) {
      v <- solid$faces[f + 1, ]
      centre <- unit_rows(matrix(colMeans(V[v, ]), 1))[1, ]
      mid <- unit_rows(matrix(V[v[1], ] + V[v[2], ], 1))[1, ]
      # The right triangle (centre, vertex, edge midpoint), a sixth of the face
      ring <- rbind(arc_points(centre, V[v[1], ], 2000),
                    arc_points(V[v[1], ], mid, 2000),
                    arc_points(mid, centre, 2000))
      ll <- hexify:::vec_lonlat(ring, ic)
      s <- projection_stages(g, ll[, 1], ll[, 2], face = f)
      sixth <- 4 * pi / n_faces / 6
      expect_equal(shoelace(s$lambert_u, s$lambert_v), sixth, tolerance = 1e-6)
      expect_equal(shoelace(s$nudge_u, s$nudge_v), sixth / r1^2, tolerance = 1e-9)
      expect_equal(shoelace(s$snyder_u, s$snyder_v), sixth, tolerance = 1e-9)
      # The whole face: Lambert's image has curved edges and the face's area
      edge <- do.call(rbind, lapply(1:3, function(k) {
        arc_points(V[v[k], ], V[v[k %% 3 + 1], ], 4000)
      }))
      ll <- hexify:::vec_lonlat(edge, ic)
      s <- projection_stages(g, ll[, 1], ll[, 2], face = f)
      expect_equal(shoelace(s$lambert_u, s$lambert_v), 4 * pi / n_faces, tolerance = 1e-6)
      expect_equal(shoelace(s$snyder_u, s$snyder_v), 4 * pi / n_faces, tolerance = 1e-9)
      # Snyder's edge is straight, Lambert's curved: the distance of the
      # edge's middle point from the chord between its ends, over the chord
      i <- c(1, 2000, 4000)
      bow <- function(u, v) {
        du <- u[i[3]] - u[i[1]]
        dv <- v[i[3]] - v[i[1]]
        abs((u[i[2]] - u[i[1]]) * dv - (v[i[2]] - v[i[1]]) * du) / (du^2 + dv^2)
      }
      expect_lt(bow(s$snyder_u, s$snyder_v), 1e-12)
      expect_gt(bow(s$lambert_u, s$lambert_v), 0.01)
    }
  }
})

test_that("Tissot's indicatrix at each step", {
  p <- random_lonlat(2000, 984)
  for (g in construction_grids()) {
    s <- projection_stages(g, p[, 1], p[, 2])
    z <- s$arc * pi / 180
    lam <- projection_distortion(g, p[, 1], p[, 2], stage = "lambert")
    nud <- projection_distortion(g, p[, 1], p[, 2], stage = "nudge")
    face <- projection_distortion(g, p[, 1], p[, 2])
    r1 <- attr(s, "r1")
    if (hexify:::grid_flattening(g) == 0) {
      # Lambert's azimuthal equal-area: 1 / cos(z / 2) across the radius,
      # cos(z / 2) along it
      expect_equal(lam$a, 1 / cos(z / 2), tolerance = 1e-12)
      expect_equal(lam$b, cos(z / 2), tolerance = 1e-12)
    }
    expect_equal(lam$areal, rep(1, nrow(p)), tolerance = 1e-12)
    expect_equal(nud$a, face$a / r1, tolerance = 1e-12)
    expect_equal(nud$b, face$b / r1, tolerance = 1e-12)
    expect_equal(nud$areal, rep(1 / r1^2, nrow(p)), tolerance = 1e-12)
    expect_equal(nud$angular, face$angular, tolerance = 1e-10)
  }
})

test_that("Snyder's azimuth agrees with Lambert's along the radii to the vertices and edge midpoints", {
  g <- hex_grid(resolution = 2, aperture = 3)
  ic <- hexify:::icosa_arg(g)
  solid <- hexify:::icosa_solid(ic)
  V <- unit_rows(solid$vertices)
  v <- solid$faces[1, ]
  centre <- unit_rows(matrix(colMeans(V[v, ]), 1))[1, ]
  mid <- unit_rows(matrix(V[v[1], ] + V[v[2], ], 1))[1, ]
  ray <- rbind(arc_points(centre, V[v[1], ], 50)[-1, ], arc_points(centre, mid, 50)[-1, ])
  ll <- hexify:::vec_lonlat(ray, ic)
  s <- projection_stages(g, ll[, 1], ll[, 2], face = 0L)
  turn <- (s$az_prime - s$az + 180) %% 360 - 180
  expect_lt(max(abs(turn)), 1e-9)
})

test_that("projection_stages checks its arguments", {
  g <- hex_grid(resolution = 2, aperture = 3)
  expect_error(projection_stages(hex_grid(resolution = 2, aperture = 3, projection = "fuller"),
                                 0, 0), "isea")
  expect_error(projection_stages(hex_grid(resolution = 1, type = "h3"), 0, 0), "H3")
  expect_error(projection_stages(g, 0, c(0, 1)), "same length")
  expect_error(projection_stages(g, 0, 0, face = 20), "face")
  expect_error(projection_stages(g, c(0, 1), c(0, 1), face = c(1, 2, 3)), "face")
  expect_error(projection_distortion(hex_grid(resolution = 2, aperture = 3, projection = "ivea"),
                                     0, 0, stage = "nudge"), "isea")
  expect_equal(nrow(projection_stages(g, numeric(0), numeric(0))), 0)
})

test_that("a face drawn at each step lies where projection_stages() puts it", {
  g <- hex_grid(resolution = 2, aperture = 3)
  r1 <- attr(projection_stages(g, 0, 0), "r1")
  snyder <- hexify:::stage_outline(g, 0L, "snyder")
  face <- hexify:::stage_outline(g, 0L, "face")
  edge <- sqrt(16 * pi / (sqrt(3) * 20))
  # The plane triangle's corners at R' tan(g) from the centre
  g_deg <- hexify:::solid_info("icosahedron")$g_deg
  expect_equal(max(sqrt(rowSums(snyder^2))), edge / sqrt(3), tolerance = 1e-9)
  expect_equal(max(sqrt(rowSums(face^2))), r1 * tan(g_deg * pi / 180), tolerance = 1e-9)
  lambert <- hexify:::stage_outline(g, 0L, "lambert")
  expect_equal(max(sqrt(rowSums(lambert^2))), 2 * sin(g_deg * pi / 360), tolerance = 1e-9)
  walls <- hexify:::stage_walls(g, 0L, "nudge")
  inside <- walls[stats::complete.cases(walls), , drop = FALSE] * r1
  tri <- sweep(inside / edge, 2, c(0.5, 1 / (2 * sqrt(3))), "+")
  expect_true(all(hexify:::in_triangle(tri, hexify:::FACE_TRIANGLE, tol = 1e-6)))
  e <- hexify:::stage_ellipses(hex_grid(resolution = 2, aperture = 3, projection = "fuller"),
                               0L, "face")
  expect_gt(max(attr(e, "angular")), 0)
  expect_error(hexify:::stage_walls(hex_grid(resolution = 2, aperture = 3, projection = "ivea"),
                                    0L, "lambert"), "isea")
})
