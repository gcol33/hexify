# Distortion of the face projections (#88)

# The point at arc offset `eps` (in face-plane units) from face 0's centre
# towards its first vertex, in lon/lat
near_centre_lonlat <- function(g, eps = 1e-7) {
  ll <- cpp_face_xy_to_ll(icosa_arg(g), 0.5, 1 / (2 * sqrt(3)) + eps, 0L)
  c(ll[["lon"]], ll[["lat"]])
}

test_that("Snyder's projection reproduces Snyder's (1992) Table 1 at the face centre", {
  g <- hex_grid(resolution = 1, aperture = 3)
  p <- near_centre_lonlat(g)
  d <- projection_distortion(g, p[1], p[2])
  expect_equal(round(d$angular, 2), 17.27)
  expect_equal(round(d$a, 3), 1.163)
  expect_equal(round(d$b, 3), 0.860)
})

test_that("Snyder's angular deformation is largest at the face centre", {
  g <- hex_grid(resolution = 1, aperture = 3)
  ic <- icosa_arg(g)
  n <- 120
  ij <- expand.grid(i = 0:n, j = 0:n)
  ij <- ij[ij$i + ij$j <= n, ]
  x <- (ij$i + 0.5 * ij$j) / n
  y <- ij$j * sqrt(3) / 2 / n
  s <- cpp_face_tri_tissot(ic, 0L, x, y)
  peak <- projection_distortion(g, near_centre_lonlat(g)[1], near_centre_lonlat(g)[2])
  expect_lte(max(2 * asin((s$a - s$b) / (s$a + s$b)) * 180 / pi), peak$angular + 1e-6)
})

test_that("Snyder's projection and IVEA keep areas on every solid", {
  set.seed(88)
  lon <- runif(500, -180, 180)
  lat <- asin(runif(500, -1, 1)) * 180 / pi
  for (proj in c("isea", "ivea")) for (poly in c("icosahedron", "octahedron", "tetrahedron")) {
    g <- hex_grid(resolution = 1, aperture = 4, polyhedron = poly, projection = proj)
    d <- projection_distortion(g, lon, lat)
    expect_lt(max(abs(d$areal - 1)), 1e-12)
    expect_true(all(d$a >= d$b))
  }
})

test_that("Fuller's areal scale averages 1 over the sphere", {
  g <- hex_grid(resolution = 1, aperture = 3, projection = "fuller")
  # Fibonacci points: equal areas of the sphere
  n <- 20000
  k <- seq_len(n) - 0.5
  lat <- asin(1 - 2 * k / n) * 180 / pi
  lon <- ((k * 180 * (3 - sqrt(5))) %% 360) - 180
  d <- projection_distortion(g, lon, lat)
  expect_equal(mean(d$areal), 1, tolerance = 1e-3)
  expect_gt(sd(d$areal), 0.01)
  expect_gt(max(d$angular), 5)
})

test_that("the exact scale factors match finite differences", {
  finite_scale <- function(g, lon, lat, h = 1e-6) {
    ic <- icosa_arg(g)
    face <- cpp_which_face(ic, lon, lat)
    at <- function(lo, la) {
      unname(cpp_lonlat_to_face_solid(ic, face, lo, la)[1, c("tx", "ty")])
    }
    hd <- h * 180 / pi
    east <- (at(lon + hd / cos(lat * pi / 180), lat) -
               at(lon - hd / cos(lat * pi / 180), lat)) / (2 * h)
    north <- (at(lon, lat + hd) - at(lon, lat - hd)) / (2 * h)
    sv <- svd(cbind(east, north) * face_plane_edge(ic))$d
    c(a = sv[1], b = sv[2])
  }
  pts <- rbind(c(23.4, 41.2), c(-120.5, -12.3), c(77.7, 66.6), c(5.1, -50.2))
  for (proj in c("isea", "fuller", "ivea")) {
    g <- hex_grid(resolution = 1, aperture = 3, projection = proj)
    d <- projection_distortion(g, pts[, 1], pts[, 2])
    for (i in seq_len(nrow(pts))) {
      fd <- finite_scale(g, pts[i, 1], pts[i, 2])
      expect_equal(c(d$a[i], d$b[i]), unname(fd), tolerance = 1e-6)
    }
  }
})

test_that("projection_distortion() refuses H3 and bad points", {
  expect_error(projection_distortion(hex_grid(resolution = 3, type = "h3"), 0, 0),
               "H3")
  g <- hex_grid(resolution = 1, aperture = 3)
  expect_error(projection_distortion(g, 1:2, 1), "same length")
  expect_error(projection_distortion(g, NA_real_, 0), "finite")
})

# Every segment of a path on one face, as rows (face, tx, ty), crosses no
# crease of the projection strictly between its ends: the lines from the face
# centre to the corners, and on IVEA to the edge midpoints as well
crosses_crease <- function(face, tx, ty, path = rep(1, length(tx)),
                           projection = "isea") {
  centre <- c(0.5, 1 / (2 * sqrt(3)))
  corner <- rbind(c(0.5, sqrt(3) / 2), c(0, 0), c(1, 0))
  if (projection == "ivea") {
    corner <- rbind(corner, c(0.25, sqrt(3) / 4), c(0.5, 0), c(0.75, sqrt(3) / 4))
  }
  n <- length(tx)
  same <- face[-1] == face[-n] & path[-1] == path[-n]
  bad <- FALSE
  for (k in seq_len(nrow(corner))) {
    e <- corner[k, ] - centre
    side <- e[1] * (ty - centre[2]) - e[2] * (tx - centre[1])
    s0 <- side[-n]
    s1 <- side[-1]
    hit <- same & s0 * s1 < 0 & abs(s0) > 1e-9 & abs(s1) > 1e-9
    if (!any(hit)) next
    # where the segment crosses the line, along the line from the centre
    t <- s0[hit] / (s0[hit] - s1[hit])
    px <- tx[-n][hit] + t * (tx[-1][hit] - tx[-n][hit])
    py <- ty[-n][hit] + t * (ty[-1][hit] - ty[-n][hit])
    r <- ((px - centre[1]) * e[1] + (py - centre[2]) * e[2]) / sum(e^2)
    bad <- bad || any(r > 0 & r < 1)
  }
  bad
}

test_that("cell boundaries on the faces have a vertex on every crease they cross", {
  for (proj in c("isea", "ivea")) {
    g <- hex_grid(resolution = 2, aperture = 3, projection = proj)
    p <- grid_surface_paths(g, NULL, step = 0.5)
    expect_false(crosses_crease(p[, "face"], p[, "tx"], p[, "ty"], p[, "cell"], proj),
                 label = proj)
  }
})

test_that("densified cell boundaries have a vertex on every crease they cross", {
  for (proj in c("isea", "ivea")) {
    g <- hex_grid(resolution = 2, aperture = 3, projection = proj)
    ic <- icosa_arg(g)
    rings <- isea_cell_rings(seq_len(n_cells(g)), g@resolution, g@aperture, ic)
    hits <- vapply(rings, function(r) {
      n <- nrow(r)
      mid <- (unit_vec(r[-n, 1], r[-n, 2]) + unit_vec(r[-1, 1], r[-1, 2]))
      face <- point_faces(mid / sqrt(rowSums(mid^2)), ic)
      any(vapply(seq_len(n - 1L), function(i) {
        t <- cpp_lonlat_to_face_solid(ic, face[i], r[i + 0:1, 1], r[i + 0:1, 2])
        crosses_crease(c(face[i], face[i]), t[, "tx"], t[, "ty"], projection = proj)
      }, logical(1)))
    }, logical(1))
    expect_false(any(hits), label = proj)
  }
})

test_that("the face mesh is cut along the creases", {
  for (proj in c("isea", "ivea")) {
    g <- hex_grid(resolution = 1, aperture = 3, projection = proj)
    m <- cpp_globe_faces(icosa_arg(g), 0.2)
    idx <- matrix(m$index + 1L, ncol = 3L, byrow = TRUE)
    tri <- matrix(m$tri, ncol = 2L, byrow = TRUE)
    face <- m$item[idx[, 1]] - 1L
    bad <- vapply(seq_len(nrow(idx)), function(k) {
      crosses_crease(rep(face[k], 4L), tri[idx[k, c(1:3, 1)], 1], tri[idx[k, c(1:3, 1)], 2],
                     projection = proj)
    }, logical(1))
    expect_false(any(bad), label = proj)
  }
})

test_that("the plot method draws distortion and Tissot's indicatrix", {
  pdf(NULL)
  on.exit(dev.off())
  g <- hex_grid(resolution = 2, aperture = 3)
  expect_invisible(plot(g, surface = "net", land = FALSE, distortion = "angular",
                        tissot = 30, graticule = 15, cells = numeric(0)))
  expect_invisible(plot(g, surface = "solid", land = FALSE, distortion = "areal",
                        tissot = TRUE))
  expect_invisible(plot(g, land = FALSE, distortion = "angular"))
  gf <- hex_grid(resolution = 2, aperture = 3, projection = "fuller")
  expect_invisible(plot(gf, surface = "net", land = FALSE, distortion = "areal"))
  expect_error(plot(g, tissot = TRUE, land = FALSE), "net")
  expect_error(plot(hex_grid(resolution = 2, type = "h3"), distortion = "angular",
                    land = FALSE), "H3")
})

test_that("distortion classes cover the values", {
  s <- distortion_scale("angular", c(0.2, 17.3))
  expect_equal(range(s$breaks), c(0, 18))
  expect_length(s$col, 18)
  s <- distortion_scale("areal", c(0.95, 1.07))
  expect_equal(range(s$breaks), c(0.93, 1.07))
})
