# tests/testthat/test-net-layout.R
# Net layouts: the PLANE net, Van de Sande's Gosper World on the octahedron,
# the land-joined net, net_project(), and the "gosper" orientation.

octa_gosper <- function(resolution = 2) {
  hex_grid(resolution = resolution, aperture = 4, polyhedron = "octahedron",
           orientation = "gosper")
}

piece_area <- function(p) {
  r <- p$region
  abs(det(cbind(r[2, ] - r[1, ], r[3, ] - r[1, ]))) / 2
}

# Points spread evenly over the sphere (a Fibonacci lattice), as lon/lat
sphere_points <- function(n) {
  i <- seq_len(n) - 0.5
  cbind(lon = ((i * 180 * (3 - sqrt(5))) %% 360) - 180, lat = asin(1 - 2 * i / n) * 180 / pi)
}

# =============================================================================
# The PLANE layout
# =============================================================================

test_that("the plane layout puts every face where DGGRID's PLANE layout does", {
  grids <- list(hex_grid(resolution = 2, aperture = 3),
                hex_grid(resolution = 2, aperture = 3, orientation = c(30, 20, 40)),
                hex_grid(resolution = 2, aperture = 4, polyhedron = "octahedron"),
                hex_grid(resolution = 2, aperture = 3, orientation = "dymaxion",
                         projection = "fuller"))
  for (g in grids) {
    net <- net_layout(g, "plane")
    icosa <- hexify:::icosa_arg(g)
    P <- hexify:::grid_surface_paths(g, NULL, 0.05)
    ref <- hexify:::cpp_icosa_tri_to_plane(icosa, as.integer(P[, "face"]), P[, "tx"], P[, "ty"])
    for (p in net$pieces) {
      at <- P[, "face"] == p$face
      xy <- hexify:::place_points(P[at, c("tx", "ty"), drop = FALSE], p)
      expect_equal(unname(xy), unname(cbind(ref$plane_x[at], ref$plane_y[at])),
                   tolerance = 1e-12)
    }
  }
})

test_that("cell boundaries on the plane net are the PLANE coordinates of the paths", {
  g <- hex_grid(resolution = 2, aperture = 3)
  P <- hexify:::grid_surface_paths(g, NULL, 0.05)
  seg <- hexify:::net_segments(P, net_layout(g, "plane"))
  n <- nrow(P)
  s <- which(P[-1L, "cell"] == P[-n, "cell"] & P[-1L, "face"] == P[-n, "face"])
  ref <- hexify:::cpp_icosa_tri_to_plane(numeric(0), as.integer(P[, "face"]), P[, "tx"], P[, "ty"])
  old <- cbind(ref$plane_x[s], ref$plane_y[s], ref$plane_x[s + 1L], ref$plane_y[s + 1L])
  key <- function(m) do.call(order, as.data.frame(round(m, 9)))
  expect_equal(seg[key(seg), ], unname(old[key(old), ]), tolerance = 1e-12)
})

# =============================================================================
# Pieces
# =============================================================================

test_that("every piece is placed by a rotation and a translation", {
  g <- octa_gosper()
  for (name in c("plane", "gosper", "gosper_flower", "gosper_land", "land")) {
    for (p in net_layout(g, name)$pieces) {
      expect_equal(crossprod(p$A), diag(2), tolerance = 1e-12, info = name)
      expect_equal(det(p$A), 1, tolerance = 1e-12, info = name)
    }
  }
})

test_that("the tree layouts hold each face once", {
  for (g in list(octa_gosper(), hex_grid(resolution = 1, aperture = 3))) {
    n_faces <- nrow(hexify:::icosa_solid(hexify:::icosa_arg(g))$faces)
    names <- if (n_faces == 8L) {
      c("plane", "gosper", "gosper_land", "land")
    } else {
      c("plane", "rhombic", "land")
    }
    for (name in names) {
      net <- net_layout(g, name)
      area <- tapply(vapply(net$pieces, piece_area, numeric(1)),
                     vapply(net$pieces, `[[`, numeric(1), "face"), sum)
      expect_equal(as.vector(area), rep(sqrt(3) / 4, n_faces), tolerance = 1e-12, info = name)
      # Every point of the sphere off the cuts has one place on the map.
      ll <- sphere_points(2000)
      xy <- net_project(net, ll[, 1], ll[, 2])
      expect_equal(tabulate(xy$point, nrow(ll)), rep(1L, nrow(ll)), info = name)
    }
  }
})

test_that("the Gosper World is four regular hexagons of one face and three thirds", {
  net <- net_layout(octa_gosper(), "gosper")
  expect_length(net$pieces, 16L)
  tiles <- split(net$pieces, vapply(net$pieces, `[[`, numeric(1), "tile"))
  expect_length(tiles, 4L)
  for (tile in tiles) {
    whole <- vapply(tile, function(p) all(startsWith(p$labels, "v")), logical(1))
    expect_equal(sum(whole), 1L)
    # The hexagon's corners: the middle face's three vertices and the centres
    # of its three neighbours, all 1 / sqrt(3) face edges from the middle
    # face's centre.
    mid <- tile[[which(whole)]]
    ctr <- colMeans(hexify:::place_points(mid$region, mid))
    corners <- unique(round(do.call(rbind, lapply(tile, function(p) {
      hexify:::place_points(p$region, p)
    })), 9))
    expect_equal(nrow(corners), 6L)
    r <- sqrt(rowSums(sweep(corners, 2, ctr)^2))
    expect_equal(r, rep(1 / sqrt(3), 6), tolerance = 1e-8)
  }
})

test_that("the Gosper tree has nine cuts and the flower three inside it", {
  g <- octa_gosper()
  e <- hexify:::net_edges(net_layout(g, "gosper"))
  expect_equal(sum(!e$joined), 18L)
  # Every cut runs from a vertex of the octahedron to a face centre.
  cut <- e[!e$joined, ]
  expect_true(all(xor(startsWith(cut$a, "v"), startsWith(cut$b, "v"))))
  expect_true(all(table(cut$key) == 2L))

  f <- hexify:::net_edges(net_layout(g, "gosper_flower"))
  # The flower's outline is 18 hexagon sides; three cuts inside it are
  # drawn from both sides.
  expect_equal(sum(!f$joined), 24L)
})

test_that("both sides of a cut are the same points of the sphere", {
  g <- octa_gosper()
  for (name in c("gosper", "gosper_flower", "gosper_land", "land")) {
    net <- net_layout(g, name)
    e <- hexify:::net_edges(net)
    cut <- e[!e$joined, ]
    for (key in unique(cut$key)) {
      sides <- cut[cut$key == key, ]
      lonlat <- lapply(seq_len(nrow(sides)), function(i) {
        s <- sides[i, ]
        p <- net$pieces[[s$piece]]
        # Run each side from its smaller label, so t = 0.3 is one point.
        ends <- rbind(c(s$x0, s$y0), c(s$x1, s$y1))
        if (s$a > s$b) ends <- ends[2:1, ]
        xy <- ends[1, ] + 0.3 * (ends[2, ] - ends[1, ])
        tri <- drop(solve(p$A, xy - p$b))
        hexify:::cpp_face_xy_to_ll(net$icosa, tri[1], tri[2], p$face)
      })
      for (ll in lonlat[-1]) {
        expect_equal(ll, lonlat[[1]], tolerance = 1e-9, info = paste(name, key))
      }
    }
  }
})

test_that("a cell crossing a joined side of the Gosper World continues across it", {
  g <- octa_gosper(3)
  net <- net_layout(g, "gosper")
  seg <- hexify:::net_segments(hexify:::grid_surface_paths(g, NULL, 0.02), net)
  e <- hexify:::net_edges(net)
  inner <- e[e$joined & !e$solid_edge, ]
  # Ends of cell segments on a joined side between two thirds of a face
  ends <- rbind(seg[, 1:2], seg[, 3:4])
  for (i in seq_len(nrow(inner))) {
    a <- c(inner$x0[i], inner$y0[i])
    d <- c(inner$x1[i], inner$y1[i]) - a
    off <- abs(d[1] * (ends[, 2] - a[2]) - d[2] * (ends[, 1] - a[1])) / sqrt(sum(d^2))
    along <- ((ends[, 1] - a[1]) * d[1] + (ends[, 2] - a[2]) * d[2]) / sum(d^2)
    on <- ends[off < 1e-9 & along > 1e-6 & along < 1 - 1e-6, , drop = FALSE]
    # Each crossing point is reached from both sides.
    k <- table(paste(round(on[, 1], 7), round(on[, 2], 7)))
    expect_true(all(k %% 2L == 0L), info = i)
  }
})

# =============================================================================
# Area and distortion
# =============================================================================

test_that("the Gosper World on Snyder's projection keeps areas to one constant", {
  g <- octa_gosper()
  net <- net_layout(g, "gosper")
  ll <- sphere_points(300)
  r <- pi / 180
  ratio <- vapply(seq_len(nrow(ll)), function(i) {
    # A small spherical triangle around the point
    lo <- ll[i, 1] + c(0, 0.01, 0) / cos(ll[i, 2] * r)
    la <- ll[i, 2] + c(0, 0, 0.01)
    if (any(abs(la) > 89.9)) return(NA_real_)
    xy <- net_project(net, lo, la)
    if (nrow(xy) != 3L || length(unique(xy$piece)) != 1L) return(NA_real_)
    plane <- abs(det(cbind(c(xy$x[2] - xy$x[1], xy$y[2] - xy$y[1]),
                           c(xy$x[3] - xy$x[1], xy$y[3] - xy$y[1])))) / 2
    V <- hexify:::unit_vec(lo, la)
    sphere <- 2 * atan(abs(sum(V[1, ] * hexify:::cross3(V[2, ], V[3, ]))) /
                         (1 + sum(V[1, ] * V[2, ]) + sum(V[2, ] * V[3, ]) +
                            sum(V[3, ] * V[1, ])))
    plane / sphere
  }, numeric(1))
  ratio <- ratio[!is.na(ratio)]
  expect_gt(length(ratio), 250L)
  # A face edge of 1 holds sqrt(3)/4 of map for 4 pi / 8 of sphere.
  expect_equal(ratio, rep((sqrt(3) / 4) / (pi / 2), length(ratio)), tolerance = 1e-4)
})

test_that("Snyder's projection reproduces his Table 1 distortion", {
  # Snyder (1992), Table 1: maximum angular deformation omega and the
  # largest and smallest scale factors a and b. On each face omega peaks on
  # the radii from the centre to the vertices, where the map has a cusp, so
  # the points sit just off one radius.
  table1 <- list(icosahedron = c(omega = 17.27, a = 1.163, b = 0.860),
                 octahedron = c(omega = 34.45, a = 1.357, b = 0.737))
  r <- pi / 180
  for (polyhedron in names(table1)) {
    icosa <- hexify:::projection_icosa("isea", polyhedron)
    fc <- hexify:::net_faces(icosa)[[1]]
    n_faces <- nrow(hexify:::icosa_solid(icosa)$faces)
    k <- sqrt((4 * pi / n_faces) / (sqrt(3) / 4))
    ctr <- colMeans(fc$tri)
    d <- fc$tri[1, ] - ctr
    perp <- c(-d[2], d[1]) / sqrt(sum(d^2))
    s <- seq(0.002, 0.998, by = 0.002)
    t <- cbind(ctr[1] + s * d[1] + 1e-5 * perp[1], ctr[2] + s * d[2] + 1e-5 * perp[2])
    ll <- t(vapply(seq_len(nrow(t)), function(i) {
      hexify:::cpp_face_xy_to_ll(icosa, t[i, 1], t[i, 2], 0L)
    }, numeric(2)))
    f <- function(lo, la) {
      hexify:::cpp_lonlat_to_face_solid(icosa, 0L, lo, la)[, c("tx", "ty"), drop = FALSE]
    }
    h <- 1e-7
    de <- (f(ll[, 1] + h, ll[, 2]) - f(ll[, 1] - h, ll[, 2])) / (2 * h * r * cos(ll[, 2] * r))
    dn <- (f(ll[, 1], ll[, 2] + h) - f(ll[, 1], ll[, 2] - h)) / (2 * h * r)
    sv <- t(vapply(seq_len(nrow(ll)), function(i) svd(k * cbind(de[i, ], dn[i, ]))$d,
                   numeric(2)))
    omega <- 2 * asin((sv[, 1] - sv[, 2]) / (sv[, 1] + sv[, 2])) / r
    expect_equal(sv[, 1] * sv[, 2], rep(1, nrow(sv)), tolerance = 1e-5, info = polyhedron)
    expect_equal(max(sv[, 1]), table1[[polyhedron]][["a"]], tolerance = 1e-3,
                 info = polyhedron)
    expect_equal(min(sv[, 2]), table1[[polyhedron]][["b"]], tolerance = 1e-3,
                 info = polyhedron)
    expect_lt(abs(max(omega) - table1[[polyhedron]][["omega"]]), 0.05)
  }
})

# =============================================================================
# The land-joined net
# =============================================================================

test_that("the land net is a spanning tree of whole faces without overlap", {
  grids <- list(hex_grid(resolution = 1, aperture = 3),
                hex_grid(resolution = 1, aperture = 3, orientation = "dymaxion",
                         projection = "fuller"),
                octa_gosper(1))
  for (g in grids) {
    net <- net_layout(g, "land")
    n <- length(net$pieces)
    e <- hexify:::net_edges(net)
    expect_equal(sum(e$joined) / 2, n - 1)
    polys <- sf::st_sfc(lapply(net$pieces, function(p) {
      xy <- hexify:::place_points(p$region, p)
      sf::st_polygon(list(rbind(xy, xy[1, ])))
    }))
    expect_equal(as.numeric(sf::st_area(sf::st_union(polys))), n * sqrt(3) / 4,
                 tolerance = 1e-9)
  }
})

test_that("no cut of the land net crosses more land than the joins it replaces", {
  # Cycle property of a maximum spanning tree: every edge left out carries
  # no more land than any edge on the tree path between its two faces.
  g <- hex_grid(resolution = 1, aperture = 3, orientation = "dymaxion",
                projection = "fuller")
  icosa <- hexify:::icosa_arg(g)
  E <- hexify:::icosa_solid(icosa)$edges
  w <- hexify:::edge_land_share(icosa, TRUE)
  e <- hexify:::net_edges(net_layout(g, "land"))
  key <- ifelse(paste0("v", E[, "v1"]) < paste0("v", E[, "v2"]),
                paste0("v", E[, "v1"], " v", E[, "v2"]),
                paste0("v", E[, "v2"], " v", E[, "v1"]))
  joined <- key %in% e$key[e$joined]
  expect_equal(sum(joined), 19L)
  tree_path <- function(from, to) {
    prev <- rep(NA_integer_, 20)
    prev[from] <- 0L
    queue <- from
    while (length(queue) > 0L) {
      f <- queue[1]
      queue <- queue[-1]
      for (j in which(joined & (E[, "f1"] == f | E[, "f2"] == f))) {
        n <- if (E[j, "f1"] == f) E[j, "f2"] else E[j, "f1"]
        if (is.na(prev[n])) {
          prev[n] <- j
          queue <- c(queue, n)
        }
      }
    }
    path <- integer(0)
    f <- to
    while (f != from) {
      j <- prev[f]
      path <- c(path, j)
      f <- if (E[j, "f1"] == f) E[j, "f2"] else E[j, "f1"]
    }
    path
  }
  for (j in which(!joined)) {
    expect_true(all(w[tree_path(E[j, "f1"], E[j, "f2"])] >= w[j]), info = j)
  }
  expect_gt(sum(w[joined]), 0)
})

test_that("gosper_land joins the four tiles along the most land", {
  g <- octa_gosper()
  net <- net_layout(g, "gosper_land")
  e <- hexify:::net_edges(net)
  tile <- vapply(net$pieces, `[[`, numeric(1), "tile")
  # The twelve hexagon sides run from a vertex to a face centre; each lies
  # between two tiles.
  hs <- e[startsWith(e$a, "c") != startsWith(e$b, "c"), ]
  sides <- unique(hs$key)
  expect_length(sides, 12L)
  ends <- lapply(sides, function(k) sort(unique(tile[hs$piece[hs$key == k]])))
  expect_true(all(lengths(ends) == 2L))
  solid <- hexify:::icosa_solid(net$icosa)
  num <- function(k, p) {
    as.integer(sub(p, "", regmatches(k, regexpr(paste0(p, "[0-9]+"), k))))
  }
  w <- hexify:::arc_land_share(solid$vertices[num(sides, "v"), , drop = FALSE],
                               solid$normals[num(sides, "c") + 1L, , drop = FALSE], TRUE)
  joined <- sides %in% hs$key[hs$joined]
  expect_equal(sum(!e$joined), 2L * sum(!joined))
  tree_weight <- function(pool) {
    sets <- utils::combn(which(pool), 3L)
    max(apply(sets, 2, function(j) {
      links <- do.call(rbind, ends[j])
      reach <- links[1, ]
      for (i in 1:3) {
        reach <- unique(c(reach, links[apply(links, 1, function(l) any(l %in% reach)), ]))
      }
      if (length(reach) == 4L) sum(w[j]) else -Inf
    }))
  }
  # The joined sides hold a tree of the four tiles; where it joins two of the
  # three sides around a face centre, the third closes by itself. No choice
  # of three sides that joins the tiles carries more land than that tree.
  expect_equal(tree_weight(joined), tree_weight(rep(TRUE, 12L)))
})
test_that("the rhombic layout is DGGAL's 5 x 6 space with y down the page", {
  g <- hex_grid(resolution = 2, aperture = 3)
  net <- net_layout(g, "rhombic")
  for (p in net$pieces) {
    d <- matrix(hexify:::RHOMBIC_5X6[p$face + 1L, ], 3L, 2L, byrow = TRUE)
    expect_equal(unname(hexify:::place_points(p$region, p)), cbind(d[, 1], -d[, 2]),
                 tolerance = 1e-12)
    # A unit-edge triangle onto half a unit square, the same shear on every face
    expect_equal(det(p$A), 2 / sqrt(3), tolerance = 1e-12)
  }
  e <- hexify:::net_edges(net)
  # One staircase: 19 joined face edges, its outline the 11 cut ones twice
  expect_equal(sum(e$joined) / 2, 19)
  expect_equal(sum(!e$joined), 22L)
  ll <- sphere_points(1000)
  xy <- net_project(net, ll[, 1], ll[, 2])
  expect_true(all(xy$x >= -1e-9 & xy$x <= 5 + 1e-9 & xy$y <= 1e-9 & xy$y >= -6 - 1e-9))
  expect_error(net_layout(octa_gosper(), "rhombic"), "icosahedron")
})

# =============================================================================
# net_project() and arguments
# =============================================================================

test_that("net_project places a point where the plot draws its face", {
  g <- octa_gosper()
  net <- net_layout(g, "gosper_flower")
  xy <- net_project(net, c(16.37, -74.0, 0), c(48.21, 40.71, 90))
  expect_named(xy, c("point", "piece", "x", "y"))
  # The middle tile holds the north pole, placed at the origin.
  pole <- xy[xy$point == 3L, ]
  expect_equal(c(pole$x[1], pole$y[1]), c(0, 0), tolerance = 1e-9)
  expect_equal(nrow(net_project(net, numeric(0), numeric(0))), 0L)
})

test_that("the Gosper layouts lie flat-topped with north near up at the centre", {
  g <- octa_gosper()
  for (centre in list(c(10, 50), c(-120, -30), c(0, 90))) {
    net <- net_layout(g, "gosper", centre = centre)
    north <- if (centre[2] > 89) c(centre[1], 89.99) else centre + c(0, 0.01)
    xy <- net_project(net, c(centre[1], north[1]), c(centre[2], north[2]))
    a <- xy[xy$point == 1L, ][1, ]
    b <- xy[xy$point == 2L, ][1, ]
    expect_equal(c(a$x, a$y), c(0, 0), tolerance = 1e-9)
    # North, or at the pole the meridian of `lon` pointing down
    up <- atan2(b$x - a$x, b$y - a$y) * 180 / pi
    if (centre[2] > 89) up <- (up + 360) %% 360 - 180
    expect_lte(abs(up), 30 + 1e-9)
    # Every hexagon side is horizontal or at 60 degrees to it.
    e <- hexify:::net_edges(net)
    side <- e[!e$joined, ]
    angle <- (atan2(side$y1 - side$y0, side$x1 - side$x0) * 180 / pi) %% 60
    expect_lt(max(pmin(angle, 60 - angle)), 1e-9)
  }
})

test_that("mirror builds the hexagons on the other four faces", {
  g <- octa_gosper()
  mid <- function(net) {
    sort(vapply(Filter(function(p) all(startsWith(p$labels, "v")), net$pieces),
                `[[`, numeric(1), "face"))
  }
  a <- mid(net_layout(g, "gosper"))
  b <- mid(net_layout(g, "gosper", mirror = TRUE))
  expect_length(a, 4L)
  expect_setequal(c(a, b), 0:7)
})

test_that("net_layout checks its arguments", {
  ico <- hex_grid(resolution = 1, aperture = 3)
  expect_error(net_layout(ico, "gosper"), "octahedron")
  expect_error(net_layout(hex_grid(resolution = 0, type = "h3")), "ISEA grid")
  expect_error(net_layout(ico, "plane", centre = c(0, 0)), "Gosper layout")
  expect_error(net_layout(ico, "plane", mirror = TRUE), "Gosper layout")
  expect_error(net_layout(ico, "plane", land = FALSE), "keep joined")
  expect_error(net_layout(octa_gosper(), "gosper", centre = c(0, 100)), "c\\(lon, lat\\)")
  expect_error(net_project(list(), 0, 0), "hexify_net")
  expect_output(print(net_layout(octa_gosper(), "gosper")), "16 pieces in 4 tiles")
})

# =============================================================================
# The "gosper" orientation
# =============================================================================

test_that("orientation = \"gosper\" puts the poles at octahedron edge midpoints", {
  g <- octa_gosper()
  V <- hexify:::icosa_solid(hexify:::icosa_arg(g))$vertices
  # The poles are 45 degrees from the two ends of an edge and 90 from the rest.
  pole_angles <- sort(acos(pmin(1, V %*% c(0, 0, 1))) * 180 / pi)
  expect_equal(pole_angles, c(45, 45, 90, 90, 135, 135), tolerance = 1e-9)

  ll <- hexify:::vec_lonlat(V)
  pts <- sf::st_as_sfc(lapply(seq_len(nrow(ll)), function(i) sf::st_point(ll[i, ])),
                       crs = 4326)
  old <- suppressMessages(sf::sf_use_s2(TRUE))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  land <- sf::st_union(sf::st_make_valid(sf::st_geometry(hexify_world)))
  expect_gt(min(as.numeric(sf::st_distance(pts, land))), 600e3)

  expect_error(hex_grid(resolution = 1, orientation = "gosper"), "octahedron")
})
