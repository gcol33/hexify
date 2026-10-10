# Grids on the octahedron and the tetrahedron, and Snyder's projection on
# both.

octa_grids <- list(
  list(ap = 3, res = 3), list(ap = 4, res = 3), list(ap = 7, res = 2),
  list(ap = 7, res = 3), list(ap = "4/3", res = 3), list(ap = c(4, 7, 3), res = 3)
)

# The tetrahedron's second quad holds the far edge of its box in place of the
# near edge it shares with the first; resolution 1 puts a single cell on each
# of those edges, and resolution 4 of aperture 3 a Class I grid with many.
tetra_grids <- c(octa_grids, list(list(ap = 4, res = 1), list(ap = 3, res = 1),
                                  list(ap = 7, res = 1), list(ap = 3, res = 4)))

# Each solid with its grids, its vertex count and the sides of a vertex cell
small_solids <- list(
  octahedron = list(grids = octa_grids, n_verts = 6, sides = 4),
  tetrahedron = list(grids = tetra_grids, n_verts = 4, sides = 3)
)

solid_grid <- function(s, solid) {
  hex_grid(resolution = s$res, aperture = s$ap, polyhedron = solid)
}

octa_grid <- function(s) solid_grid(s, "octahedron")

test_that("the solids carry Snyder's constants and their topology", {
  ico <- hexify:::solid_info("icosahedron")
  oct <- hexify:::solid_info("octahedron")
  tet <- hexify:::solid_info("tetrahedron")
  expect_equal(c(ico$n_faces, ico$n_verts, ico$n_diamonds), c(20, 12, 10))
  expect_equal(c(oct$n_faces, oct$n_verts, oct$n_diamonds), c(8, 6, 4))
  expect_equal(c(tet$n_faces, tet$n_verts, tet$n_diamonds), c(4, 4, 2))
  expect_true(ico$has_grid)
  expect_true(oct$has_grid)
  expect_true(tet$has_grid)
  expect_equal(oct$valence, rep(4, 6))
  expect_equal(tet$valence, rep(3, 4))
  # Snyder (1992), Table 1
  expect_equal(ico$g_deg, 37.37736814, tolerance = 1e-9)
  expect_equal(oct$g_deg, 54.73561032, tolerance = 1e-9)
  expect_equal(tet$g_deg, 70.52877937, tolerance = 1e-9)
  expect_equal(c(ico$G_deg, oct$G_deg, tet$G_deg), c(36, 45, 60))
})

test_that("hex_grid() builds octahedral and tetrahedral grids", {
  g <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron")
  expect_identical(g@polyhedron, "octahedron")
  expect_identical(hex_grid(resolution = 3)@polyhedron, "icosahedron")
  expect_output(print(g), "octahedron")
  expect_equal(n_cells(g), 4 * 4^3 + 2)
  t4 <- hex_grid(resolution = 3, aperture = 4, polyhedron = "tetrahedron")
  expect_identical(t4@polyhedron, "tetrahedron")
  expect_output(print(t4), "tetrahedron")
  for (ap in c(3, 4, 7)) {
    expect_equal(n_cells(hex_grid(resolution = 3, aperture = ap,
                                  polyhedron = "tetrahedron")), 2 * ap^3 + 2)
  }
  expect_error(hex_grid(resolution = 3, polyhedron = "octahedron", projection = "fuller"),
               "icosahedron")
  expect_error(hex_grid(resolution = 3, polyhedron = "tetrahedron", projection = "fuller"),
               "icosahedron")
  expect_error(hex_grid(resolution = 3, type = "h3", polyhedron = "octahedron"))
})

test_that("octahedral and tetrahedral cells tile the sphere with equal areas", {
  for (solid in names(small_solids)) {
    sp <- small_solids[[solid]]
    for (s in sp$grids) {
      g <- solid_grid(s, solid)
      ids <- seq_len(n_cells(g))
      sides <- hexify:::isea_cell_sides(ids, g)
      expect_equal(sum(sides == sp$sides), sp$n_verts)
      expect_true(all(sides[sides != sp$sides] == 6))
      area <- cell_area(ids, g)
      expect_equal(sum(area), hexify:::body_surface_km2(hexify:::grid_radius_km(g)),
                   tolerance = 1e-12)
      lv <- hexify:::isea_levels(g@aperture, g@resolution)
      sa <- hexify:::cpp_cell_solid_angle(hexify:::icosa_arg(g), as_cell_id(ids),
                                           lv$resolution, lv$aperture, lv$ap_seq, 1e-4)
      hex <- 4 * pi / (length(ids) - 2)
      expect_equal(sa[sides == 6] / hex, rep(1, sum(sides == 6)), tolerance = 1e-5)
      expect_equal(sa[sides == sp$sides] / hex, rep(sp$sides / 6, sp$n_verts),
                   tolerance = 1e-5)
    }
  }
})

test_that("every octahedral and tetrahedral cell corner is shared by three cells", {
  for (solid in names(small_solids)) {
    for (s in small_solids[[solid]]$grids) {
      g <- solid_grid(s, solid)
      ids <- seq_len(n_cells(g))
      rings <- hexify:::isea_cell_rings(ids, g@resolution, g@aperture,
                                        hexify:::icosa_arg(g), tolerance = 0)
      # A corner on a pole is written twice, once per meridian reaching it
      key <- unlist(lapply(rings, function(r) {
        v <- hexify:::unit_vec(r[-nrow(r), 1], r[-nrow(r), 2])
        unique(paste(round(v[, 1], 7), round(v[, 2], 7), round(v[, 3], 7)))
      }))
      expect_true(all(table(key) == 3))
    }
  }
})

test_that("octahedral and tetrahedral cells round-trip, their neighbours symmetric", {
  for (solid in names(small_solids)) {
    for (s in small_solids[[solid]]$grids) {
      g <- solid_grid(s, solid)
      ids <- seq_len(n_cells(g))
      ctr <- cell_to_lonlat(ids, g)
      expect_equal(lonlat_to_cell(ctr$lon_deg, ctr$lat_deg, g), as_cell_id(ids))
      sides <- hexify:::isea_cell_sides(ids, g)
      nb <- get_neighbors(ids, g)
      expect_equal(lengths(nb), sides)
      for (k in ids) {
        expect_true(all(vapply(as.integer(nb[[k]]), function(m) k %in% nb[[m]],
                               logical(1))))
      }
      # Each wall of a cell faces one neighbour
      lv <- hexify:::isea_levels(g@aperture, g@resolution)
      w <- hexify:::cpp_cell_walls(hexify:::icosa_arg(g), as_cell_id(ids), lv$resolution,
                                   lv$aperture, lv$ap_seq, 1e-6, TRUE)
      expect_equal(nrow(w$walls), sum(sides))
    }
  }
})

test_that("tetrahedral cells hold the points of their share of the sphere", {
  # A uniform sample falls in each cell as often as its area says: a hexagon
  # 1 / (N - 2) of the sphere, a vertex triangle half of that.
  set.seed(11)
  n <- 60000
  lon <- runif(n, -180, 180)
  lat <- asin(runif(n, -1, 1)) * 180 / pi
  for (s in list(list(ap = 4, res = 2), list(ap = 3, res = 3), list(ap = 7, res = 2))) {
    g <- solid_grid(s, "tetrahedron")
    ids <- seq_len(n_cells(g))
    count <- tabulate(as.integer(lonlat_to_cell(lon, lat, g)), length(ids))
    sides <- hexify:::isea_cell_sides(ids, g)
    expected <- n / (length(ids) - 2) * ifelse(sides == 3, 0.5, 1)
    expect_lt(max(abs(count - expected) / sqrt(expected)), 4.5)
  }
})

test_that("octahedral and tetrahedral cells keep their parents and children", {
  for (solid in names(small_solids)) {
    for (ap in list(3, 4, 7, "4/3")) {
      g <- hex_grid(resolution = 3, aperture = ap, polyhedron = solid)
      ids <- seq_len(n_cells(g))
      parent <- get_parent(ids, g)
      gp <- hexify:::grid_at_resolution(g, 2)
      expect_identical(gp@polyhedron, solid)
      expect_true(all(parent >= 1 & parent <= n_cells(gp)))
      kids <- get_children(seq_len(n_cells(gp)), gp)
      expect_true(all(vapply(seq_along(kids),
                             function(p) all(get_parent(kids[[p]], g) == p), logical(1))))
      expect_setequal(as.integer(hexify:::cell_id_unlist(kids)), ids)
      idx <- cell_to_index(ids, g)
      expect_false(anyDuplicated(idx) > 0)
    }
    for (ap in c(3, 4, 7)) {
      g <- hex_grid(resolution = 3, aperture = ap, polyhedron = solid)
      ids <- seq_len(n_cells(g))
      it <- hexify:::index_type_for_aperture(g@aperture)
      a <- hexify:::aperture_to_int(g@aperture)
      expect_equal(hexify:::isea_index_to_cells(cell_to_index(ids, g), a, it,
                                                hexify:::icosa_arg(g)), ids)
    }
  }
})

test_that("a tetrahedral cell's parent holds its centre or borders the cell that does", {
  for (ap in c(3, 4, 7)) {
    g <- hex_grid(resolution = 4, aperture = ap, polyhedron = "tetrahedron")
    gp <- hexify:::grid_at_resolution(g, 3)
    ids <- seq_len(n_cells(g))
    parent <- get_parent(ids, g)
    ctr <- cell_to_lonlat(ids, g)
    holder <- lonlat_to_cell(ctr$lon_deg, ctr$lat_deg, gp)
    nb <- get_neighbors(holder, gp)
    near <- parent == holder |
      vapply(seq_along(ids), function(k) parent[k] %in% nb[[k]], logical(1))
    expect_true(all(near))
    # Aperture 7 nests: the Z7 parent is the cell holding the centre
    if (ap == 7) expect_equal(parent, holder)
  }
})

test_that("cells on the tetrahedron's folded edges carry their own index", {
  # Quad 2 holds the far edge of its box: those cells are written under 06 at
  # their place on the near edge, and read back to themselves.
  g <- hex_grid(resolution = 3, aperture = 4, polyhedron = "tetrahedron")
  ids <- seq_len(n_cells(g))
  q <- hexify:::grid_quad_ij(ids, g)
  far <- q$quad == 2 & q$j == 8
  expect_equal(sum(far), 7)
  idx <- cell_to_index(ids, g)
  expect_true(all(substr(idx[far], 1, 2) == "06"))
  expect_false(any(substr(idx[!far], 1, 2) == "06"))
  expect_equal(hexify:::isea_index_to_cells(idx[far], 4L, "zorder", hexify:::icosa_arg(g)),
               as_cell_id(ids[far]))
})

test_that("octahedral and tetrahedral vertex cells compact like any cell", {
  for (solid in names(small_solids)) {
    for (ap in c(3, 4, 7)) {
      g <- hex_grid(resolution = 3, aperture = ap, polyhedron = solid)
      ids <- seq_len(n_cells(g))
      expect_false(any(is_pentagon(ids, g)))
      idx <- cell_to_index(ids, g)
      packed <- hex_compact(idx, g)
      expect_lt(length(packed), length(idx))
      expect_setequal(hex_uncompact(packed, g, 3), idx)
    }
  }
  g <- hex_grid(resolution = 3, aperture = 4)
  expect_equal(sum(is_pentagon(seq_len(n_cells(g)), g)), 12)
})

test_that("Snyder's projection inverts on the octahedron and the tetrahedron", {
  set.seed(7)
  lon <- runif(500, -180, 180)
  lat <- asin(runif(500, -1, 1)) * 180 / pi
  for (solid in c("octahedron", "tetrahedron")) {
    n_faces <- nrow(hexify_face_centers(solid))
    expect_equal(n_faces, c(octahedron = 8, tetrahedron = 4)[[solid]])
    faces <- integer(0)
    for (k in seq_along(lon)) {
      f <- hexify_forward(lon[k], lat[k], polyhedron = solid)
      faces <- c(faces, f[["face"]])
      expect_equal(hexify_which_face(lon[k], lat[k], polyhedron = solid), f[["face"]])
      ll <- hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]], f[["face"]],
                           polyhedron = solid)
      gap <- row_angle(hexify:::unit_vec(ll[[1]], ll[[2]]), hexify:::unit_vec(lon[k], lat[k]))
      expect_lt(gap, 1e-11)
    }
    expect_setequal(unique(faces), seq_len(n_faces) - 1)
    expect_error(hexify_forward(0, 0, projection = "fuller", polyhedron = solid))
  }
})

test_that("Snyder's projection on the tetrahedron is equal-area", {
  # Equal areas on the sphere map to equal areas on the face: a face holds a
  # quarter of the sphere, so a uniform sample falls in each face equally often
  # and spreads uniformly over the face triangle, three quarters of whose area
  # lies below half its height.
  set.seed(3)
  n <- 20000
  lon <- runif(n, -180, 180)
  lat <- asin(runif(n, -1, 1)) * 180 / pi
  f <- t(vapply(seq_len(n), function(k) hexify_forward(lon[k], lat[k], polyhedron = "tetrahedron"),
                numeric(3)))
  expect_equal(as.numeric(table(f[, 1])) / n, rep(0.25, 4), tolerance = 0.05)
  expect_equal(mean(f[, 3] < sqrt(3) / 4), 0.75, tolerance = 0.02)
})

test_that("the octahedron has 6 unit vertices, 8 faces and 12 equal edges", {
  s <- hexify:::icosa_solid(hexify:::standard_icosa("octahedron"))
  expect_equal(dim(s$vertices), c(6L, 3L))
  expect_equal(rowSums(s$vertices^2), rep(1, 6), tolerance = 1e-12)
  expect_equal(dim(s$faces), c(8L, 3L))
  expect_equal(nrow(s$edges), 12L)
  len <- sqrt(rowSums((s$vertices[s$edges[, "v1"], ] - s$vertices[s$edges[, "v2"], ])^2))
  expect_equal(len, rep(sqrt(2), 12), tolerance = 1e-12)
})

test_that("octahedral and tetrahedral cell boundaries lie on the faces and on the sphere", {
  for (solid in names(small_solids)) {
    s <- hexify:::icosa_solid(hexify:::standard_icosa(solid))
    inradius <- sqrt(sum(colMeans(s$vertices[s$faces[1, ], ])^2))
    for (ap in c(3, 4, 7)) {
      g <- hex_grid(resolution = 2, aperture = ap, polyhedron = solid)
      P <- hexify:::grid_surface_paths(g, NULL, 0.02)
      Q <- P[, c("solid_x", "solid_y", "solid_z")]
      S <- P[, c("sphere_x", "sphere_y", "sphere_z")]
      expect_equal(rowSums(Q * s$normals[P[, "face"] + 1, ]), rep(inradius, nrow(P)),
                   tolerance = 1e-7)
      expect_equal(rowSums(S^2), rep(1, nrow(P)), tolerance = 1e-12)
      expect_setequal(unique(P[, "cell"]), seq_len(n_cells(g)))
    }
  }
})

test_that("octahedral and tetrahedral grids plot on every surface and on the globe", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (solid in names(small_solids)) {
    g <- hex_grid(resolution = 2, aperture = 4, polyhedron = solid)
    for (surface in c("sphere", "solid", "net")) {
      expect_identical(plot(g, surface = surface, land = FALSE), g)
    }
  }
  skip_if_not_installed("htmlwidgets")
  for (solid in names(small_solids)) {
    g <- hex_grid(resolution = 2, aperture = 4, polyhedron = solid)
    w <- hex_globe(g, values = seq_len(n_cells(g)), surface = "solid", land = FALSE)
    # The Grid uniform names the solid's faces, and its levels hold every cell
    uniform <- readBin(hexify:::cpp_base64_decode(w$x$grid$uniform), "integer", n = 1296,
                       size = 4, endian = "little")
    expect_equal(uniform[16], c(octahedron = 8, tetrahedron = 4)[[solid]])
    expect_equal(w$x$grid$given[1], as.numeric(n_cells(g)))
  }
})

test_that("the tetrahedron and the octahedron unfold into a net", {
  for (solid in c("octahedron", "tetrahedron")) {
    p <- hexify_lonlat_to_plane(c(10, -120, 75), c(20, -45, 60), polyhedron = solid)
    expect_true(all(is.finite(unlist(p[c("plane_x", "plane_y")]))))
  }
})
