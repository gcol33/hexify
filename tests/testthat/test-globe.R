# tests/testthat/test-globe.R
# Meshes and paths on the icosahedron and the sphere behind hex_globe(), and
# the widget itself.

mesh_points <- function(v) matrix(v, ncol = 3, byrow = TRUE)

# The bytes of base64 text.
b64_bytes <- function(s) {
  v <- match(strsplit(gsub("=", "", s), "")[[1]],
             c(LETTERS, letters, 0:9, "+", "/")) - 1L
  bits <- vapply(v, function(n) as.integer(intToBits(n))[6:1], integer(6))
  bits <- as.vector(bits)
  bits <- matrix(bits[seq_len(length(bits) %/% 8 * 8)], nrow = 8)
  packBits(as.logical(bits[8:1, ]), "raw")
}

# Areas of the sphere triangles of a mesh (Van Oosterom and Strackee), summed
# per item.
mesh_item_area <- function(m) {
  P <- mesh_points(m$sphere)
  i <- matrix(m$index + 1L, ncol = 3, byrow = TRUE)
  A <- P[i[, 1], , drop = FALSE]
  B <- P[i[, 2], , drop = FALSE]
  C <- P[i[, 3], , drop = FALSE]
  BxC <- cbind(B[, 2] * C[, 3] - B[, 3] * C[, 2],
               B[, 3] * C[, 1] - B[, 1] * C[, 3],
               B[, 1] * C[, 2] - B[, 2] * C[, 1])
  a <- 2 * atan2(abs(rowSums(A * BxC)),
                 1 + rowSums(A * B) + rowSums(B * C) + rowSums(C * A))
  tapply(a, m$item[i[, 1]], sum)
}

# =============================================================================
# Meshes
# =============================================================================

test_that("the faces mesh covers the sphere once, on the unit sphere", {
  m <- hexify:::cpp_globe_faces(numeric(0), 0.1)
  expect_equal(sum(mesh_item_area(m)), 4 * pi, tolerance = 1e-9)
  expect_equal(sqrt(rowSums(mesh_points(m$sphere)^2)),
               rep(1, length(m$item)), tolerance = 1e-12)
  expect_true(all(sqrt(rowSums(mesh_points(m$solid)^2)) <= 1 + 1e-12))
  expect_equal(sort(unique(m$item)), 1:20)
})

test_that("no edge of a refined mesh is longer than the spacing", {
  m <- hexify:::cpp_globe_faces(numeric(0), 0.1)
  S <- mesh_points(m$solid)
  i <- matrix(m$index + 1L, ncol = 3, byrow = TRUE)
  # A face edge of the inscribed icosahedron is 1.0515 radii long.
  edge <- sqrt(rowSums((S[i[, 1], ] - S[i[, 2], ])^2))
  expect_lte(max(edge), 0.1 * 1.0515 + 1e-9)
})

test_that("H3 cells tile the sphere", {
  cells <- h3_all_cells(0)
  area <- mesh_item_area(hexify:::h3_surface_mesh(cells, 0.05))
  expect_length(area, length(cells))
  expect_equal(sum(area), 4 * pi, tolerance = 1e-6)
})

test_that("the faces mesh carries each vertex's triangle coordinates", {
  m <- hexify:::cpp_globe_faces(numeric(0), 0.25)
  tri <- matrix(m$tri, ncol = 2, byrow = TRUE)
  face <- m$item - 1L
  for (k in seq(1, nrow(tri), by = 7)) {
    expect_equal(as.vector(hexify:::cpp_face_tri_to_solid(numeric(0), face[k], tri[k, 1], tri[k, 2])),
                 mesh_points(m$solid)[k, ], tolerance = 1e-12)
  }
})

test_that("land triangles keep the land's area", {
  land <- hexify:::surface_land(TRUE)
  polys <- hexify:::sfc_polygons(land)
  m <- hexify:::cpp_globe_polygons(numeric(0), polys, 0.05)
  area <- mesh_item_area(m)
  expect_length(area, length(polys))
  ref <- vapply(polys, function(p) {
    as.numeric(sf::st_area(sf::st_sfc(sf::st_polygon(p), crs = 4326)))
  }, numeric(1))
  # s2 measures on a sphere of 6371.0088 km.
  expect_equal(as.vector(area) * 6371008.8^2, ref[as.integer(names(area))],
               tolerance = 1e-5)
})

test_that("a polygon with a hole leaves the hole open", {
  outer <- cbind(c(0, 20, 20, 0, 0), c(0, 0, 20, 20, 0))
  hole <- cbind(c(5, 5, 15, 15, 5), c(5, 15, 15, 5, 5))
  full <- mesh_item_area(hexify:::cpp_globe_polygons(numeric(0), list(list(outer)), 0.05))
  holed <- mesh_item_area(hexify:::cpp_globe_polygons(numeric(0), list(list(outer, hole)), 0.05))
  inner <- mesh_item_area(hexify:::cpp_globe_polygons(numeric(0), list(list(hole)), 0.05))
  expect_equal(unname(holed), unname(full - inner), tolerance = 1e-9)
})

test_that("a polygon wider than a hemisphere is refused", {
  ring <- cbind(c(-170, 0, 170, 0, -170), c(0, -80, 0, 80, 0))
  expect_error(hexify:::cpp_globe_polygons(numeric(0), list(list(ring)), 0.05),
               "88 degrees")
})

# =============================================================================
# Paths
# =============================================================================

test_that("paths on the faces keep both surfaces and stay on their faces", {
  m <- hexify:::cpp_sphere_paths_on_faces(numeric(0), c(0, 90, 179, -170), c(0, 10, 20, 60),
                                          rep(1L, 4), 0.02)
  expect_equal(sqrt(rowSums(m[, c("sphere_x", "sphere_y", "sphere_z")]^2)),
               rep(1, nrow(m)), tolerance = 1e-12)
  solid <- hexify:::icosa_solid(numeric(0))
  S <- m[, c("solid_x", "solid_y", "solid_z")]
  n <- solid$normals[m[, "face"] + 1L, ]
  h <- rowSums(solid$normals * solid$vertices[solid$faces[, 1], ])[m[, "face"] + 1L]
  expect_equal(rowSums(S * n), h, tolerance = 1e-9)
})

test_that("face edges lie on the edges of the solid", {
  m <- hexify:::edge_surface_paths(0.05, numeric(0))
  solid <- hexify:::icosa_solid(numeric(0))
  V <- solid$vertices
  E <- solid$edges
  S <- m[, c("solid_x", "solid_y", "solid_z")]
  # Distance of each point to the segment of its own edge.
  k <- m[, "cell"]
  a <- V[E[k, "v1"], ]
  b <- V[E[k, "v2"], ]
  t <- pmin(pmax(rowSums((S - a) * (b - a)) / rowSums((b - a)^2), 0), 1)
  expect_lt(max(sqrt(rowSums((S - (a + t * (b - a)))^2))), 1e-9)
})

# =============================================================================
# Widget
# =============================================================================

test_that("buffers are base64 of little-endian 32-bit words", {
  expect_equal(hexify:::cpp_base64_buffer(c(1, 2.5), "f32"), "AACAPwAAIEA=")
  expect_equal(hexify:::cpp_base64_buffer(c(1, 7), "u32"), "AQAAAAcAAAA=")
  expect_equal(hexify:::cpp_base64_buffer(numeric(0), "f32"), "")
  expect_error(hexify:::cpp_base64_buffer(1, "f64"), "type")
})

test_that("values map onto the ramp, NA below it", {
  pos <- hexify:::ramp_position(c(0, 5, 10, NA, 20), 1:5, c(0, 10))
  expect_equal(pos, c(0, 0.5, 1, -1, 1))
  expect_equal(hexify:::ramp_position(c(3, 3), 1:2, NULL), c(0.5, 0.5))
  expect_error(hexify:::ramp_position(1:3, 1:2, NULL), "one per cell")
})

test_that("a palette is a named hcl palette or several colours", {
  expect_length(hexify:::ramp_colours("viridis"), 256)
  expect_length(hexify:::ramp_colours("Blue-Red 3"), 256)
  expect_length(hexify:::ramp_colours(c("white", "#FF000080")), 256)
  expect_error(hexify:::ramp_colours("not a palette"), "palette")
})

test_that("hex_globe builds a widget with every layer", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 2, aperture = 3)
  ids <- seq_len(grid_n_cells(g))
  w <- hex_globe(g, values = ids, surface = "icosahedron",
                 projection = "perspective", distance = 2, tilt = 20)
  expect_s3_class(w, "htmlwidget")
  x <- w$x
  expect_true(nzchar(x$shader))
  expect_equal(x$fold, 0)
  expect_true(x$foldable)
  for (layer in c("surface", "land")) {
    expect_true(x[[layer]]$n_index > 0, info = layer)
  }
  expect_true(nzchar(x$surface$tri))
  for (layer in c("land_lines", "edge_lines")) {
    expect_true(x[[layer]]$n_segment > 0, info = layer)
  }
  expect_equal(x$grid$n_keys, 0)
  expect_true(x$grid$dense)
  expect_true(x$grid$all)
  expect_equal(x$grid$ramp_map, c(1, 1 / (length(ids) - 1), 0))
  expect_length(x$projection$faces, 20 * 16)
  expect_length(x$projection$edges, 12 * 8)
  expect_null(x$cells)
  expect_null(x$grid_lines)
  expect_equal(x$camera$distance, 2)
  expect_length(x$palette, 4 * 256)
})

test_that("an ISEA grid is sent as its frame and its cells' IDs, not outlines", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 10, aperture = 3)
  x <- hex_globe(g, land = FALSE)$x
  expect_equal(x$grid$n_keys, 0)
  expect_null(x$grid$keys)
  expect_true(x$grid$all)
  expect_equal(x$grid$per_quad, c(0, 3^10))
  expect_lt(object.size(x), 2e6)

  cells <- c(5e5, 17, 3, 590492)
  x <- hex_globe(g, values = c(1, NA, 3, 4), cells = cells, land = FALSE)$x
  expect_false(x$grid$all)
  expect_equal(x$grid$n_keys, 4)
  expect_false(x$grid$dense)
  # Sorted IDs as high and low words, values in the same order.
  words <- readBin(b64_bytes(x$grid$keys), "integer", n = 8, size = 4,
                   endian = "little")
  expect_equal(words, c(0, 3, 0, 17, 0, 5e5, 0, 590492))
  raw <- readBin(b64_bytes(x$grid$values), "numeric", n = 4, size = 4,
                 endian = "little")
  expect_equal(raw[-2], c(3, 1, 4))
  expect_true(is.nan(raw[2]))
  expect_error(hex_globe(g, cells = c(1, 590493)), "cell IDs")
})

test_that("a grid finer than 32-bit floats resolve is refused", {
  skip_if_not_installed("htmlwidgets")
  expect_silent(hexify:::globe_grid(hex_grid(resolution = 16, aperture = 7),
                                    NULL, NULL, NULL))
  expect_error(hex_globe(hex_grid(resolution = 17, aperture = 7)), "32-bit")
  expect_error(hex_globe(hex_grid(resolution = 24, aperture = 4)), "32-bit")
})

test_that("values map onto the ramp the shader reads", {
  expect_equal(hexify:::ramp_map(c(2, NA, 6), NULL), c(2, 0.25, 0))
  expect_equal(hexify:::ramp_map(c(3, 3), NULL), c(3, 0, 0.5))
  expect_equal(hexify:::ramp_map(NA_real_, NULL), c(0, 1, 0))
  expect_error(hexify:::ramp_map(1, c(1, NA)), "two numbers")
})

test_that("cell IDs above 2^32 split into two words", {
  expect_equal(hexify:::split_u64(c(5, 2^32 + 7, 3 * 2^40)),
               c(0, 5, 1, 7, 3 * 2^8, 0))
})

test_that("the frame names the substrate sublattice of each grid", {
  f <- hexify:::cpp_globe_frame(3L, 3L, integer(0))
  expect_equal(f[c("dim", "index", "c")], list(dim = 9, index = 3, c = 2))
  expect_equal(f$generator, c(2, 1))
  f <- hexify:::cpp_globe_frame(3L, 7L, integer(0))
  expect_equal(f[c("dim", "index", "c")], list(dim = 49, index = 7, c = 5))
  f <- hexify:::cpp_globe_frame(4L, 4L, integer(0))
  expect_equal(f[c("dim", "index", "per_quad")], list(dim = 16, index = 1, per_quad = 256))
  f <- hexify:::cpp_globe_frame(0L, 0L, c(4L, 3L, 7L))
  expect_equal(f$index, 21)
  expect_equal(f$per_quad, 21)
})

test_that("hex_globe leaves out layers that are switched off", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 1, aperture = 4)
  x <- hex_globe(g, land = FALSE, face_edges = FALSE)$x
  expect_null(x$land)
  expect_null(x$land_lines)
  expect_null(x$edge_lines)
  expect_null(x$cells)
  expect_null(x$grid$values)
  expect_null(x$camera$distance)
})

test_that("hex_globe draws an H3 grid on the sphere only", {
  skip_if_not_installed("htmlwidgets")
  h3 <- hex_grid(resolution = 0, type = "h3")
  x <- hex_globe(h3, values = seq_along(h3_all_cells(0)), land = FALSE)$x
  expect_false(x$foldable)
  expect_true(x$cells$n_index > 0)
  expect_error(hex_globe(h3, surface = "icosahedron"), "ISEA grid")
})

test_that("hex_globe checks its arguments", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 1, aperture = 3)
  expect_error(hex_globe(g, values = 1:3), "one per cell")
  expect_error(hex_globe(g, tilt = 10), "perspective")
  expect_error(hex_globe(g, palette = "nope"), "palette")
})

test_that("hex_globe_png saves the view as a PNG of the asked size", {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if(is.null(suppressMessages(chromote::find_chrome())), "no Chromium browser")
  g <- hex_grid(resolution = 2, aperture = 3)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file))
  expect_equal(hex_globe_png(hex_globe(g), file, width = 200, height = 150,
                             scale = 2), file)
  head <- readBin(file, "raw", 24)
  expect_equal(head[2:4], charToRaw("PNG"))
  size <- function(b) sum(as.integer(b) * 256^(3:0))
  expect_equal(c(size(head[17:20]), size(head[21:24])), c(400, 300))
})

# =============================================================================
# The shader's cells against lonlat_to_cell()
# =============================================================================

skip_without_gpu_browser <- function() {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if(is.null(suppressMessages(chromote::find_chrome())), "no Chromium browser")
}

f32 <- function(x) readBin(writeBin(x, raw(), size = 4), "double", size = 4,
                           n = length(x))
unitize <- function(m) m / sqrt(rowSums(m^2))
xyz_lonlat <- function(p) {
  cbind(atan2(p[, 2], p[, 1]) * 180 / pi,
        atan2(p[, 3], sqrt(p[, 1]^2 + p[, 2]^2)) * 180 / pi)
}

# Points as the GPU reads them (32-bit floats): a third spread over the
# sphere, a third within a few cells of the face edges (which hold the quad
# seams) and a third within a few cells of the 12 vertices.
shader_probe <- function(n, grid) {
  s <- hexify:::icosa_solid(hexify:::orient_arg(grid))
  V <- unitize(s$vertices)
  E <- s$edges
  spread <- 3 * sqrt(4 * pi / grid_n_cells(grid))
  m <- n %/% 3
  z <- runif(m, -1, 1)
  phi <- runif(m, -pi, pi)
  k <- sample(nrow(E), m, TRUE)
  t <- runif(m)
  near <- rbind(V[E[k, 1], ] * (1 - t) + V[E[k, 2], ] * t,
                V[sample(12, n - 2 * m, TRUE), ])
  near <- unitize(near) + matrix(rnorm(3 * nrow(near)), ncol = 3) * spread
  p <- unitize(rbind(cbind(sqrt(1 - z^2) * cos(phi), sqrt(1 - z^2) * sin(phi), z),
                     near))
  matrix(f32(as.vector(p)), ncol = 3)
}

# Which of the points p (rows) have `cell` as lonlat_to_cell() of some point
# within `eps` radians: a disagreement that 32-bit rounding explains.
within_reach <- function(p, cell, grid, eps) {
  if (nrow(p) == 0) return(logical(0))
  a <- unitize(p)
  u <- cbind(-a[, 2], a[, 1], 0)
  polar <- abs(a[, 3]) > 0.9
  u[polar, ] <- cbind(0, -a[polar, 3], a[polar, 2])
  u <- unitize(u)
  w <- cbind(a[, 2] * u[, 3] - a[, 3] * u[, 2], a[, 3] * u[, 1] - a[, 1] * u[, 3],
             a[, 1] * u[, 2] - a[, 2] * u[, 1])
  th <- seq(0, 2 * pi, length.out = 65)[-65]
  hit <- rep(FALSE, nrow(a))
  for (t in th) {
    ll <- xyz_lonlat(unitize(a + eps * (cos(t) * u + sin(t) * w)))
    hit <- hit | lonlat_to_cell(ll[, 1], ll[, 2], grid) == cell
  }
  hit
}

# Grids from the coarsest to the finest the shader's 32-bit floats resolve.
# The shader places a point to a few millionths of a radian; a disagreement
# must be explained within 1e-5 radians (the largest seen is 5e-6, on face
# edges, about a pixel at the closest zoom). Past these grids the cells are
# within a few times that error, and near a vertex the shader can land on a
# centre the numbering gives to the neighbouring quad (4/7 at resolution 14,
# cells 1e-5 radians wide, does).
GLOBE_AGREEMENT_CASES <- c(
  lapply(c(0, 1, 2, 5, 10, 15, 18), function(r) list(ap = 3, res = r)),
  lapply(c(0, 1, 3, 8, 14), function(r) list(ap = 4, res = r)),
  lapply(c(0, 1, 2, 3, 6, 10), function(r) list(ap = 7, res = r)),
  lapply(c(2, 7, 12, 18), function(r) list(ap = "4/3", res = r)),
  lapply(c(3, 9, 12), function(r) list(ap = "4/7", res = r)),
  list(list(ap = c(4, 4, 7, 3), res = 4), list(ap = c(3, 7, 4, 7, 3, 4), res = 6)),
  list(list(ap = 3, res = 10, orient = c(-40, 20, 33)),
       list(ap = 7, res = 6, orient = c(0, 90, 0)),
       list(ap = "4/3", res = 12, orient = c(150, -60, 200)))
)
GLOBE_AGREEMENT_N <- 1e6
GLOBE_AGREEMENT_EPS <- 1e-5

test_that("the shader finds the cell lonlat_to_cell() finds", {
  skip_without_gpu_browser()
  set.seed(79)
  for (spec in GLOBE_AGREEMENT_CASES) {
    grid <- hex_grid(resolution = spec$res, aperture = spec$ap,
                     orientation = if (is.null(spec$orient)) "standard" else spec$orient)
    p <- shader_probe(GLOBE_AGREEMENT_N, grid)
    ll <- xyz_lonlat(p)
    ref <- lonlat_to_cell(ll[, 1], ll[, 2], grid)
    gpu <- hexify:::globe_shader_cells(grid, p)
    bad <- which(gpu != ref)
    unexplained <- bad[!within_reach(p[bad, , drop = FALSE], gpu[bad], grid,
                                     GLOBE_AGREEMENT_EPS)]
    expect_length(unexplained, 0)
    if (length(unexplained)) {
      print(cbind(lon = ll[unexplained, 1], lat = ll[unexplained, 2],
                  ref = ref[unexplained], gpu = gpu[unexplained]))
    }
  }
})

test_that("the cell under the centre of the view is the one looked at", {
  skip_without_gpu_browser()
  grid <- hex_grid(resolution = 5, aperture = 4)
  globe <- hex_globe(grid, values = seq_len(grid_n_cells(grid)) / 2,
                     center = c(lon = 15.3, lat = 32.1), land = FALSE)
  got <- hexify:::globe_in_chrome(globe, 201, 201, 1, 60, function(session, read) {
    read("document.querySelector('.hexify-globe').hexGlobe.pick(100, 100)")
  })
  id <- lonlat_to_cell(15.3, 32.1, grid)
  expect_equal(got$id, id)
  expect_equal(got$value, id / 2)
})

test_that("hex_globe_png takes only a globe", {
  expect_error(hex_globe_png(list(), tempfile()), "hex_globe")
})
