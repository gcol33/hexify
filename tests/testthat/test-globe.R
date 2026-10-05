# tests/testthat/test-globe.R
# Meshes and paths on the icosahedron and the sphere behind hex_globe(), and
# the widget itself.

mesh_points <- function(v) matrix(v, ncol = 3, byrow = TRUE)

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
  m <- hexify:::cpp_globe_faces(0.1)
  expect_equal(sum(mesh_item_area(m)), 4 * pi, tolerance = 1e-9)
  expect_equal(sqrt(rowSums(mesh_points(m$sphere)^2)),
               rep(1, length(m$item)), tolerance = 1e-12)
  expect_true(all(sqrt(rowSums(mesh_points(m$solid)^2)) <= 1 + 1e-12))
  expect_equal(sort(unique(m$item)), 1:20)
})

test_that("no edge of a refined mesh is longer than the spacing", {
  m <- hexify:::cpp_globe_faces(0.1)
  S <- mesh_points(m$solid)
  i <- matrix(m$index + 1L, ncol = 3, byrow = TRUE)
  # A face edge of the inscribed icosahedron is 1.0515 radii long.
  edge <- sqrt(rowSums((S[i[, 1], ] - S[i[, 2], ])^2))
  expect_lte(max(edge), 0.1 * 1.0515 + 1e-9)
})

test_that("cells tile the sphere and each fill matches its boundary", {
  for (spec in list(list(ap = 3, res = 3), list(ap = 4, res = 2),
                    list(ap = 7, res = 2), list(ap = 3, res = 0))) {
    g <- hex_grid(resolution = spec$res, aperture = spec$ap)
    ids <- seq_len(grid_n_cells(g))
    m <- hexify:::grid_surface_mesh(g, ids, 0.05)
    area <- mesh_item_area(m)
    expect_length(area, length(ids))
    expect_equal(sum(area), 4 * pi, tolerance = 1e-9,
                 info = paste("aperture", spec$ap, "res", spec$res))
  }
})

test_that("a mixed aperture sequence fills the sphere", {
  g <- hex_grid(resolution = 3, aperture = "4/3")
  ids <- seq_len(grid_n_cells(g))
  area <- mesh_item_area(hexify:::grid_surface_mesh(g, ids, 0.05))
  expect_equal(sum(area), 4 * pi, tolerance = 1e-9)
})

test_that("a cell crossing a face edge twice is filled once", {
  # Cell 25 of aperture 3 res 3 leaves face 0 for face 5 and for face 1 and
  # comes back each time; its fill matches the area its boundary encloses.
  g <- hex_grid(resolution = 3, aperture = 3)
  m <- hexify:::grid_surface_mesh(g, 25, 0.05)
  ref <- unname(cell_area(25, g)) / hexify:::grid_radius_km(g)^2
  expect_equal(unname(sum(mesh_item_area(m))), ref, tolerance = 0.01)
})

test_that("land triangles keep the land's area", {
  land <- hexify:::surface_land(TRUE)
  polys <- hexify:::sfc_polygons(land)
  m <- hexify:::cpp_globe_polygons(polys, 0.05)
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
  full <- mesh_item_area(hexify:::cpp_globe_polygons(list(list(outer)), 0.05))
  holed <- mesh_item_area(hexify:::cpp_globe_polygons(list(list(outer, hole)), 0.05))
  inner <- mesh_item_area(hexify:::cpp_globe_polygons(list(list(hole)), 0.05))
  expect_equal(unname(holed), unname(full - inner), tolerance = 1e-9)
})

test_that("a polygon wider than a hemisphere is refused", {
  ring <- cbind(c(-170, 0, 170, 0, -170), c(0, -80, 0, 80, 0))
  expect_error(hexify:::cpp_globe_polygons(list(list(ring)), 0.05),
               "88 degrees")
})

# =============================================================================
# Paths
# =============================================================================

test_that("paths on the faces keep both surfaces and stay on their faces", {
  m <- hexify:::cpp_sphere_paths_on_faces(c(0, 90, 179, -170), c(0, 10, 20, 60),
                                          rep(1L, 4), 0.02)
  expect_equal(sqrt(rowSums(m[, c("sphere_x", "sphere_y", "sphere_z")]^2)),
               rep(1, nrow(m)), tolerance = 1e-12)
  solid <- hexify:::icosa_solid()
  S <- m[, c("solid_x", "solid_y", "solid_z")]
  n <- solid$normals[m[, "face"] + 1L, ]
  h <- rowSums(solid$normals * solid$vertices[solid$faces[, 1], ])[m[, "face"] + 1L]
  expect_equal(rowSums(S * n), h, tolerance = 1e-9)
})

test_that("face edges lie on the edges of the solid", {
  m <- hexify:::edge_surface_paths(0.05)
  solid <- hexify:::icosa_solid()
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
  for (layer in c("surface", "land", "cells")) {
    expect_true(x[[layer]]$n_index > 0, info = layer)
  }
  for (layer in c("grid_lines", "land_lines", "edge_lines")) {
    expect_true(x[[layer]]$n_segment > 0, info = layer)
  }
  expect_equal(x$camera$distance, 2)
  expect_length(x$palette, 4 * 256)
})

test_that("hex_globe leaves out layers that are switched off", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 1, aperture = 4)
  x <- hex_globe(g, land = FALSE, face_edges = FALSE)$x
  expect_null(x$land)
  expect_null(x$land_lines)
  expect_null(x$edge_lines)
  expect_null(x$cells)
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

test_that("hex_globe_png takes only a globe", {
  expect_error(hex_globe_png(list(), tempfile()), "hex_globe")
})
