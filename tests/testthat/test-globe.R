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
  w <- hex_globe(g, values = ids, surface = "solid",
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
  expect_true(x$grid$all)
  expect_true(x$grid$valued)
  expect_false(x$grid$keyed)
  expect_equal(x$grid$ramp_map, c(1, 1 / (length(ids) - 1), 0))
  # A level per resolution from the grid's down to 0, each whole
  expect_equal(nrow(x$grid$levels), 3)
  expect_equal(x$grid$given, c(92, 32, 12))
  expect_length(b64_bytes(x$grid$uniform), 5184)
  expect_equal(vapply(x$grid$textures, `[[`, integer(1), "binding"), 6:8)
  expect_null(x$cells)
  expect_null(x$grid_lines)
  expect_equal(x$camera$distance, 2)
  expect_length(x$palette, 4 * 256)
})

test_that("an ISEA grid is sent as its frame and its cells' table, not outlines", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 10, aperture = 3)
  x <- hex_globe(g, land = FALSE)$x
  expect_true(x$grid$all)
  expect_null(x$grid$levels)
  expect_equal(vapply(x$grid$textures, function(t) prod(t$size), numeric(1)), c(1, 1, 1))
  expect_lt(object.size(x), 2e6)
  uniform <- readBin(b64_bytes(x$grid$uniform), "integer", n = 1296, size = 4,
                     endian = "little")
  expect_equal(uniform[17:19], c(0, 3^10, 0))

  # Four cells of 590,492 go into a perfect hash: their own slots, and the
  # padding of the quads their slots border.
  cells <- c(5e5, 17, 3, 590492)
  x <- hex_globe(g, values = c(1, NA, 3, 4), cells = cells, land = FALSE)$x
  expect_false(x$grid$all)
  expect_true(x$grid$keyed)
  expect_equal(x$grid$given[1], 4)
  words <- readBin(b64_bytes(x$grid$textures[[1]]$data), "integer", n = 1e5, size = 4,
                   endian = "little")
  stored <- words[words != -1L]
  float <- function(v) readBin(writeBin(v, raw(), size = 4), "integer", size = 4)
  expect_true(all(c(float(1), float(3), float(4), 0x7fc00000L) %in% stored))
  expect_error(hex_globe(g, cells = c(1, 590493)), "cell IDs")
  expect_error(hex_globe(g, cells = c(3, 3)), "repeat")
})

# A coarser level laid out slot by slot gathers each cell's children from the
# finer block through a stencil; packed by the hash, each given child pushes
# itself into the coarse cells nearest it. Both give every level the same
# cells, with the same words to the 16-bit step of their place and cover.
test_that("the two ways of building coarser levels agree cell by cell", {
  for (spec in list(list(ap = 3, res = 6), list(ap = 4, res = 5), list(ap = 7, res = 4),
                    list(ap = "4/7", res = 5), list(ap = c(4, 3, 7, 4), res = 4),
                    list(ap = 3, res = 5, poly = "octahedron"),
                    list(ap = 3, res = 6, orient = c(-40, 20, 33)),
                    list(ap = 3, res = 6, frac = 0.3), list(ap = 7, res = 4, frac = 0.5))) {
    g <- hex_grid(resolution = spec$res, aperture = spec$ap,
                  polyhedron = if (is.null(spec$poly)) "icosahedron" else spec$poly,
                  orientation = if (is.null(spec$orient)) "standard" else spec$orient)
    icosa <- hexify:::icosa_arg(g)
    n <- as.numeric(n_cells(g))
    set.seed(93)
    ids <- if (is.null(spec$frac)) seq_len(n) else sort(sample(n, round(spec$frac * n)))
    vals <- sin(ids / 37) + ifelse(ids %% 11 == 0, NA, 0)
    levels <- lapply(g@resolution:0, function(r) hexify:::isea_levels(g@aperture, r))
    build <- function(layout) {
      hexify:::cpp_globe_table(icosa, levels, bit64::as.integer64(ids), vals, FALSE,
                               c(-1, 0.5, 0), layout, TRUE)
    }
    slots <- build("slots")
    hashed <- build("hash")
    expect_false(slots$keyed)
    expect_true(hashed$keyed)
    expect_equal(slots$given, hashed$given)
    for (k in seq_along(slots$cells)[-1]) {
      a <- slots$cells[[k]]
      b <- hashed$cells[[k]]
      expect_equal(sort(a$idx), sort(b$idx))
      wa <- a$word[order(a$idx)]
      wb <- b$word[order(b$idx)]
      expect_lte(max(abs(wa %/% 65536 - wb %/% 65536), 0), 1)
      expect_lte(max(abs(wa %% 65536 - wb %% 65536), 0), 1)
    }
  }
})

# Every slot of a table laid out slot by slot holds the value of the cell
# owning its point, as cpp_quad_ij_to_cell() moves a point into its quad (for
# apertures that store substrate coordinates, all but 7), and absent where
# no quad owns it; every coarser level covers each of its cells whole.
test_that("the cells' table holds every cell where the shader reads it", {
  for (spec in list(list(ap = 3, res = 5), list(ap = 4, res = 4), list(ap = 3, res = 0),
                    list(ap = "4/3", res = 5), list(ap = c(4, 3, 4), res = 3),
                    list(ap = 4, res = 3, poly = "octahedron"),
                    list(ap = 3, res = 4, orient = c(-40, 20, 33)))) {
    g <- hex_grid(resolution = spec$res, aperture = spec$ap,
                  polyhedron = if (is.null(spec$poly)) "icosahedron" else spec$poly,
                  orientation = if (is.null(spec$orient)) "standard" else spec$orient)
    icosa <- hexify:::icosa_arg(g)
    lv <- hexify:::isea_levels(g@aperture, g@resolution)
    n <- as.numeric(n_cells(g))
    value <- function(id) id %% 1000003 + 0.5
    levels <- lapply(g@resolution:0, function(r) hexify:::isea_levels(g@aperture, r))
    tb <- hexify:::cpp_globe_table(icosa, levels, bit64::integer64(0), value(seq_len(n)),
                                   TRUE, c(0, 1, 0))
    expect_false(tb$keyed)
    expect_equal(tb$given, vapply(levels, function(l) {
      as.numeric(hexify:::cpp_globe_frame(icosa, l$resolution, l$aperture, l$ap_seq)$n_cells)
    }, numeric(1)))
    words <- readBin(b64_bytes(tb$values), "integer", n = 1e7, size = 4, endian = "little")
    L <- tb$levels[1, ]
    quads <- nrow(hexify:::icosa_solid(icosa)$vertices) - 2
    slot <- seq_len(L[["hp"]] * L[["wp"]] * quads) - 1
    q <- slot %/% (L[["hp"]] * L[["wp"]]) + 1
    row <- (slot %% (L[["hp"]] * L[["wp"]])) %/% L[["wp"]]
    col <- slot %% L[["wp"]]
    u <- row - L[["pr"]]
    v <- L[["index"]] * (col - L[["pc"]]) + (L[["c"]] * u) %% L[["index"]]
    owner <- as.numeric(hexify:::cpp_quad_ij_to_cell(icosa, as.integer(q), u, v,
                                                     lv$resolution, lv$aperture, lv$ap_seq))
    expected <- ifelse(is.na(owner), -1L,
                       readBin(writeBin(value(owner), raw(), size = 4), "integer",
                               n = length(owner), size = 4))
    expect_equal(words[slot + 1], expected)
    # Coarser levels: the low 16 bits are the cover, whole everywhere
    total <- sum(tb$levels[, "hp"] * tb$levels[, "wp"] * quads)
    coarse <- if (total > length(slot)) words[(length(slot) + 1):total]
    coarse <- coarse[coarse != -1L]
    if (length(coarse)) expect_true(all(bitwAnd(coarse, 65535L) == 65534L))
  }
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
  f <- hexify:::cpp_globe_frame(numeric(0), 3L, 3L, integer(0))
  expect_equal(f[c("dim", "index", "c")], list(dim = 9, index = 3, c = 2))
  expect_equal(f$generator, c(2, 1))
  f <- hexify:::cpp_globe_frame(numeric(0), 3L, 7L, integer(0))
  expect_equal(f[c("dim", "index", "c")], list(dim = 49, index = 7, c = 5))
  f <- hexify:::cpp_globe_frame(numeric(0), 4L, 4L, integer(0))
  expect_equal(f[c("dim", "index", "per_quad")], list(dim = 16, index = 1, per_quad = 256))
  f <- hexify:::cpp_globe_frame(numeric(0), 2L, 0L, c(4L, 3L, 7L))
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
  expect_error(hex_globe(h3, surface = "solid"), "ISEA grid")
})

test_that("hex_globe checks its arguments", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 1, aperture = 3)
  expect_error(hex_globe(g, values = 1:3), "one per cell")
  expect_error(hex_globe(g, tilt = 10), "perspective")
  expect_error(hex_globe(g, palette = "nope"), "palette")
  expect_error(hex_globe(g, smooth = NA), "smooth")
  expect_error(hex_globe(hex_grid(resolution = 0, type = "h3"), smooth = TRUE), "ISEA")
})

test_that("hex_globe_png in Chrome saves the view as a PNG of the asked size", {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if(is.null(suppressMessages(chromote::find_chrome())), "no Chromium browser")
  g <- hex_grid(resolution = 2, aperture = 3)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file))
  expect_equal(hex_globe_png(hex_globe(g), file, width = 200, height = 150,
                             scale = 2, renderer = "chrome"), file)
  head <- readBin(file, "raw", 24)
  expect_equal(head[2:4], charToRaw("PNG"))
  size <- function(b) sum(as.integer(b) * 256^(3:0))
  expect_equal(c(size(head[17:20]), size(head[21:24])), c(400, 300))

  # The page reads the PNG back: the globe fills the centre, the corner is
  # off it.
  png64 <- paste(sprintf("%02x", as.integer(readBin(file, "raw", file.size(file)))),
                 collapse = "")
  alpha <- hexify:::globe_in_chrome(hex_globe(g), 50, 50, 1, 60, function(session, read) {
    read(paste0(
      "(async () => {",
      "  const hex = '", png64, "';",
      "  const bytes = new Uint8Array(hex.length / 2);",
      "  for (let i = 0; i < bytes.length; i++) bytes[i] = parseInt(hex.substr(2 * i, 2), 16);",
      "  const img = await createImageBitmap(new Blob([bytes], { type: 'image/png' }));",
      "  const c = document.createElement('canvas'); c.width = img.width; c.height = img.height;",
      "  const x = c.getContext('2d'); x.drawImage(img, 0, 0);",
      "  return [x.getImageData(200, 150, 1, 1).data[3], x.getImageData(2, 2, 1, 1).data[3]];",
      "})()"))
  })
  expect_equal(unlist(alpha), c(255, 0))
})

# =============================================================================
# The globe drawn through wgpu (hexglobe)
# =============================================================================

f32_words <- function(bytes) readBin(bytes, "double", n = length(bytes) / 4,
                                     size = 4, endian = "little")

test_that("a scene holds the widget's layers in its drawing order", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 2, aperture = 3)
  x <- hex_globe(g, values = seq_len(n_cells(g)), na_fill = "grey")$x
  s <- hexify:::globe_scene(x, 120, 80, 2)
  expect_equal(c(s$width, s$height), c(240L, 160L))
  expect_identical(s$shader, x$shader)
  expect_equal(vapply(s$layers, `[[`, integer(1), "kind"), c(0L, 0L, 2L, 3L, 3L))
  # The ocean and the grid are drawn on one faces mesh, uploaded once.
  expect_identical(s$layers[[1]]$pos, s$layers[[3]]$pos)
  expect_length(s$layers[[3]]$tri, length(s$layers[[3]]$pos) / 3)
  expect_equal(s$layers[[1]]$count, x$surface$n_index)
  # Layer uniforms: colour, lift, line width in device pixels, shaded,
  # coloured by value.
  u <- f32_words(s$layers[[3]]$uniform)
  expect_length(u, 12)
  expect_equal(u[1:8], c(x$style$grid_border, hexify:::GLOBE_LIFT$cells,
                         2 * x$style$grid_lwd, 1, 1), tolerance = 1e-6)
  coast <- f32_words(s$layers[[4]]$uniform)
  expect_equal(coast[5:6], c(hexify:::GLOBE_LIFT$coast, 2 * x$style$land_lwd),
               tolerance = 1e-6)
  expect_length(s$layers[[3]]$grid, 5184)
  grid <- readBin(s$layers[[3]]$grid, "integer", n = 1296, size = 4, endian = "little")
  # every cell drawn, values given, no smooth fill, 20 faces; a level per
  # resolution, laid out slot by slot
  expect_equal(grid[13:16], c(1, 1, 0, 20))
  expect_equal(grid[17:20], c(0, 9, 3, 0))
  expect_equal(f32_words(s$layers[[3]]$grid)[25:28], x$style$na_fill, tolerance = 1e-6)
  # Level 0 is the grid's own frame: a quad side of 3, the aligned lattice
  expect_equal(grid[785:787], c(3, 1, 0))
  expect_equal(vapply(s$layers[[3]]$textures, `[[`, integer(1), "binding"), 6:8)
  expect_type(s$layers[[3]]$textures[[1]]$data, "raw")
})

test_that("the scene's camera is the plot() method's view", {
  skip_if_not_installed("htmlwidgets")
  g <- hex_grid(resolution = 1, aperture = 4)
  x <- hex_globe(g, center = c(lon = 40, lat = -10), rotation = 20, land = FALSE)$x
  cam <- f32_words(hexify:::globe_camera_uniform(x, 300, 200))
  expect_length(cam, 28)
  view <- hexify:::surface_view(c(lon = 40, lat = -10), rotation = 20)
  expect_equal(cam[c(1:3, 5:7, 9:11)], as.vector(view$cam), tolerance = 1e-6)
  expect_equal(cam[13:24], c(0, 0, 0, 0, 0, 0, 1.02, 1, 300, 200, 0, 1), tolerance = 1e-6)
  expect_equal(cam[25:28], c(view$light, 1), tolerance = 1e-6)

  x <- hex_globe(g, projection = "perspective", distance = 2, tilt = 25,
                 surface = "solid", land = FALSE)$x
  cam <- f32_words(hexify:::globe_camera_uniform(x, 300, 300))
  view <- hexify:::surface_view(c(lon = 15, lat = 32), distance = 2, tilt = 25)
  frame <- hexify:::view_frame(view, NA)
  reach <- sqrt(sum(view$eye^2))
  expect_equal(cam[13:24], c(view$eye, 1, frame, view$scale, 300, 300, view$near,
                             reach + 1.5), tolerance = 1e-6)
  expect_equal(cam[28], 0)
})

test_that("unsigned words above 2^31 keep their bits", {
  b <- hexify:::le32(c(0, 2^31, 2^32 - 1, 7), "u32")
  expect_equal(as.integer(b), c(0, 0, 0, 0, 0, 0, 0, 128, 255, 255, 255, 255, 7, 0, 0, 0))
  expect_equal(hexify:::le32(c(-1, -2^31, 5), "i32"),
               writeBin(c(-1L, NA_integer_, 5L), raw(), size = 4L, endian = "little"))
})

skip_without_gpu_browser <- function() {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if(is.null(suppressMessages(chromote::find_chrome())), "no Chromium browser")
}

skip_without_wgpu <- function() {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("hexglobe")
  skip_if(inherits(try(hexglobe::gpu_adapter(), silent = TRUE), "try-error"),
          "no graphics adapter for wgpu")
}

test_that("hex_globe_png through wgpu saves the view as a PNG of the asked size", {
  skip_without_wgpu()
  g <- hex_grid(resolution = 2, aperture = 3)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file))
  expect_equal(hex_globe_png(hex_globe(g), file, width = 200, height = 150,
                             scale = 2, renderer = "wgpu"), file)
  img <- hexglobe::read_png(file)
  expect_equal(c(img$width, img$height), c(400L, 300L))
  alpha <- matrix(as.integer(img$rgba[seq(4, length(img$rgba), by = 4)]), 400)
  expect_equal(alpha[c(201, 3), c(151, 3)][c(1, 4)], c(255, 0))
})

# The two renderers draw one widget with one shader; what is left between
# them is the graphics cards' rasterisation and anti-aliasing of edges.
test_that("wgpu and Chrome draw the same pixels", {
  skip_without_wgpu()
  skip_without_gpu_browser()
  g <- hex_grid(resolution = 3, aperture = 3)
  cells <- seq_len(n_cells(g))
  globes <- list(
    hex_globe(g),
    hex_globe(g, surface = "solid", center = "pacific"),
    hex_globe(g, values = cell_to_lonlat(cells, g)$lat_deg, cells = cells,
              palette = "Blue-Red 3"),
    hex_globe(g, projection = "perspective", distance = 2.2, tilt = 30,
              rotation = 15),
    hex_globe(hex_grid(resolution = 1, type = "h3"),
              values = seq_along(h3_all_cells(1)), land = FALSE),
    hex_globe(hex_grid(resolution = 2, aperture = 4, polyhedron = "octahedron"),
              surface = "solid")
  )
  files <- replicate(2, tempfile(fileext = ".png"))
  on.exit(unlink(files))
  for (w in globes) {
    hex_globe_png(w, files[1], 240, 200, renderer = "wgpu")
    hex_globe_png(w, files[2], 240, 200, renderer = "chrome")
    a <- hexglobe::read_png(files[1])
    b <- hexglobe::read_png(files[2])
    expect_equal(c(a$width, a$height), c(b$width, b$height))
    d <- apply(matrix(abs(as.integer(a$rgba) - as.integer(b$rgba)), 4), 2, max)
    expect_lt(mean(d > 2), 0.01)
    expect_lt(mean(d > 16), 0.002)
  }
})

test_that("base64 text decodes to its bytes", {
  for (x in list(0, c(1, 2, 3), c(4294967295, 7, 65536, 1))) {
    text <- hexify:::cpp_base64_buffer(x, "u32")
    bytes <- hexify:::cpp_base64_decode(text)
    expect_equal(readBin(bytes, "integer", n = length(x), size = 4, endian = "little"),
                 ifelse(x > 2^31 - 1, x - 2^32, x))
  }
  expect_equal(hexify:::cpp_base64_decode("TWE="), charToRaw("Ma"))
  expect_equal(hexify:::cpp_base64_decode("TQ=="), charToRaw("M"))
  expect_equal(hexify:::cpp_base64_decode(""), raw(0))
})

# =============================================================================
# The shader's cells against lonlat_to_cell()
# =============================================================================

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
  s <- hexify:::icosa_solid(hexify:::icosa_arg(grid))
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
       list(ap = "4/3", res = 12, orient = c(150, -60, 200))),
  list(list(ap = 3, res = 10, proj = "fuller"),
       list(ap = 4, res = 8, proj = "fuller"),
       list(ap = 7, res = 5, proj = "fuller"),
       list(ap = "4/3", res = 12, proj = "fuller", orient = c(-40, 20, 33))),
  list(list(ap = 3, res = 10, proj = "ivea"),
       list(ap = 4, res = 8, proj = "ivea"),
       list(ap = 7, res = 5, proj = "ivea"),
       list(ap = "4/3", res = 12, proj = "ivea", orient = c(-40, 20, 33)))
)
GLOBE_AGREEMENT_N <- 1e6
GLOBE_AGREEMENT_EPS <- 1e-5

# A value every cell can hold exactly as a 32-bit float, and that float's bits
cell_value <- function(id) id %% 1000003 + 0.5
float_word <- function(x) {
  w <- readBin(writeBin(x, raw(), size = 4), "integer", n = length(x), size = 4)
  ifelse(w < 0, w + 2^32, w)
}

# The table the shader reads holds each cell's value where the shader finds
# the cell. Grids up to a few million cells carry a value on every cell, laid
# out slot by slot; finer ones carry values on a sample of the cells the
# probes hit, packed by the perfect hash, and a cell outside the sample reads
# as absent. Either way the word read for a point is the one of the cell the
# shader numbers there, though the read takes no step into the quad that owns
# the cell.
test_that("the shader finds the cell lonlat_to_cell() finds, and reads its value", {
  skip_without_gpu_browser()
  set.seed(79)
  for (spec in GLOBE_AGREEMENT_CASES) {
    grid <- hex_grid(resolution = spec$res, aperture = spec$ap,
                     orientation = if (is.null(spec$orient)) "standard" else spec$orient,
                     projection = if (is.null(spec$proj)) "isea" else spec$proj)
    p <- shader_probe(GLOBE_AGREEMENT_N, grid)
    ll <- xyz_lonlat(p)
    ref <- lonlat_to_cell(ll[, 1], ll[, 2], grid)
    n <- as.numeric(n_cells(grid))
    if (n <= 4e6) {
      cells <- NULL
      gpu <- hexify:::globe_shader_cells(grid, p, values = cell_value(seq_len(n)))
    } else {
      cells <- unique(ref[seq(1, length(ref), by = 2)])
      gpu <- hexify:::globe_shader_cells(grid, p, cells = cells,
                                         values = cell_value(as.numeric(cells)))
    }
    bad <- which(gpu$id != ref)
    unexplained <- bad[!within_reach(p[bad, , drop = FALSE], gpu$id[bad], grid,
                                     GLOBE_AGREEMENT_EPS)]
    expect_length(unexplained, 0)
    if (length(unexplained)) {
      print(cbind(lon = ll[unexplained, 1], lat = ll[unexplained, 2],
                  ref = ref[unexplained], gpu = gpu$id[unexplained]))
    }
    given <- is.null(cells) | gpu$id %in% as.numeric(cells)
    expected <- ifelse(given, float_word(cell_value(gpu$id)), 2^32 - 1)
    wrong <- which(gpu$word != expected)
    expect_length(wrong, 0)
    if (length(wrong)) {
      print(cbind(lon = ll[wrong, 1], lat = ll[wrong, 2], id = gpu$id[wrong],
                  word = gpu$word[wrong], expected = expected[wrong]))
    }
    expect_equal(gpu$flags, ifelse(given, 7, 1))
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
  expect_equal(got$id, as.numeric(id))
  expect_equal(got$value, as.numeric(id) / 2)
})

test_that("hex_globe_png takes only a globe", {
  expect_error(hex_globe_png(list(), tempfile()), "hex_globe")
})
