# Hex9 (Griffin 2026): the shifted-aperture-9 grid of the octahedron

hex9_ids <- function(n) bit64::as.integer64(seq_len(n))

test_that("aperture 9 builds Hex9 on the octahedron", {
  g <- hex_grid(resolution = 3, aperture = 9)
  expect_equal(g@polyhedron, "octahedron")
  expect_equal(g@aperture, "9")
  expect_equal(n_cells(g), 12 * 9^3)
  expect_error(hex_grid(resolution = 3, aperture = 9, polyhedron = "icosahedron"),
               "octahedron only")
  expect_error(hex_grid(resolution = 19, aperture = 9), "64-bit")
  expect_silent(hex_grid(resolution = 18, aperture = 9))
  expect_error(hex_grid(resolution = 3, aperture = c(9, 3)), "Aperture sequence")
})

test_that("area_km2 picks the Hex9 resolution", {
  surface <- body_surface_km2(EARTH_RADIUS_KM)
  for (L in 0:6) {
    g <- hex_grid(area_km2 = surface / (12 * 9^L), aperture = 9)
    expect_equal(g@resolution, L)
    expect_equal(g@area_km2, surface / (12 * 9^L))
  }
})

test_that("every cell centre falls in its own cell and IDs run 1..n", {
  for (L in 0:3) {
    g <- hex_grid(resolution = L, aperture = 9)
    ids <- hex9_ids(n_cells(g))
    cc <- cell_to_lonlat(ids, g)
    expect_identical(lonlat_to_cell(cc$lon_deg, cc$lat_deg, g), ids)
  }
})

test_that("labels name every cell once and read back", {
  for (L in 0:3) {
    g <- hex_grid(resolution = L, aperture = 9)
    ids <- hex9_ids(n_cells(g))
    lab <- cell_to_index(ids, g)
    expect_false(anyDuplicated(lab) > 0)
    expect_true(all(nchar(lab) == L + 3L))
    back <- cpp_hex9_parse_label(lab)
    expect_identical(back$cell_id, ids)
    expect_true(all(back$resolution == L))
    body <- sub("\\..*$", "", lab)
    expect_identical(cpp_hex9_parse_label(body)$cell_id, ids)
  }
  expect_true(is.na(cpp_hex9_parse_label("43.5")$cell_id))
  expect_true(is.na(cpp_hex9_parse_label("c3")$cell_id))
})

test_that("every cell has nine children whose parent it is", {
  for (L in 0:2) {
    g <- hex_grid(resolution = L, aperture = 9)
    gc <- hex_grid(resolution = L + 1, aperture = 9)
    ids <- hex9_ids(n_cells(g))
    kids <- get_children(ids, g)
    expect_true(all(lengths(kids) == 9L))
    flat <- cell_id_unlist(kids)
    expect_identical(sort(flat), hex9_ids(n_cells(gc)))
    expect_identical(get_parent(flat, gc), rep(ids, each = 9))
  }
  g <- hex_grid(resolution = 1, aperture = 9)
  two <- get_children(bit64::as.integer64(5), g, levels = 2)[[1]]
  g3 <- hex_grid(resolution = 3, aperture = 9)
  expect_length(two, 81)
  expect_true(all(get_parent(two, g3, levels = 2) == 5))
})

test_that("twelve cells, two at each vertex, have five neighbours", {
  for (L in 0:2) {
    g <- hex_grid(resolution = L, aperture = 9)
    ids <- hex9_ids(n_cells(g))
    nb <- get_neighbors(ids, g)
    expect_equal(sum(lengths(nb) == 5L), 12L)
    expect_equal(sum(lengths(nb) == 6L), n_cells(g) - 12)
    pairs <- data.frame(a = rep(as.character(ids), lengths(nb)),
                        b = as.character(cell_id_unlist(nb)))
    back <- paste(pairs$b, pairs$a) %in% paste(pairs$a, pairs$b)
    expect_true(all(back))
  }
})

test_that("the grid's digits agree with libhex9 on the octahedron", {
  pts <- utils::read.csv(test_path("data", "libhex9_points.csv"),
                         colClasses = "character")
  x <- as.numeric(pts$x); y <- as.numeric(pts$y); z <- as.numeric(pts$z)
  for (col in grep("^L", names(pts), value = TRUE)) {
    L <- as.integer(sub("L", "", col))
    ids <- cpp_hex9_octahedron_cell(x, y, z, L)
    expect_identical(cpp_hex9_label(ids, L), pts[[col]], info = col)
  }
})

test_that("neighbours, parents and children agree with libhex9", {
  ref <- utils::read.csv(test_path("data", "libhex9_cells_L1.csv"),
                         colClasses = "character")
  g1 <- hex_grid(resolution = 1, aperture = 9)
  ids <- cpp_hex9_parse_label(ref$label)$cell_id
  expect_equal(length(ids), 108L)
  nb <- lapply(get_neighbors(ids, g1), function(x) sort(cell_to_index(x, g1)))
  expect_identical(nb, lapply(strsplit(ref$neighbours, ";"), sort))
  g0 <- hex_grid(resolution = 0, aperture = 9)
  expect_identical(cell_to_index(get_parent(ids, g1), g0), ref$parent)
  g2 <- hex_grid(resolution = 2, aperture = 9)
  kids <- lapply(get_children(ids, g1), function(x) sort(cell_to_index(x, g2)))
  expect_identical(kids, lapply(strsplit(ref$children, ";"), sort))
})

test_that("ancestors and owned cells two levels apart are libhex9's", {
  # hex9_owned_cells() of 35 level-1 cells at level 3
  # (paper/bench/hex9_libhex9_dump.cpp, mode owned 1 2 40 7)
  ref <- utils::read.csv(test_path("data", "libhex9_owned.csv"),
                         colClasses = "character")
  g1 <- hex_grid(resolution = 1, aperture = 9)
  g3 <- hex_grid(resolution = 3, aperture = 9)
  ids <- cpp_hex9_parse_label(ref$label)$cell_id
  owned <- lapply(get_children(ids, g1, levels = 2), function(x) sort(cell_to_index(x, g3)))
  expect_identical(owned, lapply(strsplit(ref$owned, ";"), sort))
  flat <- cell_id_unlist(get_children(ids, g1, levels = 2))
  expect_identical(get_parent(flat, g3, levels = 2), rep(ids, each = 81))
})

test_that("owned cells partition a level, and differ from children's children", {
  g1 <- hex_grid(resolution = 1, aperture = 9)
  g3 <- hex_grid(resolution = 3, aperture = 9)
  ids <- hex9_ids(n_cells(g1))
  owned <- cell_id_unlist(get_children(ids, g1, levels = 2))
  expect_identical(sort(owned), hex9_ids(n_cells(g3)))
  g2 <- hex_grid(resolution = 2, aperture = 9)
  fine <- hex9_ids(n_cells(g3))
  lineage <- get_parent(get_parent(fine, g3), g2)
  # libhex9 counts the two relations apart on one cell in nine
  expect_equal(mean(lineage != get_parent(fine, g3, levels = 2)), 1 / 9)
})

test_that("Hex9 cells are equal-area on the equal-area projections", {
  for (projection in c("isea", "ivea")) {
    g <- hex_grid(resolution = 2, aperture = 9, projection = projection)
    ids <- hex9_ids(n_cells(g))
    a <- cell_area(ids, g)
    expect_equal(unname(a), rep(g@area_km2, length(ids)))
    expect_equal(sum(a), body_surface_km2(EARTH_RADIUS_KM))
    expect_false(any(is_pentagon(ids, g)))
    sr <- cpp_cell_solid_angle(icosa_arg(g), ids, 2L, 9L, integer(0),
                               CELL_WALL_TOLERANCE)
    expect_equal(sr, rep(4 * pi / length(ids), length(ids)), tolerance = 1e-6)
  }
})

test_that("Hex9 polygons tile the sphere", {
  skip_if_not_installed("sf")
  g <- hex_grid(resolution = 1, aperture = 9)
  cells <- grid_global(g, wrap_dateline = FALSE)
  expect_equal(nrow(cells), 108L)
  expect_true(all(sf::st_is_valid(cells)))
  old <- sf::sf_use_s2(TRUE)
  on.exit(sf::sf_use_s2(old), add = TRUE)
  area <- as.numeric(sf::st_area(cells)) / 1e6
  expect_equal(sum(area), body_surface_km2(EARTH_RADIUS_KM), tolerance = 0.01)
})

test_that("cell IDs follow the solid under any orientation", {
  g <- hex_grid(resolution = 2, aperture = 9)
  r <- hex_grid(resolution = 2, aperture = 9, orientation = c(20, 70, 35))
  ids <- hex9_ids(n_cells(g))
  expect_identical(get_neighbors(ids, r), get_neighbors(ids, g))
  expect_identical(get_parent(ids, r), get_parent(ids, g))
  cc <- cell_to_lonlat(ids, r)
  expect_identical(lonlat_to_cell(cc$lon_deg, cc$lat_deg, r), ids)
})

test_that("digit strings name one cell at every level", {
  expect_true(cpp_hex9_digit_strings_unique())
})

test_that("Hex9 cells compact into their parents and back", {
  g <- hex_grid(resolution = 2, aperture = 9)
  g3 <- hex_grid(resolution = 3, aperture = 9)
  kids <- get_children(bit64::as.integer64(c(17, 500)), g)
  lab <- cell_to_index(cell_id_unlist(kids), g3)
  stray <- cell_to_index(bit64::as.integer64(8000), g3)
  out <- hex_compact(c(lab, stray), g)
  expect_setequal(out, c(cell_to_index(bit64::as.integer64(c(17, 500)), g), stray))
  expect_setequal(hex_uncompact(out, g, 3L), c(lab, stray))
})

test_that("a vertex cell meets its partner across two walls", {
  g <- hex_grid(resolution = 1, aperture = 9)
  ids <- hex9_ids(n_cells(g))
  w <- wall_metrics(ids, g)
  nb <- get_neighbors(ids, g)
  five <- ids[lengths(nb) == 5L]
  expect_length(five, 12L)
  # Each wall once: six per cell, two cells per wall
  expect_equal(nrow(w), 6L * length(ids) / 2L)
  pair <- paste(pmin(w$cell_id, w$neighbor_id), pmax(w$cell_id, w$neighbor_id))
  twice <- names(which(table(pair) == 2L))
  expect_length(twice, 6L)
  m <- cell_metrics(ids, g)
  expect_equal(nrow(m), length(ids))
})

test_that("a Hex9 grid has its own DGGRS definition and no DGGRID one", {
  def <- dggrs_definition(hex_grid(resolution = 3, aperture = 9))
  expect_equal(def$title, "Hex9")
  expect_equal(def$dggh$definition$zoneTypes, "hexagon")
  expect_equal(def$dggh$definition$refinementRatio, 9L)
  expect_error(as_dggrid(hex_grid(resolution = 3, aperture = 9)), "no aperture 9")
})
