# hex_aggregate(): values of cells divided among the coarser cells they
# overlap, by area share or wholly to the parent holding the centre.

all_cells <- function(g) hexify:::grid_cells(g)
coarser <- function(g, levels = 1L) hexify:::grid_at_resolution(g, g@resolution - levels)

test_that("the global total is the same at every level under both rules", {
  skip_on_cran()
  grids <- list(hex_grid(resolution = 4, aperture = 3),
                hex_grid(resolution = 3, aperture = 4),
                hex_grid(resolution = 3, aperture = 7),
                hex_grid(resolution = 4, aperture = c(4, 7, 3, 4)),
                hex_grid(resolution = 3, type = "h3"))
  set.seed(91)
  for (g in grids) {
    ids <- all_cells(g)
    v <- rexp(length(ids))
    for (rule in c("area", "centre")) {
      for (lv in seq_len(min(3L, g@resolution))) {
        up <- hex_aggregate(ids, v, g, levels = lv, rule = rule)
        expect_equal(sum(up$value), sum(v), tolerance = 1e-12,
                     info = paste(g@aperture, rule, lv))
        expect_equal(nrow(up), as.numeric(grid_n_cells(coarser(g, lv))),
                     info = paste(g@aperture, rule, lv))
      }
    }
  }
})

test_that("a uniform density gives every area-weighted parent its own area", {
  # Hexagons all get the same total, vertex cells their sides / 6 of it; the
  # cases cover each aperture, a family and a per-level sequence, the
  # octahedron, and the vertex-oriented projection
  grids <- list(hex_grid(resolution = 4, aperture = 3),
                hex_grid(resolution = 3, aperture = 4),
                hex_grid(resolution = 3, aperture = 7),
                hex_grid(resolution = 5, aperture = "4/3"),
                hex_grid(resolution = 4, aperture = c(3, 7, 4, 3)),
                hex_grid(resolution = 4, aperture = 3, polyhedron = "octahedron"),
                hex_grid(resolution = 4, aperture = 3, polyhedron = "tetrahedron"),
                hex_grid(resolution = 3, aperture = 4, polyhedron = "tetrahedron"),
                hex_grid(resolution = 3, aperture = 7, polyhedron = "tetrahedron"),
                hex_grid(resolution = 3, aperture = 7, projection = "ivea"))
  for (g in grids) {
    ids <- all_cells(g)
    density <- unname(cell_area(ids, g))
    for (lv in 1:2) {
      up <- hex_aggregate(ids, density, g, levels = lv)
      area <- unname(cell_area(up$cell_id, coarser(g, lv)))
      expect_equal(up$value, area, tolerance = 1e-12,
                   info = paste(g@aperture, grid_polyhedron(g), grid_projection(g), lv))
      mean_up <- hex_aggregate(ids, rep(2.5, length(ids)), g, levels = lv,
                               measure = "mean")
      expect_equal(mean_up$value, rep(2.5, nrow(mean_up)), tolerance = 1e-12)
    }
    pent <- is_pentagon(up$cell_id, coarser(g, 2L))
    if (any(pent)) {
      expect_equal(up$value[pent] / max(up$value), rep(5 / 6, sum(pent)),
                   tolerance = 1e-12)
    }
    tri <- hexify:::isea_cell_sides(up$cell_id, coarser(g, 2L)) == 3
    if (any(tri)) {
      expect_equal(up$value[tri] / max(up$value), rep(1 / 2, sum(tri)), tolerance = 1e-12)
    }
  }
})

test_that("the centre rule gives each cell wholly to get_parent()", {
  g <- hex_grid(resolution = 4, aperture = 3)
  ids <- all_cells(g)
  set.seed(3)
  v <- rpois(length(ids), 4)
  up <- hex_aggregate(ids, v, g, levels = 2, rule = "centre")
  ref <- tapply(v, as.character(get_parent(ids, g, levels = 2)), sum)
  expect_equal(up$value, unname(as.numeric(ref[as.character(up$cell_id)])))
})

test_that("the lattice shares are the shares measured on the sphere", {
  # Each cell clipped against the coarser cells it overlaps face by face and
  # the pieces measured on the sphere: the shares 1/3, 1/2, 11/12 and 1/12
  set.seed(5)
  for (case in list(list(3, 5), list(4, 4), list(7, 3),
                    list(c(4, 7, 3, 7), 4), list(7, 2, "octahedron"),
                    list(3, 4, "tetrahedron"), list(4, 3, "tetrahedron"),
                    list(7, 3, "tetrahedron"))) {
    g <- hex_grid(resolution = case[[2]], aperture = case[[1]],
                  polyhedron = if (length(case) > 2) case[[3]] else "icosahedron")
    ids <- all_cells(g)
    ids <- ids[sort(sample(length(ids), min(60L, length(ids))))]
    ov <- get_parent(ids, g, overlapping = TRUE)
    lattice <- hexify:::lattice_shares(lengths(ov), hexify:::step_aperture(g))
    measured <- hexify:::sphere_shares(ids, ov, g)
    expect_lt(max(abs(measured - lattice)), 1e-10)
  }
})

test_that("every cell around a vertex divides as on the plane", {
  # The cells around the twelve vertices, where the angle deficit is
  for (case in list(c(3, 6), c(4, 5), c(7, 4))) {
    g <- hex_grid(resolution = case[2], aperture = case[1])
    pent <- all_cells(g)[is_pentagon(all_cells(g), g)]
    near <- unique(hexify:::cell_id_unlist(get_neighbors(pent, g, k = 1,
                                                          include_self = TRUE)))
    ov <- get_parent(near, g, overlapping = TRUE)
    lattice <- hexify:::lattice_shares(lengths(ov), case[1])
    measured <- hexify:::sphere_shares(near, ov, g)
    expect_lt(max(abs(measured - lattice)), 1e-11)
  }
})

test_that("every cell around a vertex of the tetrahedron divides as on the plane", {
  # The solid turns half a turn about each vertex, a symmetry of the finer and
  # of the coarser lattice alike, so the shares there are the plane's
  for (case in list(c(3, 5), c(4, 4), c(7, 3))) {
    g <- hex_grid(resolution = case[2], aperture = case[1], polyhedron = "tetrahedron")
    ids <- all_cells(g)
    vertex <- ids[hexify:::isea_cell_sides(ids, g) == 3]
    expect_length(vertex, 4)
    near <- unique(hexify:::cell_id_unlist(get_neighbors(vertex, g, k = 1,
                                                          include_self = TRUE)))
    ov <- get_parent(near, g, overlapping = TRUE)
    lattice <- hexify:::lattice_shares(lengths(ov), case[1])
    measured <- hexify:::sphere_shares(near, ov, g)
    expect_lt(max(abs(measured - lattice)), 1e-11)
  }
})

test_that("each step of a mixed sequence takes the shares of its aperture", {
  g <- hex_grid(resolution = 4, aperture = c(4, 7, 3, 4))
  steps <- vapply(4:1, function(r) hexify:::step_aperture(coarser(g, 4L - r)), 1L)
  expect_equal(steps, c(4L, 3L, 7L, 4L))
  # ISEA43H at resolution 4 is 4,4,3,3 and at 3 is 4,3,3: the step is a 4
  expect_equal(hexify:::step_aperture(hex_grid(resolution = 4, aperture = "4/3")), 4L)
  expect_equal(hexify:::step_aperture(hex_grid(resolution = 5, aperture = "4/3")), 3L)
  for (r in 4:1) {
    gr <- coarser(g, 4L - r)
    ov <- get_parent(all_cells(gr), gr, overlapping = TRUE)
    n <- sort(unique(lengths(ov)))
    expect_equal(n, switch(as.character(steps[5L - r]), `3` = c(1L, 3L), c(1L, 2L)))
  }
})

test_that("Fuller and H3 shares are measured on the sphere", {
  skip_on_cran()
  for (g in list(hex_grid(resolution = 3, aperture = 7, projection = "fuller"),
                 hex_grid(resolution = 3, aperture = 4, projection = "fuller"),
                 hex_grid(resolution = 2, type = "h3"))) {
    ids <- all_cells(g)
    area <- unname(cell_area(ids, g))
    up <- hex_aggregate(ids, area, g)
    expect_equal(up$value, unname(cell_area(up$cell_id, coarser(g))), tolerance = 1e-8)
    expect_equal(sum(up$value), sum(area), tolerance = 1e-12)
  }
  # Fuller's shares depart from the plane's: it is not equal-area
  g <- hex_grid(resolution = 3, aperture = 7, projection = "fuller")
  ids <- all_cells(g)[1:200]
  ov <- get_parent(ids, g, overlapping = TRUE)
  measured <- hexify:::sphere_shares(ids, ov, g)
  expect_gt(max(abs(measured - hexify:::lattice_shares(lengths(ov), 7))), 1e-4)
})

test_that("means, missing values and several value columns", {
  g <- hex_grid(resolution = 3, aperture = 4)
  ids <- all_cells(g)
  set.seed(7)
  v <- data.frame(a = runif(length(ids)), b = rnorm(length(ids)))
  up <- hex_aggregate(ids, v, g)
  expect_named(up, c("cell_id", "a", "b"))
  expect_equal(up$a, hex_aggregate(ids, v$a, g)$value)
  expect_s3_class(up$cell_id, "integer64")
  expect_true(all(diff(as.numeric(up$cell_id)) > 0))

  # A mean lies between the smallest and largest value it averages
  m <- hex_aggregate(ids, v$a, g, measure = "mean")
  expect_true(all(m$value >= min(v$a) & m$value <= max(v$a)))

  # NA propagates unless na.rm; with na.rm a missing cell is an absent one
  x <- v$a
  x[1] <- NA
  expect_true(anyNA(hex_aggregate(ids, x, g)$value))
  kept <- hex_aggregate(ids, x, g, na.rm = TRUE)
  dropped <- hex_aggregate(ids[-1], x[-1], g)
  expect_equal(kept$value[match(dropped$cell_id, kept$cell_id)], dropped$value)
  mk <- hex_aggregate(ids, x, g, measure = "mean", na.rm = TRUE)
  md <- hex_aggregate(ids[-1], x[-1], g, measure = "mean")
  expect_equal(mk$value[match(md$cell_id, mk$cell_id)], md$value)
  all_na <- hex_aggregate(ids[1:3], rep(NA_real_, 3), g, measure = "mean", na.rm = TRUE)
  expect_true(all(is.na(all_na$value)))

  expect_error(hex_aggregate(ids, v$a[-1], g), "one entry")
  expect_error(hex_aggregate(c(ids[1], ids[1]), 1:2, g), "once")
  expect_error(hex_aggregate(ids, v$a, g, levels = 4), "levels")
  expect_error(hex_aggregate(ids, as.character(v$a), g), "numeric")
})

test_that("a Hex9 cell on its parents' edge splits in halves", {
  # Six of nine children inside their parent, three centred on its edges and
  # cut there into their two half-hexagons; the parent holding the mode-0
  # half comes first. Snyder's and IVEA's shares are the lattice's for every
  # cell, the twelve beside the vertices among them; Kaseorg's projection
  # ("ak") is measured, and so is Hex9's own ("akw").
  for (projection in c("isea", "ivea")) {
    g <- hex_grid(resolution = 3, aperture = 9, projection = projection)
    ids <- all_cells(g)
    ov <- get_parent(ids, g, overlapping = TRUE)
    expect_equal(as.vector(table(lengths(ov))), length(ids) * c(2, 1) / 3)
    expect_identical(hexify:::cell_id_unlist(lapply(ov, function(x) x[1])), get_parent(ids, g))
    lattice <- hexify:::lattice_shares(lengths(ov), 9L)
    expect_lt(max(abs(hexify:::sphere_shares(ids, ov, g) - lattice)), 1e-10)
    density <- unname(cell_area(ids, g))
    for (lv in 1:2) {
      up <- hex_aggregate(ids, density, g, levels = lv)
      expect_equal(up$value, unname(cell_area(up$cell_id, coarser(g, lv))), tolerance = 1e-12)
      v <- seq_along(ids) %% 5
      for (rule in c("area", "centre")) {
        expect_equal(sum(hex_aggregate(ids, v, g, levels = lv, rule = rule)$value), sum(v),
                     tolerance = 1e-12)
      }
    }
  }
  projections <- "ak"
  if (hexify:::cpp_hex9_warp_ready() || file.exists(hexify:::hex9_warp_path())) {
    projections <- c(projections, "akw")
  }
  for (projection in projections) {
    g <- hex_grid(resolution = 2, aperture = 9, projection = projection)
    ids <- all_cells(g)
    area <- unname(cell_area(ids, g))
    up <- hex_aggregate(ids, area, g)
    expect_equal(up$value, unname(cell_area(up$cell_id, coarser(g))), tolerance = 1e-6)
    expect_equal(sum(up$value), sum(area), tolerance = 1e-12)
  }
})
