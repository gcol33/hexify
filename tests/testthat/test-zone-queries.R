# Zone queries of OGC Topic 21: overlapping parents, siblings, hollow rings,
# and the commutation of binning with the hierarchy.

# Share of a cell lying in each coarser cell it overlaps, measured in a
# Lambert azimuthal equal-area plane centred on the cell
overlap_shares <- function(child, parents, g, pg) {
  ctr <- cell_to_lonlat(child, g)
  laea <- sprintf("+proj=laea +lon_0=%.10f +lat_0=%.10f +R=6371007.18", ctr$lon_deg, ctr$lat_deg)
  plane <- function(ids, grid) {
    sf::st_geometry(sf::st_transform(cell_to_sf(ids, grid, wrap_dateline = FALSE), laea))
  }
  cg <- plane(child, g)
  pp <- plane(parents, pg)
  a <- as.numeric(sf::st_area(cg))
  vapply(seq_along(parents), function(k) {
    as.numeric(sf::st_area(sf::st_intersection(cg, pp[k]))) / a
  }, numeric(1))
}

first_of <- function(ov) cell_id_unlist(lapply(ov, function(v) v[1]))

test_that("overlapping parents: three, two or one, starting with the parent", {
  for (ap in c(3, 4, 7)) {
    res <- c(`3` = 5, `4` = 4, `7` = 3)[[as.character(ap)]]
    g <- hex_grid(resolution = res, aperture = ap)
    ids <- grid_cells(g)
    ov <- get_parent(ids, g, overlapping = TRUE)
    expect_equal(first_of(ov), get_parent(ids, g))
    n <- lengths(ov)[!is_pentagon(ids, g)]
    # Per parent: aperture 3 one centred child and six on corners, each a third
    # in three parents; aperture 4 one centred and six on edges, each half in
    # two; aperture 7 one centred and six reaching 1/12 into a neighbour
    expected <- switch(as.character(ap), `3` = c(1, 3), `4` = c(1, 2), `7` = c(1, 2))
    expect_setequal(unique(n), expected)
    share <- mean(n == expected[2])
    expect_equal(share, switch(as.character(ap), `3` = 2 / 3, `4` = 3 / 4, `7` = 6 / 7),
                 tolerance = 0.02)
  }
})

test_that("each overlapping parent holds the share of the cell the lattice gives", {
  cases <- list(list(3, 5, c(1, 3) / 3), list(4, 4, c(1, 1) / 2), list(7, 3, c(1, 11) / 12))
  for (case in cases) {
    g <- hex_grid(resolution = case[[2]], aperture = case[[1]])
    pg <- hex_grid(resolution = case[[2]] - 1, aperture = case[[1]])
    cells <- lonlat_to_cell(c(10, 20, 30, 40, 50) + 0.37, c(15, 25, 35, -20, -40), g)
    ov <- get_parent(cells, g, overlapping = TRUE)
    for (k in seq_along(cells)) {
      if (length(ov[[k]]) == 1L) next
      share <- overlap_shares(cells[k], ov[[k]], g, pg)
      expect_equal(sum(share), 1, tolerance = 1e-3)
      expect_true(all(abs(share - case[[3]][1]) < 0.01 | abs(share - case[[3]][2]) < 0.01))
    }
  }
})

test_that("overlapping parents several levels up and on H3", {
  g <- hex_grid(resolution = 5, aperture = 7)
  ids <- lonlat_to_cell(c(5, 60, -120), c(40, -10, 70), g)
  ov2 <- get_parent(ids, g, levels = 2, overlapping = TRUE)
  expect_equal(first_of(ov2), get_parent(ids, g, levels = 2))
  # every cell that overlaps a grandparent's region overlaps one of the parents'
  ov1 <- get_parent(ids, g, overlapping = TRUE)
  g4 <- hex_grid(resolution = 4, aperture = 7)
  for (k in seq_along(ids)) {
    up <- unique(cell_id_unlist(get_parent(ov1[[k]], g4, overlapping = TRUE)))
    expect_true(all(ov2[[k]] %in% up))
  }

  h <- hex_grid(resolution = 6, type = "h3")
  hid <- lonlat_to_cell(c(5, 60), c(40, -10), h)
  hov <- get_parent(hid, h, overlapping = TRUE)
  expect_equal(vapply(hov, `[`, "", 1), get_parent(hid, h))
  expect_true(all(lengths(hov) %in% c(1L, 2L, 3L)))
})

test_that("siblings are the other children of the parent", {
  for (g in list(hex_grid(resolution = 4, aperture = 7),
                 hex_grid(resolution = 5, aperture = 3),
                 hex_grid(resolution = 5, type = "h3"))) {
    cells <- lonlat_to_cell(c(16.37, -70, 120), c(48.21, -30, 10), g)
    sib <- get_siblings(cells, g)
    all_kids <- get_siblings(cells, g, include_self = TRUE)
    pg <- hexify:::grid_at_resolution(g, g@resolution - 1L)
    kids <- get_children(get_parent(cells, g), pg)
    for (k in seq_along(cells)) {
      expect_false(cells[k] %in% sib[[k]])
      expect_setequal(c(sib[[k]], cells[k]), kids[[k]])
      expect_setequal(all_kids[[k]], kids[[k]])
    }
  }
  g7 <- hex_grid(resolution = 4, aperture = 7)
  expect_equal(lengths(get_siblings(lonlat_to_cell(16.37, 48.21, g7), g7)), 6L)
})

test_that("ring = TRUE keeps only the cells k hops away", {
  g <- hex_grid(resolution = 5, aperture = 4)
  cell <- lonlat_to_cell(16.37, 48.21, g)
  disk <- get_neighbors(cell, g, k = 3, distances = TRUE)[[1]]
  ring <- get_neighbors(cell, g, k = 3, ring = TRUE)[[1]]
  expect_setequal(ring, disk$cell_id[disk$ring_distance == 3])
  expect_length(ring, 18)
  rd <- get_neighbors(cell, g, k = 2, ring = TRUE, distances = TRUE)[[1]]
  expect_true(all(rd$ring_distance == 2))
  expect_equal(get_neighbors(cell, g, k = 0, ring = TRUE)[[1]], cell)
  links <- get_neighbors(cell, g, k = 2, ring = TRUE, as_sf = TRUE)
  expect_true(all(links$ring_distance == 2))

  h <- hex_grid(resolution = 6, type = "h3")
  hc <- lonlat_to_cell(16.37, 48.21, h)
  expect_length(get_neighbors(hc, h, k = 2, ring = TRUE)[[1]], 12)
})

test_that("binning then taking the parent disagrees with binning at the parent's resolution as the lattice predicts", {
  # A uniform point lands in a cell that reaches out of its parent with the
  # probability the lattice gives: aperture 3, 2/3 of the area lies in corner
  # cells and 2/3 of each lies outside the parent it is assigned to, 4/9;
  # aperture 4, 3/4 in edge cells half outside, 3/8; aperture 7, 6/7 in cells
  # 1/12 outside, 1/14. Snyder's projection is equal-area, so the shares on the
  # sphere are those of the plane, up to the cells beside the twelve vertices.
  pts <- sphere_test_points(10000)
  rate <- function(g, pg) {
    fine <- lonlat_to_cell(pts$lon, pts$lat, g)
    mean(get_parent(fine, g) != lonlat_to_cell(pts$lon, pts$lat, pg))
  }
  expect_equal(rate(hex_grid(resolution = 6, aperture = 3),
                    hex_grid(resolution = 5, aperture = 3)), 4 / 9, tolerance = 0.03)
  expect_equal(rate(hex_grid(resolution = 5, aperture = 4),
                    hex_grid(resolution = 4, aperture = 4)), 3 / 8, tolerance = 0.03)
  expect_equal(rate(hex_grid(resolution = 4, aperture = 7),
                    hex_grid(resolution = 3, aperture = 7)), 1 / 14, tolerance = 0.05)
  # A mixed sequence: the last step is aperture 3
  expect_equal(rate(hex_grid(resolution = 5, aperture = c(4, 7, 4, 4, 3)),
                    hex_grid(resolution = 4, aperture = c(4, 7, 4, 4))), 4 / 9,
               tolerance = 0.03)
})
