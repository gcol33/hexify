# Aperture-7 neighbours by IGEO7 digit arithmetic, against the neighbours the
# quad lattice gives with every quad-leaving step sent through lon/lat.

library(testthat)

z7_and_lonlat_neighbours <- function(g, ids) {
  lv <- hexify:::isea_levels(g@aperture, g@resolution)
  icosa <- hexify:::icosa_arg(g)
  list(z7 = hexify:::cpp_get_neighbors_isea(icosa, ids, lv$resolution, lv$aperture,
                                            lv$ap_seq),
       lonlat = hexify:::cpp_get_neighbors_isea_lonlat(icosa, ids, lv$resolution,
                                                       lv$aperture, lv$ap_seq))
}

expect_same_neighbours <- function(nb, info) {
  expect_identical(lengths(nb$z7), lengths(nb$lonlat), info = info)
  expect_identical(cell_id_unlist(nb$z7), cell_id_unlist(nb$lonlat), info = info)
}

test_that("digit arithmetic gives every aperture-7 cell's neighbours, resolutions 0-5", {
  for (res in 0:5) {
    g <- hex_grid(resolution = res, aperture = 7)
    ids <- as_cell_id(seq_len(n_cells(g)))
    expect_same_neighbours(z7_and_lonlat_neighbours(g, ids), sprintf("resolution %d", res))
  }
})

test_that("digit arithmetic gives every aperture-7 cell's neighbours at resolution 6", {
  skip_on_cran()
  g <- hex_grid(resolution = 6, aperture = 7)
  ids <- as_cell_id(seq_len(n_cells(g)))
  expect_same_neighbours(z7_and_lonlat_neighbours(g, ids), "resolution 6")
})

test_that("digit arithmetic follows the topology under any orientation and projection", {
  set.seed(103)
  grids <- list(
    hex_grid(resolution = 3, aperture = 7, orientation = c(-40, 20, 33)),
    hex_grid(resolution = 4, aperture = 7, projection = "fuller"),
    hex_grid(resolution = 3, aperture = 7, projection = "ivea"))
  for (g in grids) {
    ids <- as_cell_id(seq_len(n_cells(g)))
    expect_same_neighbours(z7_and_lonlat_neighbours(g, ids),
                           paste(grid_projection(g), g@resolution))
  }
  # resolution 9, cells near quad edges and vertices: the neighbours of the
  # cells of random points and of the twelve pentagons
  g <- hex_grid(resolution = 9, aperture = 7)
  lon <- runif(2000, -180, 180)
  lat <- asin(runif(2000, -1, 1)) * 180 / pi
  pent <- index_to_cell(sprintf("%02d%s", 0:11, strrep("0", 9)), g)
  ids <- c(lonlat_to_cell(lon, lat, g), pent)
  ring <- cell_id_unlist(get_neighbors(ids, g))
  expect_same_neighbours(z7_and_lonlat_neighbours(g, unique(c(ids, ring))), "resolution 9")
})

test_that("the twelve pentagons have five neighbours and every other cell six", {
  for (res in c(1L, 4L)) {
    g <- hex_grid(resolution = res, aperture = 7)
    n <- lengths(get_neighbors(as_cell_id(seq_len(n_cells(g))), g))
    pent <- index_to_cell(sprintf("%02d%s", 0:11, strrep("0", res)), g)
    expect_identical(sort(as.integer(pent)), which(n == 5L), info = sprintf("resolution %d", res))
    expect_true(all(n[-as.integer(pent)] == 6L), info = sprintf("resolution %d", res))
  }
})

test_that("adjacency from digit arithmetic is symmetric", {
  g <- hex_grid(resolution = 5, aperture = 7)
  ids <- as_cell_id(seq_len(n_cells(g)))
  nb <- get_neighbors(ids, g)
  from <- rep(ids, lengths(nb))
  to <- cell_id_unlist(nb)
  key <- function(a, b) paste(as.character(a), as.character(b))
  expect_true(all(key(to, from) %in% key(from, to)))
})
