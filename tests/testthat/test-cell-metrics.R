# test-cell-metrics.R
# Tests for cell_metrics() and wall_metrics()

gc_rad <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * asin(pmin(1, sqrt(a)))
}

ring_length_rad <- function(m) {
  n <- nrow(m)
  sum(gc_rad(m[-n, 1], m[-n, 2], m[-1, 1], m[-1, 2]))
}

test_that("an aperture-3 grid has 10 * 3^r + 2 cells (Carr et al. 1997, Table 1)", {
  for (res in 0:4) {
    m <- cell_metrics(grid = hex_grid(resolution = res, aperture = 3))
    expect_equal(nrow(m), 10 * 3^res + 2)
  }
})

test_that("every wall of a grid is measured once", {
  for (g in list(hex_grid(resolution = 3, aperture = 3),
                 hex_grid(resolution = 2, aperture = 7),
                 hex_grid(resolution = 3, aperture = 4, projection = "fuller"),
                 hex_grid(resolution = 3, aperture = "4/3"),
                 hex_grid(resolution = 3, aperture = 3, polyhedron = "octahedron"),
                 hex_grid(resolution = 3, aperture = 4, polyhedron = "tetrahedron"),
                 hex_grid(resolution = 0, aperture = 3),
                 hex_grid(resolution = 1, type = "h3"))) {
    w <- wall_metrics(grid = g)
    n <- grid_n_cells(g)
    # Every vertex cell lacks as many walls as its solid's angular deficit in
    # sixths, and these add up to 12 on any solid with triangular faces.
    expect_equal(nrow(w), (6 * n - 12) / 2)
    expect_false(anyDuplicated(paste(pmin(w$cell_id, w$neighbor_id),
                                     pmax(w$cell_id, w$neighbor_id))) > 0)
  }
})

test_that("a wall measures the same from either side", {
  g <- hex_grid(resolution = 3, aperture = 7, projection = "fuller")
  w <- cell_wall_measures(as_cell_id(seq_len(grid_n_cells(g))), g, walls = TRUE)$walls
  nbr <- as.integer(w$neighbor_id)
  key <- paste(pmin(w$cell, nbr), pmax(w$cell, nbr))
  expect_true(all(table(key) == 2L))
  spread <- function(x) tapply(x, key, function(v) diff(range(v)))
  expect_lt(max(spread(w$wall) / tapply(w$wall, key, mean)), 1e-9)
  expect_lt(max(spread(w$midpoint_offset)), 1e-9)
  expect_lt(max(spread(w$centre_distance)), 1e-12)
})

test_that("perimeters add up to twice the walls", {
  g <- hex_grid(resolution = 2, aperture = 4)
  m <- cell_metrics(grid = g)
  w <- wall_metrics(grid = g)
  expect_equal(sum(m$perimeter_km), 2 * sum(w$wall_km), tolerance = 1e-12)
})

test_that("an ISEA perimeter matches its boundary densified in lon/lat", {
  # Two densifications of the same curved walls: on the sphere (cell_metrics)
  # and in lon/lat (cell_to_sf's boundary), each to a tight tolerance. The
  # lon/lat walk stops at 2^12 pieces per edge without cutting it at face
  # edges, so where a wall bends sharply across a face edge, next to a vertex
  # of the solid, it cuts the bend and comes out short by up to about 1e-5.
  for (proj in c("isea", "fuller")) {
    g <- hex_grid(resolution = 3, aperture = 7, projection = proj)
    ids <- as_cell_id(c(1, 2, 100, 500, 1000, grid_n_cells(g)))
    lv <- isea_levels(g@aperture, g@resolution)
    rings <- cpp_cell_to_corners(icosa_arg(g), ids, lv$resolution, lv$aperture,
                                 lv$ap_seq, 1e-7)
    corners <- cpp_cell_to_corners(icosa_arg(g), ids, lv$resolution, lv$aperture,
                                   lv$ap_seq, 0)
    p <- cell_metrics(ids, g)$perimeter_km / grid_radius_km(g)
    expect_equal(p, vapply(rings, ring_length_rad, 0), tolerance = 1e-5)
    # A curved wall is longer than the great-circle arc between its corners
    expect_true(all(p > vapply(corners, ring_length_rad, 0)))
  }
})

test_that("an H3 perimeter is the great-circle length of its boundary", {
  g <- hex_grid(resolution = 4, type = "h3")
  ids <- lonlat_to_cell(c(0, 30, -120, 10.5), c(0, 60, -45, 64.7), g)
  ids <- c(ids, h3_all_cells(0L)[1:3])
  p <- cell_metrics(ids, hex_grid(resolution = 4, type = "h3"))$perimeter_km
  expected <- vapply(cpp_h3_cellToBoundary(ids), ring_length_rad, 0) *
    grid_radius_km(g)
  expect_equal(p, expected, tolerance = 1e-12)
})

test_that("compactness is at most that of the regular polygon", {
  g <- hex_grid(resolution = 6, aperture = 4)
  vertex <- c(1, 2, grid_n_cells(g))
  ids <- c(as_cell_id(vertex), lonlat_to_cell(c(0, 45, 100), c(10, 30, -60), g))
  m <- cell_metrics(ids, g)
  # A small regular k-gon scores 2 sqrt(pi A) / P, up to the sphere's
  # curvature across the cell, about 1e-4 here. Each wall of a vertex cell
  # bends where it crosses a face edge, so the vertex cells, alike by the
  # solid's symmetry, fall below the regular pentagon.
  regular <- function(k) 2 * sqrt(pi * k / 4 / tan(pi / k)) / k
  pent <- m$compactness[seq_along(vertex)]
  expect_equal(pent, rep(pent[1], length(vertex)), tolerance = 1e-9)
  expect_lt(pent[1], regular(5))
  hex <- m$compactness[-seq_along(vertex)]
  expect_true(all(hex < regular(6) + 1e-3 & hex > 0.9))
})

test_that("normalized area and ipq follow from area, perimeter and cell count", {
  g <- hex_grid(resolution = 3, aperture = 7)
  m <- cell_metrics(grid = g)
  # Snyder's projection is equal-area: the twelve pentagons have 5/6 of a
  # hexagon's area, so a hexagon has N / (N - 2) of the mean
  n <- nrow(m)
  pent <- is_pentagon(m$cell_id, g)
  expect_equal(sum(pent), 12)
  expect_equal(mean(m$normalized_area), 1, tolerance = 1e-12)
  expect_equal(range(m$normalized_area[!pent]), rep(n / (n - 2), 2), tolerance = 1e-12)
  expect_equal(range(m$normalized_area[pent]), rep(5 / 6 * n / (n - 2), 2),
               tolerance = 1e-12)
  # ipq = compactness^2 + (a / (r p))^2
  r <- grid_radius_km(g)
  omega <- m$area_km2 / hexify:::body_surface_km2(r) * 4 * pi
  expect_equal(m$ipq, m$compactness^2 + (omega / (m$perimeter_km / r))^2,
               tolerance = 1e-12)
})

test_that("a grid on another body scales its lengths and keeps its shapes", {
  earth <- hex_grid(resolution = 2, aperture = 3)
  mars <- hex_grid(resolution = 2, aperture = 3, radius_km = "mars")
  ratio <- grid_radius_km(mars) / grid_radius_km(earth)
  me <- cell_metrics(grid = earth)
  mm <- cell_metrics(grid = mars)
  expect_equal(mm$perimeter_km, me$perimeter_km * ratio, tolerance = 1e-12)
  expect_equal(mm$compactness, me$compactness, tolerance = 1e-9)
  we <- wall_metrics(grid = earth)
  wm <- wall_metrics(grid = mars)
  expect_equal(wm$center_distance_km, we$center_distance_km * ratio, tolerance = 1e-12)
  expect_equal(wm$midpoint_ratio, we$midpoint_ratio, tolerance = 1e-12)
})

test_that("walls of resolution 0 have midpoint ratio 0 by symmetry", {
  # The cells are the solid's vertex cells, and the reflection swapping two
  # adjacent vertices maps the grid onto itself and their wall onto itself
  # end for end, so the arc between the centres crosses it at its middle.
  g <- hex_grid(resolution = 0, aperture = 3)
  w <- wall_metrics(grid = g)
  expect_true(all(w$midpoint_ratio < 1e-9))
})

test_that("corner rings give the same walls as the grid's own geometry", {
  g <- hex_grid(resolution = 3, aperture = 3)
  ids <- as_cell_id(seq_len(grid_n_cells(g)))
  lv <- isea_levels(g@aperture, g@resolution)
  corners <- cpp_cell_to_corners(icosa_arg(g), ids, lv$resolution, lv$aperture,
                                 lv$ap_seq, 0)
  ctr <- cell_to_lonlat(ids, g)
  centres <- cbind(ctr$lon_deg, ctr$lat_deg)
  nb <- grid_neighbors_isea(ids, g)
  rw <- cpp_ring_walls(corners, centres, lapply(nb, function(x) centres[as.integer(x), , drop = FALSE]))
  own <- cell_wall_measures(ids, g, walls = TRUE)$walls
  expect_equal(nrow(rw), nrow(own))
  # Rows run cell by cell and neighbour by neighbour in both
  expect_equal(rw$cell, own$cell)
  expect_equal(rw$centre_distance, own$centre_distance, tolerance = 1e-9)
  # A wall that bends across a face edge is longer than the arc between its
  # corners, by up to a few percent at this resolution
  expect_true(all(own$wall >= rw$wall - 1e-12))
  expect_equal(rw$wall, own$wall, tolerance = 0.05)
})

test_that("metrics read a HexData object's cells and keep duplicates in place", {
  g <- hex_grid(area_km2 = 5000)
  df <- data.frame(lon = c(0, 0, 10), lat = c(45, 45, 50))
  hd <- hexify(df, lon = "lon", lat = "lat", grid = g)
  m <- cell_metrics(grid = hd)
  expect_equal(nrow(m), 3)
  expect_equal(m[1, ], m[2, ], ignore_attr = TRUE)
  w <- wall_metrics(grid = hd)
  expect_equal(sort(unique(w$cell_id)), sort(unique(hd@cell_id)))
})

test_that("a NA cell gives NA metrics and no walls", {
  g <- hex_grid(area_km2 = 5000)
  m <- cell_metrics(c(NA, 5), g)
  expect_true(all(is.na(m[1, -1])))
  expect_false(anyNA(m[2, ]))
  expect_equal(nrow(wall_metrics(c(NA, 5), g)), 6L)
})
