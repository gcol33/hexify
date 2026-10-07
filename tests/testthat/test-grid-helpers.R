# tests/testthat/test-grid-helpers.R
# Tests for grid helper functions
# Note: grid_global, grid_rect, grid_clip tests removed for CRAN speed

test_that("lonlat_to_cell works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
  expect_true(all(cells > 0))
})

test_that("lonlat_to_cell works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Use the HexData object as grid source
  cells <- lonlat_to_cell(lon = 5, lat = 48, grid = result)

  expect_type(cells, "double")
  expect_length(cells, 1)
})

test_that("lonlat_to_cell works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

  expect_type(cells, "double")
  expect_length(cells, 2)
})

test_that("cell_to_lonlat returns cell centers", {
  grid <- hex_grid(area_km2 = 1000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_true("lon_deg" %in% names(coords))
  expect_true("lat_deg" %in% names(coords))
  expect_length(coords$lon_deg, 2)
})

test_that("cell_to_lonlat works with mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  coords <- cell_to_lonlat(cells, grid)

  expect_type(coords, "list")
  expect_length(coords$lon_deg, 2)
})

test_that("lonlat_to_cell -> cell_to_lonlat round-trip lands in same cell", {
  grid <- hex_grid(area_km2 = 1000)

  original_lon <- c(0, 10, -5)
  original_lat <- c(45, 50, 48)

  cells1 <- lonlat_to_cell(original_lon, original_lat, grid)
  centers <- cell_to_lonlat(cells1, grid)
  cells2 <- lonlat_to_cell(centers$lon_deg, centers$lat_deg, grid)

  expect_equal(cells1, cells2)
})

test_that("cell_to_sf creates sf polygons", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)
  polys <- cell_to_sf(cells, grid)

  expect_s3_class(polys, "sf")
  expect_true("cell_id" %in% names(polys))
  expect_equal(nrow(polys), length(unique(cells)))
})

test_that("cell_to_sf works with HexData (no cell_id)", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  polys <- cell_to_sf(grid = result)

  expect_s3_class(polys, "sf")
  expect_equal(nrow(polys), length(unique(result@cell_id)))
})

test_that("cell_to_sf errors without cell_id for HexGridInfo", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(grid = grid), "cell_id required")
})

test_that("cell_to_sf errors on empty cell_id", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  expect_error(cell_to_sf(cell_id = numeric(0), grid = grid), "No valid")
  expect_error(cell_to_sf(cell_id = c(NA, NA), grid = grid), "No valid")
})

test_that("extract_grid works with HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
  g <- hexify:::extract_grid(result)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid works with legacy hexify_grid", {
  grid <- hexify_grid(area = 1000, aperture = 3)
  g <- hexify:::extract_grid(grid)

  expect_s4_class(g, "HexGridInfo")
})

test_that("extract_grid errors on invalid input", {
  expect_error(hexify:::extract_grid(list(a = 1)), "Cannot extract grid")
  expect_error(hexify:::extract_grid(data.frame()), "Cannot extract grid")
})

test_that("cell_to_sf returns valid geometries for all cells", {
  skip_if_not_installed("sf")

  # Test that polar cells (at icosahedral vertices) have valid geometries
  # Cell 1 is always quad 0 (north pole), and the last cell is quad 11 (south pole)
  grid <- hex_grid(area_km2 = 100000)
  n_cells <- 10 * (as.integer(grid@aperture)^grid@resolution) + 2

  # Test north and south pole cells
  polar_cells <- c(1, n_cells)
  polys <- cell_to_sf(polar_cells, grid)

  # All geometries must be valid
  validity <- sf::st_is_valid(polys)
  expect_true(all(validity))

  # These pentagon cells have 5 corners (6 coords with closing)
  rings <- cpp_cell_to_corners(numeric(0), as_cell_id(polar_cells), grid@resolution,
                               as.integer(grid@aperture), integer(0))
  expect_true(all(vapply(rings, nrow, integer(1)) == 6L))
})

# =============================================================================
# ANTIMERIDIAN HANDLING
# =============================================================================

test_that("cell_to_sf splits ISEA cells at antimeridian", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 100000)
  # Get a cell near the antimeridian
  cell <- lonlat_to_cell(179.5, 0, grid)
  polys <- cell_to_sf(cell, grid)

  # All coordinates must be within [-180, 180]
  coords <- sf::st_coordinates(polys)
  expect_true(all(coords[, "X"] >= -180 & coords[, "X"] <= 180))

  # Geometry must be valid
  expect_true(all(sf::st_is_valid(polys)))
})

test_that("cell_to_sf splits H3 cells at antimeridian", {
  skip_if_not_installed("sf")

  h3 <- hex_grid(resolution = 2, type = "h3")
  cell <- lonlat_to_cell(179.5, 0, h3)
  polys <- cell_to_sf(cell, h3)

  # All coordinates must be within [-180, 180]
  coords <- sf::st_coordinates(polys)
  expect_true(all(coords[, "X"] >= -180 & coords[, "X"] <= 180))

  # Geometry must be valid
  expect_true(all(sf::st_is_valid(polys)))
})

test_that("cell_to_sf keeps non-crossing cells as POLYGON", {
  skip_if_not_installed("sf")

  grid <- hex_grid(area_km2 = 10000)
  # Cell far from antimeridian
  cell <- lonlat_to_cell(10, 45, grid)
  polys <- cell_to_sf(cell, grid)

  geom_type <- as.character(sf::st_geometry_type(polys))
  expect_equal(geom_type, "POLYGON")
})

test_that("grid_rect returns every ISEA cell meeting the box", {
  skip_if_not_installed("sf")
  cases <- list(
    list(c(-74.3, 40.5, -73.7, 40.95), 3, 11),
    list(c(-74.3, 40.5, -73.7, 40.95), 4, 8),
    list(c(-74.3, 40.5, -73.7, 40.95), 7, 5),
    list(c(5, 45, 16, 55), 7, 4),
    list(c(-10, 35, 30, 60), 7, 5),
    list(c(-20, 60, 40, 75), 4, 6)
  )
  for (cs in cases) {
    bbox <- cs[[1]]
    g <- hex_grid(type = "isea", resolution = cs[[3]], aperture = cs[[2]])
    got <- grid_rect(bbox, g)
    step <- g@diagonal_km / 111.32 / 20
    pts <- expand.grid(lon = c(seq(bbox[1], bbox[3], by = step), bbox[3]),
                       lat = c(seq(bbox[2], bbox[4], by = step), bbox[4]))
    ref <- unique(lonlat_to_cell(pts$lon, pts$lat, g))
    expect_length(setdiff(ref, got$cell_id), 0L)
    expect_false(anyDuplicated(got$cell_id) > 0L)
  }
})

# =============================================================================
# Densified cell edges, neighbour links and child polygons
# =============================================================================

ring_of <- function(x) sf::st_coordinates(sf::st_geometry(x)[[1]])[, 1:2]

test_that("cell_to_sf densifies ISEA edges to the tolerance asked", {
  g <- hex_grid(resolution = 3, aperture = 3)
  cell <- lonlat_to_cell(10, 50, g)
  corners <- nrow(ring_of(cell_to_sf(cell, g, densify = 0)))
  default <- nrow(ring_of(cell_to_sf(cell, g)))
  fine <- nrow(ring_of(cell_to_sf(cell, g, densify = 1e-5)))
  expect_equal(corners, 7L)
  expect_gt(default, corners)
  expect_gt(fine, default)
  expect_equal(ring_of(cell_to_sf(cell, g, densify = 0.001)),
               ring_of(cell_to_sf(cell, g)))
  expect_error(cell_to_sf(cell, g, densify = -1), "non-negative")
})

test_that("cell_to_sf densifies H3 edges along their great circles", {
  g <- suppressMessages(hex_grid(resolution = 2, type = "h3"))
  cell <- lonlat_to_cell(10, 50, g)
  corners <- ring_of(cell_to_sf(cell, g, wrap_dateline = FALSE))
  dense <- ring_of(cell_to_sf(cell, g, wrap_dateline = FALSE, densify = 1e-4))
  expect_equal(nrow(corners), 7L)
  expect_gt(nrow(dense), nrow(corners))
  V <- hexify:::unit_vec(corners[, 1], corners[, 2])
  normals <- t(vapply(1:6, function(i) {
    n <- c(V[i, 2] * V[i + 1, 3] - V[i, 3] * V[i + 1, 2],
           V[i, 3] * V[i + 1, 1] - V[i, 1] * V[i + 1, 3],
           V[i, 1] * V[i + 1, 2] - V[i, 2] * V[i + 1, 1])
    n / sqrt(sum(n^2))
  }, numeric(3)))
  D <- hexify:::unit_vec(dense[, 1], dense[, 2])
  off <- apply(abs(D %*% t(normals)), 1, min)
  expect_lt(max(off), 1e-10)
})

piece_km <- function(ring, radius_km) {
  p <- ring[-nrow(ring), , drop = FALSE] * pi / 180
  q <- ring[-1L, , drop = FALSE] * pi / 180
  h <- sin((q[, 2] - p[, 2]) / 2)^2 +
    cos(p[, 2]) * cos(q[, 2]) * sin((q[, 1] - p[, 1]) / 2)^2
  2 * asin(sqrt(pmin(1, h))) * radius_km
}

test_that("cell_to_sf(max_km) bounds every drawn piece in km", {
  g <- hex_grid(resolution = 3, aperture = 3)
  # Cell 82 has an 880 km edge the relative tolerance leaves as one piece
  default <- ring_of(cell_to_sf(82, g, wrap_dateline = FALSE))
  capped <- ring_of(cell_to_sf(82, g, wrap_dateline = FALSE, max_km = 50))
  expect_gt(max(piece_km(default, grid_radius_km(g))), 800)
  expect_lte(max(piece_km(capped, grid_radius_km(g))), 50)
  # The cap only halves pieces further, so every default point stays
  key <- function(r) paste(signif(r[, 1], 12), signif(r[, 2], 12))
  expect_true(all(key(default) %in% key(capped)))

  h <- suppressMessages(hex_grid(resolution = 2, type = "h3"))
  hc <- lonlat_to_cell(10, 50, h)
  expect_lte(max(piece_km(ring_of(cell_to_sf(hc, h, wrap_dateline = FALSE,
                                             max_km = 20)), grid_radius_km(h))), 20)

  g7 <- hex_grid(resolution = 2, aperture = 7)
  c7 <- lonlat_to_cell(10, 50, g7)
  for (shape in c("gosper", "descendants")) {
    ring <- ring_of(cell_to_sf(c7, g7, shape = shape, depth = 1,
                               wrap_dateline = FALSE, max_km = 20))
    expect_lte(max(piece_km(ring, grid_radius_km(g7))), 20)
  }

  expect_error(cell_to_sf(82, g, max_km = 0), "positive")
  expect_error(cell_to_sf(82, g, max_km = c(1, 2)), "positive")
})

test_that("cell_to_sf(max_km) reads km on the grid's own sphere", {
  earth <- hex_grid(resolution = 3, aperture = 3)
  moon <- hex_grid(resolution = 3, aperture = 3, radius_km = "moon")
  ratio <- grid_radius_km(moon) / grid_radius_km(earth)
  expect_equal(ring_of(cell_to_sf(82, moon, wrap_dateline = FALSE,
                                  max_km = 50 * ratio)),
               ring_of(cell_to_sf(82, earth, wrap_dateline = FALSE,
                                  max_km = 50)))
})

test_that("get_neighbors(as_sf = TRUE) links each cell to its neighbours", {
  g <- hex_grid(resolution = 4, aperture = 4)
  cells <- lonlat_to_cell(c(10, -60), c(50, -15), g)
  links <- get_neighbors(cells, g, k = 2, as_sf = TRUE)
  rings <- get_neighbors(cells, g, k = 2, distances = TRUE)
  expect_s3_class(links, "sf")
  expect_equal(links$cell_id, rep(cells, vapply(rings, nrow, integer(1))))
  expect_equal(links$neighbor_id, cell_id_unlist(lapply(rings, `[[`, "cell_id")))
  expect_equal(links$ring_distance, unlist(lapply(rings, `[[`, "ring_distance")))
  ends <- t(vapply(sf::st_geometry(links), function(l) {
    xy <- sf::st_coordinates(l)
    c(xy[1, 1:2], xy[nrow(xy), 1:2])
  }, numeric(4)))
  from <- cell_to_lonlat(links$cell_id, g)
  to <- cell_to_lonlat(links$neighbor_id, g)
  expect_equal(unname(ends), unname(cbind(from$lon_deg, from$lat_deg,
                                          to$lon_deg, to$lat_deg)),
               tolerance = 1e-9)
})

test_that("get_neighbors(as_sf = TRUE) splits a link at the antimeridian", {
  g <- hex_grid(resolution = 3, aperture = 3)
  cell <- lonlat_to_cell(179.9, 10, g)
  links <- get_neighbors(cell, g, as_sf = TRUE)
  lon <- sf::st_coordinates(links)[, 1]
  expect_true(all(lon >= -180 & lon <= 180))
})

test_that("get_children(as_sf = TRUE) returns each child with its parent", {
  for (ap in list(4, 3, "4/3")) {
    g <- hex_grid(resolution = 2, aperture = ap)
    parents <- c(1, 5, 5)
    kids <- get_children(parents, g, levels = 2, as_sf = TRUE)
    ids <- get_children(parents, g, levels = 2)
    expect_equal(kids$parent_id, rep(as_cell_id(parents), lengths(ids)), info = ap)
    expect_equal(kids$cell_id, cell_id_unlist(ids), info = ap)
    expect_false(any(sf::st_is_empty(kids)), info = ap)
  }
  h <- suppressMessages(hex_grid(resolution = 1, type = "h3"))
  hc <- lonlat_to_cell(c(10, -60), c(50, -15), h)
  hk <- get_children(hc, h, as_sf = TRUE)
  expect_equal(nrow(hk), length(unlist(get_children(hc, h))))
})
