# ISEA cell IDs are 64-bit integers, exact at every resolution a grid takes

finest <- list(
  list(aperture = 3, resolution = 30L),
  list(aperture = 4, resolution = 29L),
  list(aperture = 7, resolution = 21L)
)

random_points <- function(n) {
  set.seed(1)
  list(lon = stats::runif(n, -180, 180),
       lat = asin(stats::runif(n, -1, 1)) * 180 / pi)
}

great_circle_km <- function(lon1, lat1, lon2, lat2, radius_km) {
  r <- pi / 180
  h <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * radius_km * asin(sqrt(h))
}

test_that("as_cell_id() reads integer64, whole numbers and digit strings", {
  big <- bit64::as.integer64("9007199254740993")
  expect_identical(hexify:::as_cell_id(big), big)
  expect_identical(hexify:::as_cell_id("9007199254740993"), big)
  expect_identical(hexify:::as_cell_id(c(1, NA, 12)),
                   bit64::as.integer64(c(1, NA, 12)))
  expect_identical(hexify:::as_cell_id(5L), bit64::as.integer64(5))
  expect_identical(hexify:::as_cell_id(NA), bit64::as.integer64(NA))
})

test_that("as_cell_id() refuses what it cannot read exactly", {
  expect_error(hexify:::as_cell_id(1.5), "whole-number")
  expect_error(hexify:::as_cell_id(2^53), "2\\^53")
  expect_error(hexify:::as_cell_id("12a"), "whole-number")
  expect_error(hexify:::as_cell_id(list(1)), "cell IDs")
})

test_that("ISEA cell IDs are returned as integer64", {
  g <- hex_grid(resolution = 5, aperture = 3)
  id <- lonlat_to_cell(c(16.373, 0), c(48.2, 0), g)
  expect_s3_class(id, "integer64")
  expect_s3_class(get_parent(id, g), "integer64")
  expect_s3_class(get_neighbors(id, g)[[1]], "integer64")
  expect_s3_class(get_children(id, g)[[1]], "integer64")
  expect_s3_class(cell_to_sf(id, g)$cell_id, "integer64")
  expect_s3_class(hexify(data.frame(lon = 16.373, lat = 48.2), lon = "lon",
                         lat = "lat", grid = g)@cell_id, "integer64")
})

test_that("cell IDs given as doubles or strings name the same cells", {
  g <- hex_grid(resolution = 5, aperture = 3)
  id <- lonlat_to_cell(16.373, 48.2, g)
  expect_equal(cell_to_lonlat(as.numeric(id), g), cell_to_lonlat(id, g))
  expect_equal(cell_to_lonlat(as.character(id), g), cell_to_lonlat(id, g))
})

test_that("NA coordinates give NA cell IDs", {
  g <- hex_grid(resolution = 5, aperture = 7)
  id <- lonlat_to_cell(c(10, NA), c(NA, 20), g)
  expect_true(all(is.na(id)))
})

test_that("points round-trip at the finest resolution of each aperture", {
  p <- random_points(400)
  for (case in finest) {
    g <- hex_grid(resolution = case$resolution, aperture = case$aperture)
    id <- lonlat_to_cell(p$lon, p$lat, g)
    expect_false(anyNA(id))
    ctr <- cell_to_lonlat(id, g)

    # A cell's centre lies in the cell itself
    expect_identical(lonlat_to_cell(ctr$lon_deg, ctr$lat_deg, g), id)

    # and every point within about one hexagon circumradius of its cell's
    # centre, the hexagon of the grid's mean cell area
    d <- great_circle_km(p$lon, p$lat, ctr$lon_deg, ctr$lat_deg,
                         hexify:::grid_radius_km(g))
    circumradius <- sqrt(2 * g@area_km2 / (3 * sqrt(3)))
    expect_lt(max(d / circumradius), 1.25,
              label = sprintf("aperture %s, resolution %d", case$aperture,
                              case$resolution))
  }
})

test_that("the largest cell ID of each finest grid decodes", {
  # The last cell is the vertex quad's single cell, at the solid's vertex
  # opposite vertex 0
  vert_lon <- ISEA_VERT0_LON_DEG - 180
  vert_lat <- -ISEA_VERT0_LAT_DEG
  for (case in finest) {
    g <- hex_grid(resolution = case$resolution, aperture = case$aperture)
    n <- hexify:::isea_cell_count(g@aperture, g@resolution)
    ctr <- cell_to_lonlat(n, g)
    expect_equal(c(ctr$lon_deg, ctr$lat_deg), c(vert_lon, vert_lat),
                 tolerance = 1e-9)
    expect_identical(lonlat_to_cell(vert_lon, vert_lat, g), n)
    expect_error(cell_to_lonlat(n + 1L, g), "cell_id must be")
  }
})

test_that("neighbours are mutual at the finest resolution of each aperture", {
  p <- random_points(20)
  for (case in finest) {
    g <- hex_grid(resolution = case$resolution, aperture = case$aperture)
    id <- lonlat_to_cell(p$lon, p$lat, g)
    nb <- get_neighbors(id, g)
    for (k in seq_along(id)) {
      expect_length(nb[[k]], 6L)
      back <- get_neighbors(nb[[k]], g)
      expect_true(all(vapply(back, function(b) id[k] %in% b, logical(1))))
    }
  }
})

test_that("parent and children agree at the finest resolution of each aperture", {
  p <- random_points(20)
  for (case in finest) {
    g <- hex_grid(resolution = case$resolution, aperture = case$aperture)
    id <- lonlat_to_cell(p$lon, p$lat, g)
    parent <- get_parent(id, g)
    coarse <- hex_grid(resolution = case$resolution - 1L, aperture = case$aperture)
    kids <- get_children(parent, coarse)
    expect_true(all(vapply(seq_along(id), function(k) id[k] %in% kids[[k]],
                           logical(1))))
  }
})

test_that("grids whose cell IDs pass 2^63 - 1 are refused", {
  expect_s4_class(hex_grid(resolution = 30, aperture = 3), "HexGridInfo")
  expect_error(hex_grid(resolution = 30, aperture = 4),
               "finest resolution is 29")
  expect_error(hex_grid(resolution = 22, aperture = 7),
               "finest resolution is 21")
  expect_s4_class(hex_grid(resolution = 30, aperture = 4,
                           polyhedron = "octahedron"), "HexGridInfo")
  expect_error(hex_grid(resolution = 22, aperture = 7, polyhedron = "octahedron"),
               "finest resolution is 21")
  # 2 * 7^22 + 2 cells fit, 2 * 7^23 + 2 do not
  expect_s4_class(hex_grid(resolution = 22, aperture = 7,
                           polyhedron = "tetrahedron"), "HexGridInfo")
  expect_error(hex_grid(resolution = 23, aperture = 7, polyhedron = "tetrahedron"),
               "finest resolution is 22")
  expect_error(hex_grid(resolution = 25, aperture = "4/7"),
               "finest resolution is 24")
  expect_error(hex_grid(resolution = 22, aperture = rep(7, 22)),
               "than 64-bit cell IDs can number")
  g <- hex_grid(resolution = 21, aperture = 7)
  expect_error(get_children(lonlat_to_cell(0, 0, g), g), "finest resolution is 21")
})

test_that("a target area finer than the finest resolution clamps to it", {
  g <- hex_grid(area_km2 = 1e-20, aperture = 7)
  expect_identical(g@resolution, 21L)
})
