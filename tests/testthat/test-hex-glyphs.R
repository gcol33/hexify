# Ray glyphs and neighbour triangles (#90)

test_that("rays point down at the lower limit, level in the middle and up at the upper", {
  g <- hex_grid(resolution = 4, aperture = 3)
  cell <- lonlat_to_cell(10, 30, g)
  cells <- rep(cell, 3)
  ctr <- cell_to_lonlat(cell, g)
  rays <- hex_rays(cell, 0, g, limits = c(-1, 1))
  tip <- sf::st_coordinates(rays)[2, 1:2]
  # A level ray leaves the centre due east, along a great circle
  la <- ctr$lat_deg * pi / 180
  la2 <- tip[[2]] * pi / 180
  dlon <- (tip[[1]] - ctr$lon_deg) * pi / 180
  bearing <- atan2(sin(dlon) * cos(la2), cos(la) * sin(la2) - sin(la) * cos(la2) * cos(dlon))
  expect_equal(bearing * 180 / pi, 90, tolerance = 1e-9)

  low <- sf::st_coordinates(hex_rays(cell, -1, g, limits = c(-1, 1)))[2, 1:2]
  high <- sf::st_coordinates(hex_rays(cell, 1, g, limits = c(-1, 1)))[2, 1:2]
  expect_lt(unname(low[2]), ctr$lat_deg)
  expect_gt(unname(high[2]), ctr$lat_deg)
  expect_equal(unname(low[1]), ctr$lon_deg, tolerance = 1e-9)
})

test_that("rays reach the given share of the way to a wall", {
  g <- hex_grid(resolution = 4, aperture = 3)
  cell <- lonlat_to_cell(-40, 10, g)
  rays <- hex_rays(cell, 0, g, limits = c(-1, 1), length = 1)
  xy <- sf::st_coordinates(rays)[, 1:2]
  d <- sf::st_distance(sf::st_sfc(sf::st_point(xy[1, ]), crs = 4326),
                       sf::st_sfc(sf::st_point(xy[2, ]), crs = 4326))
  r <- cell_inradius(cell, g) * grid_radius_km(g) * 1000
  expect_equal(as.numeric(d), r, tolerance = 5e-3)
})

test_that("a second variable points left and bounds give arcs", {
  g <- hex_grid(resolution = 4, aperture = 3)
  cells <- lonlat_to_cell(c(0, 20), c(10, 20), g)
  rays <- hex_rays(cells, c(0, 1), g, lower = c(-0.5, 0.5), upper = c(0.5, 1),
                   value2 = c(0.3, NA), limits = c(-1, 1))
  expect_equal(sort(table(paste(rays$side, rays$part))[c("left ray", "right arc", "right ray")]),
               sort(c("left ray" = 1L, "right arc" = 2L, "right ray" = 2L)),
               ignore_attr = TRUE)
  left <- sf::st_coordinates(rays[rays$side == "left", ])
  ctr <- cell_to_lonlat(cells[1], g)
  expect_lt(left[2, 1], ctr$lon_deg)
  expect_error(hex_rays(cells, c(0, 1), g, lower = c(0, 0)), "both bounds")
  expect_error(hex_rays(cells, 1, g), "one value per cell")
})

test_that("cells split into one triangle per wall, towards each neighbour", {
  g <- hex_grid(resolution = 3, aperture = 3)
  pent <- which(is_pentagon(seq_len(n_cells(g)), g))[1]
  hex <- lonlat_to_cell(10, 30, g)
  cells <- c(hex, pent)
  tri <- hex_triangles(cells, g)
  expect_equal(as.vector(table(tri$cell_id)[as.character(cells)]), c(6L, 5L))
  for (k in seq_along(cells)) {
    expect_setequal(tri$neighbor_id[tri$cell_id == cells[k]],
                    get_neighbors(cells[k], g)[[1]])
  }
  # The triangles of a cell cover it
  old <- sf::sf_use_s2(TRUE)
  on.exit(sf::sf_use_s2(old), add = TRUE)
  a_tri <- tapply(as.numeric(sf::st_area(tri)), tri$cell_id, sum)
  a_cell <- as.numeric(sf::st_area(cell_to_sf(cells, g)))
  expect_equal(as.vector(a_tri[as.character(cells)]), a_cell, tolerance = 1e-3)
})

test_that("triangles carry the change towards each neighbour", {
  g <- hex_grid(resolution = 3, aperture = 3)
  cell <- lonlat_to_cell(10, 30, g)
  nb <- get_neighbors(cell, g)[[1]]
  cells <- c(cell, nb[1:3])
  value <- c(10, 11, 13, NA)
  tri <- hex_triangles(cells, g, value = value)
  own <- tri[tri$cell_id == cell, ]
  expect_equal(own$change[match(nb[1:3], own$neighbor_id)], c(1, 3, NA))
  expect_true(all(is.na(own$change[!own$neighbor_id %in% cells])))
})
