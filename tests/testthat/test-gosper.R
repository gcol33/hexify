skip_if_not_installed("sf")

every_cell <- function(g) seq_len(10 * as.numeric(g@aperture)^g@resolution + 2)

# Each fine-cell centre lies in exactly one shape, and that of the expected cell
expect_tiling <- function(shapes, fine, expected_cell) {
  centre <- cell_to_lonlat(every_cell(fine), fine)
  pts <- sf::st_as_sf(centre, coords = c("lon_deg", "lat_deg"), crs = 4326)
  hits <- sf::st_intersects(pts, shapes)
  expect_true(all(lengths(hits) == 1L))
  if (!is.null(expected_cell)) {
    expect_equal(shapes$cell_id[unlist(hits)], expected_cell)
  }
}

test_that("a Gosper island keeps its cell's area", {
  for (ap in c(7, 4, 3)) {
    g <- hex_grid(resolution = 1, aperture = ap)
    ids <- every_cell(g)
    islands <- cell_to_sf(ids, g, shape = "gosper", depth = 3, densify = 1e-4,
                          wrap_dateline = FALSE)

    expect_equal(islands$cell_id, ids)
    expect_true(all(sf::st_is_valid(islands)))
    area <- as.numeric(sf::st_area(islands)) / 1e6
    expect_lt(max(abs(area / cell_area(ids, g) - 1)), 1e-4)
  }
})

test_that("the Gosper islands of a grid tile the sphere", {
  g <- hex_grid(resolution = 1, aperture = 7)
  islands <- cell_to_sf(every_cell(g), g, shape = "gosper", depth = 2,
                        wrap_dateline = FALSE)
  expect_tiling(islands, hex_grid(resolution = 3, aperture = 7), NULL)
})

test_that("a Gosper edge turns the same way at every step", {
  g <- hex_grid(resolution = 2, aperture = 7)
  cell <- lonlat_to_cell(10, 30, g)
  ring <- function(depth) {
    sf::st_coordinates(cell_to_sf(cell, g, shape = "gosper", depth = depth,
                                  densify = 0))[, 1:2]
  }
  # 3^depth segments per edge, six edges, ring closed
  expect_equal(nrow(ring(1)), 6 * 3 + 1)
  expect_equal(nrow(ring(3)), 6 * 27 + 1)

  # The island's corners stay the hexagon's, every 3^depth-th vertex
  hex <- sf::st_coordinates(cell_to_sf(cell, g, densify = 0))[, 1:2]
  r3 <- ring(3)
  expect_equal(unname(r3[seq(1, nrow(r3), by = 27), ]), unname(hex), tolerance = 1e-7)
})

test_that("Gosper islands follow the grid's orientation and projection", {
  for (g in list(hex_grid(resolution = 1, aperture = 7, projection = "fuller"),
                 hex_grid(resolution = 1, aperture = 7, orientation = c(10, 50, 30)))) {
    islands <- cell_to_sf(every_cell(g), g, shape = "gosper", depth = 2,
                          wrap_dateline = FALSE)
    expect_true(all(sf::st_is_valid(islands)))
    fine <- grid_at_resolution(g, 3)
    expect_tiling(islands, fine, NULL)
  }
})

test_that("a descendant outline covers its cell's area", {
  g <- hex_grid(resolution = 1, aperture = 7)
  ids <- every_cell(g)
  outlines <- cell_to_sf(ids, g, shape = "descendants", depth = 2,
                         densify = 1e-4, wrap_dateline = FALSE)

  expect_equal(outlines$cell_id, ids)
  expect_true(all(sf::st_is_valid(outlines)))
  area <- as.numeric(sf::st_area(outlines)) / 1e6
  expect_lt(max(abs(area / cell_area(ids, g) - 1)), 1e-4)
})

test_that("descendant outlines tile the sphere along the hierarchy", {
  g <- hex_grid(resolution = 1, aperture = 7)
  outlines <- cell_to_sf(every_cell(g), g, shape = "descendants", depth = 2,
                         wrap_dateline = FALSE)
  fine <- hex_grid(resolution = 3, aperture = 7)
  expect_tiling(outlines, fine, get_parent(every_cell(fine), fine, levels = 2))
})

test_that("a descendant outline is the union of its children's outlines", {
  g <- hex_grid(resolution = 2, aperture = 7)
  parent <- lonlat_to_cell(16, 48, g)
  outline <- cell_to_sf(parent, g, shape = "descendants", depth = 2)

  child_grid <- hex_grid(resolution = 3, aperture = 7)
  kids <- get_children(parent, g)[[1]]
  child_outlines <- cell_to_sf(kids, child_grid, shape = "descendants",
                               depth = 1)
  expect_length(kids, 7L)

  # Centres of the resolution-5 cells around the parent fall in the outline
  # exactly when they fall in a child's outline
  fine <- hex_grid(resolution = 5, aperture = 7)
  near <- unique(unlist(get_children(c(parent, unlist(get_neighbors(parent, g))),
                                     g, levels = 3)))
  centre <- cell_to_lonlat(near, fine)
  pts <- sf::st_as_sf(centre, coords = c("lon_deg", "lat_deg"), crs = 4326)
  in_parent <- lengths(sf::st_intersects(pts, outline)) > 0
  in_child <- lengths(sf::st_intersects(pts, child_outlines)) > 0
  expect_gt(sum(in_parent), 0)
  expect_equal(in_parent, in_child)
})

test_that("shapes crossing the antimeridian are split there", {
  g <- hex_grid(resolution = 2, aperture = 7)
  cell <- lonlat_to_cell(180, -16, g)
  for (shape in c("gosper", "descendants")) {
    split <- cell_to_sf(cell, g, shape = shape, depth = 2)
    lons <- sf::st_coordinates(split)[, 1]
    expect_true(all(lons >= -180 & lons <= 180))
    expect_true(any(lons < -170) && any(lons > 170))
  }
})

test_that("cell shapes reject grids they are not defined on", {
  h3 <- hex_grid(resolution = 2, type = "h3")
  expect_error(cell_to_sf("821c07fffffffff", h3, shape = "gosper"), "ISEA")
  expect_error(cell_to_sf(1, hex_grid(resolution = 2, aperture = 3),
                          shape = "descendants"), "aperture 7")
  expect_error(cell_to_sf(1, hex_grid(resolution = 2, aperture = "4/7"),
                          shape = "descendants", depth = 1), "aperture 7")
  g <- hex_grid(resolution = 2, aperture = 7)
  expect_error(cell_to_sf(1, g, shape = "gosper", depth = 0), "positive whole")
  expect_error(cell_to_sf(1, g, shape = "descendants", depth = 29),
               "maximum resolution")
})

test_that("plots and sf exports draw cell shapes", {
  df <- data.frame(lon = c(16.4, 2.35), lat = c(48.2, 48.9))
  x <- hexify(df, lon = "lon", lat = "lat", resolution = 3, aperture = 7)

  polys <- sf::st_as_sf(x, geometry = "polygon", shape = "gosper", depth = 2)
  expect_gt(nrow(sf::st_coordinates(polys)), 2 * 6 * 9)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(x, basemap = FALSE, shape = "gosper", depth = 2))
  expect_invisible(plot(x, basemap = FALSE, shape = "descendants", depth = 2))
})
