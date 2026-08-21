# tests/testthat/test-cell-geometry.R
# Geometry of a cell boundary: orientation, pentagons, and the two shapes the
# same boundary is returned in.

gc_distance <- function(a, b) {
  r <- pi / 180
  lon1 <- a[, 1] * r; lat1 <- a[, 2] * r
  lon2 <- b[, 1] * r; lat2 <- b[, 2] * r
  acos(pmin(1, pmax(-1, sin(lat1) * sin(lat2) +
                       cos(lat1) * cos(lat2) * cos(lon2 - lon1))))
}

lonlat_matrix <- function(df) {
  as.matrix(df[, c("lon_deg", "lat_deg")])
}

test_that("the corner list and the polygon table are the same boundary", {
  for (ap in c(3, 4, 7)) {
    cells <- as.numeric(c(1, 5, 40, 10 * ap^3 + 2))

    rings <- cpp_cell_to_corners(cells, 3L, ap)
    table <- cpp_cell_to_polygon(cells, 3L, ap)

    for (k in seq_along(cells)) {
      rows <- table[table$hex_id == cells[k], ]
      expect_equal(nrow(rows), nrow(rings[[k]]))
      expect_equal(rows$lon, unname(rings[[k]][, "lon"]))
      expect_equal(rows$lat, unname(rings[[k]][, "lat"]))
      expect_equal(rows$order, seq_len(nrow(rows)))
    }
  }
})

test_that("corners fall between the neighbouring cell centres", {
  # A corner is one circumradius from its own centre and the same from the two
  # neighbours sharing it. Turning the lattice by 30 degrees would aim the
  # corners at the neighbouring centres instead, which sit sqrt(3) - 1 of a
  # circumradius away.
  for (ap in c(3, 4, 7)) {
    for (res in 2:4) {
      g <- hex_grid(resolution = res, aperture = ap)
      cells <- round(seq(3, 10 * ap^res, length.out = 12))
      cells <- cells[!is_pentagon(cells, g)]

      centres <- lonlat_matrix(cell_to_lonlat(cells, g))
      rings <- cpp_cell_to_corners(as.numeric(cells), as.integer(res),
                                   as.integer(ap))

      ratios <- numeric(0)
      for (k in seq_along(cells)) {
        neighbours <- get_neighbors(cells[k], g)[[1]]
        neighbours <- neighbours[!is.na(neighbours)]
        nb_centres <- lonlat_matrix(cell_to_lonlat(neighbours, g))

        corners <- rings[[k]]
        corners <- corners[seq_len(nrow(corners) - 1), , drop = FALSE]

        for (v in seq_len(nrow(corners))) {
          one <- matrix(corners[v, ], nrow = 1)
          own <- gc_distance(one, matrix(centres[k, ], nrow = 1))
          near <- min(gc_distance(
            matrix(corners[v, ], nrow = nrow(nb_centres), ncol = 2, byrow = TRUE),
            nb_centres
          ))
          ratios <- c(ratios, near / own)
        }
      }

      expect_gt(min(ratios), 0.78)
    }
  }
})

test_that("a cell at an icosahedral vertex is a pentagon", {
  for (ap in c(3, 4, 7)) {
    for (res in 2:3) {
      g <- hex_grid(resolution = res, aperture = ap)
      n <- 10 * ap^res + 2
      pentagons <- seq_len(n)[is_pentagon(seq_len(n), g)]
      expect_length(pentagons, 12)

      rings <- cpp_cell_to_corners(as.numeric(pentagons), as.integer(res),
                                   as.integer(ap))
      for (ring in rings) {
        expect_equal(nrow(ring), 6L)
        expect_equal(ring[1, ], ring[6, ])
      }
    }
  }
})

test_that("every pentagon corner is shared with a neighbour", {
  # Five quads meet at an icosahedral vertex, so one of the six sectors around
  # the cell there is the icosahedron's angular deficit. Dropping the wrong
  # corner leaves one that no neighbour reaches.
  for (ap in c(3, 4, 7)) {
    for (res in 2:3) {
      g <- hex_grid(resolution = res, aperture = ap)
      n <- 10 * ap^res + 2
      pentagons <- seq_len(n)[is_pentagon(seq_len(n), g)]

      centres <- lonlat_matrix(cell_to_lonlat(pentagons, g))
      rings <- cpp_cell_to_corners(as.numeric(pentagons), as.integer(res),
                                   as.integer(ap))

      for (k in seq_along(pentagons)) {
        neighbours <- get_neighbors(pentagons[k], g)[[1]]
        neighbours <- neighbours[!is.na(neighbours)]
        expect_length(neighbours, 5L)

        nb_centres <- lonlat_matrix(cell_to_lonlat(neighbours, g))
        corners <- rings[[k]][1:5, , drop = FALSE]

        ratios <- vapply(1:5, function(v) {
          own <- gc_distance(matrix(corners[v, ], nrow = 1),
                             matrix(centres[k, ], nrow = 1))
          near <- min(gc_distance(
            matrix(corners[v, ], nrow = nrow(nb_centres), ncol = 2, byrow = TRUE),
            nb_centres
          ))
          near / own
        }, numeric(1))

        expect_lt(max(ratios) / min(ratios) - 1, 0.01)
      }
    }
  }
})

test_that("cell_to_sf follows the aperture of each level", {
  skip_if_not_installed("sf")

  lon <- c(0, 30, -100, 140, -20)
  lat <- c(45, -10, 60, -35, 5)

  for (ap in c("3", "4", "7", "4/3", "4/7", "7/4")) {
    g <- hex_grid(resolution = 4, aperture = ap)
    cells <- lonlat_to_cell(lon = lon, lat = lat, grid = g)
    centres <- cell_to_lonlat(cells, g)

    polygons <- cell_to_sf(cells, g, wrap_dateline = FALSE)
    points <- sf::st_as_sf(centres, coords = c("lon_deg", "lat_deg"), crs = 4326)

    hit <- suppressMessages(
      sf::st_intersects(points, sf::st_geometry(polygons), sparse = FALSE)
    )
    expect_true(all(diag(hit)))
  }
})

test_that("a mixed sequence and its first aperture give different polygons", {
  skip_if_not_installed("sf")

  cells <- 100:104

  mixed <- hex_grid(resolution = 4, aperture = "4/3")
  pure <- hex_grid(resolution = 4, aperture = 4)

  mixed_xy <- sf::st_coordinates(cell_to_sf(cells, mixed, wrap_dateline = FALSE))
  pure_xy <- sf::st_coordinates(cell_to_sf(cells, pure, wrap_dateline = FALSE))

  expect_false(isTRUE(all.equal(mixed_xy, pure_xy)))
})
