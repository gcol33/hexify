# Grid orientation: where an ISEA grid's icosahedron sits on the sphere.
#
# A rotation of the sphere carries the standard icosahedron to any other
# orientation, and the grid with it. So a point rotated with the icosahedron
# lands in the same cell, a cell's centre and corners rotate with it, and the
# hierarchy and neighbours, which read cell IDs alone, do not change.

unit_rows <- function(lon, lat) {
  r <- pi / 180
  cbind(cos(lat * r) * cos(lon * r), cos(lat * r) * sin(lon * r), sin(lat * r))
}

# Angle between matching rows of two sets of unit vectors, in degrees
row_angle <- function(A, B) {
  cr <- cbind(A[, 2] * B[, 3] - A[, 3] * B[, 2], A[, 3] * B[, 1] - A[, 1] * B[, 3],
              A[, 1] * B[, 2] - A[, 2] * B[, 1])
  atan2(sqrt(rowSums(cr^2)), rowSums(A * B)) * 180 / pi
}

rows_lonlat <- function(P) {
  cbind(atan2(P[, 2], P[, 1]), asin(pmax(-1, pmin(1, P[, 3])))) * 180 / pi
}

# The rotation taking the standard icosahedron to orientation `o`, read from
# the two sets of vertices; rows map as P %*% t(R).
rotation_to <- function(o) {
  A <- cpp_icosa_solid(unname(ISEA_ORIENTATION))$vertices
  B <- cpp_icosa_solid(unname(o))$vertices
  s <- svd(t(A) %*% B)
  R <- s$v %*% t(s$u)
  expect_equal(A %*% t(R), B, tolerance = 1e-12)
  R
}

ORIENTS <- list(c(-40, 20, 33), c(0, 90, 0), c(150, -60, 200))

test_that("a grid takes the standard orientation unless told otherwise", {
  g <- hex_grid(resolution = 4)
  expect_identical(g@orientation, ISEA_ORIENTATION)
  expect_identical(summary(g)$orientation, ISEA_ORIENTATION)
  expect_false(any(grepl("Orientation", capture.output(print(g)))))
})

test_that("a numeric orientation is stored normalised", {
  g <- hex_grid(resolution = 4, orientation = c(190, 10, -30))
  expect_equal(unname(g@orientation), c(-170, 10, 330))
  expect_named(g@orientation, c("vert0_lon", "vert0_lat", "azimuth"))

  named <- hex_grid(resolution = 4,
                    orientation = c(azimuth = 5, vert0_lat = -20, vert0_lon = 40))
  expect_equal(unname(named@orientation), c(40, -20, 5))

  expect_true(any(grepl("Orientation", capture.output(print(g)))))
  expect_equal(as.list(g)$orientation, g@orientation)
})

test_that("bad orientations are refused", {
  expect_error(hex_grid(resolution = 4, orientation = c(0, 95, 0)), "vert0_lat")
  expect_error(hex_grid(resolution = 4, orientation = c(0, 10)), "three finite")
  expect_error(hex_grid(resolution = 4, orientation = c(0, NA, 0)), "three finite")
  expect_error(hex_grid(resolution = 4, orientation = c(a = 0, b = 1, c = 2)), "names")
  expect_error(hex_grid(resolution = 4, orientation = "sideways"), "orientation must")
  expect_error(hex_grid(resolution = 4, orientation = "region"), "needs a region")
  expect_error(hex_grid(resolution = 4, region = c(0, 0)), "region applies")
  expect_error(hex_grid(resolution = 4, type = "h3", orientation = c(0, 90, 0)),
               "H3 fixes its own orientation")
  expect_length(hex_grid(resolution = 4, type = "h3")@orientation, 0L)
})

test_that("points rotated with the icosahedron land in the same cells", {
  set.seed(80)
  lon <- runif(3000, -180, 180)
  lat <- asin(runif(3000, -1, 1)) * 180 / pi
  P <- unit_rows(lon, lat)

  for (o in ORIENTS) {
    R <- rotation_to(o)
    Q <- rows_lonlat(P %*% t(R))
    for (ap in list(3, 7, "4/3", c(4, 4, 7, 3, 4))) {
      res <- if (length(ap) > 1) length(ap) else 5
      g0 <- hex_grid(resolution = res, aperture = ap)
      g1 <- hex_grid(resolution = res, aperture = ap, orientation = o)

      ids <- lonlat_to_cell(lon, lat, g0)
      expect_identical(lonlat_to_cell(Q[, 1], Q[, 2], g1), ids)

      cells <- unique(ids)[1:50]
      c0 <- cell_to_lonlat(cells, g0)
      c1 <- cell_to_lonlat(cells, g1)
      expect_equal(unit_rows(c1[[1]], c1[[2]]),
                   unit_rows(c0[[1]], c0[[2]]) %*% t(R), tolerance = 1e-10)

      r0 <- isea_cell_rings(cells, g0@resolution, g0@aperture, orient_arg(g0), 0)
      r1 <- isea_cell_rings(cells, g1@resolution, g1@aperture, orient_arg(g1), 0)
      for (k in seq_along(cells)) {
        expect_equal(unit_rows(r1[[k]][, 1], r1[[k]][, 2]),
                     unit_rows(r0[[k]][, 1], r0[[k]][, 2]) %*% t(R),
                     tolerance = 1e-10)
      }
    }
  }
})

test_that("the hierarchy and neighbours do not depend on the orientation", {
  for (ap in list(3, 4, 7, "4/3")) {
    g0 <- hex_grid(resolution = 5, aperture = ap)
    g1 <- hex_grid(resolution = 5, aperture = ap, orientation = ORIENTS[[1]])
    cells <- c(1, 2, 7, 100, 333)
    expect_identical(get_parent(cells, g1), get_parent(cells, g0))
    expect_identical(get_children(cells[1:2], g1), get_children(cells[1:2], g0))
    expect_identical(grid_neighbors_isea(cells, g1), grid_neighbors_isea(cells, g0))
    if (!is_mixed_aperture(g0@aperture)) {
      expect_identical(cell_to_index(cells, g1), cell_to_index(cells, g0))
    }
  }
})

test_that("an oriented grid leaves no orientation behind", {
  lon <- c(0, 10, -120, 170, 45)
  lat <- c(0, 45, -30, 80, -70)
  g0 <- hex_grid(resolution = 6)
  g1 <- hex_grid(resolution = 6, orientation = ORIENTS[[2]])
  ids0 <- lonlat_to_cell(lon, lat, g0)
  qij0 <- hexify_lonlat_to_quad_ij(lon[1], lat[1], 6, 3)

  ids1 <- lonlat_to_cell(lon, lat, g1)
  expect_false(identical(ids1, ids0))
  expect_identical(lonlat_to_cell(lon, lat, g0), ids0)
  lonlat_to_cell(lon, lat, g1)
  expect_identical(hexify_lonlat_to_quad_ij(lon[1], lat[1], 6, 3), qij0)

  # The default orientation of the grid-less functions does not reach a grid
  hexify_build_icosa(0, 90, 0)
  on.exit(hexify_build_icosa(), add = TRUE)
  expect_identical(lonlat_to_cell(lon, lat, g0), ids0)
  expect_identical(lonlat_to_cell(lon, lat, g1), ids1)
})

test_that("the globe's face frames rotate with the icosahedron", {
  o <- ORIENTS[[3]]
  R <- rotation_to(o)
  f0 <- matrix(cpp_globe_projection(numeric(0))$faces, nrow = 16)
  f1 <- matrix(cpp_globe_projection(unname(o))$faces, nrow = 16)
  for (rows in list(1:3, 5:7, 9:11)) {
    expect_equal(t(f1[rows, ]), t(f0[rows, ]) %*% t(R), tolerance = 1e-12)
  }
  expect_identical(f1[13:16, ], f0[13:16, ])
})

test_that("a region orientation puts the centre on the middle of an edge", {
  phi <- (1 + sqrt(5)) / 2
  for (centre in list(c(10, 46.5), c(-75, -10), c(0, 90), c(170, -89))) {
    g <- hex_grid(resolution = 5, orientation = "region", region = centre)
    s <- icosa_solid(orient_arg(g))
    C <- drop(unit_rows(centre[1], centre[2]))
    # DGGRID's placement constants hold the edge midpoint to about 1e-7 degrees
    to_vertex <- sort(row_angle(s$vertices, matrix(C, 12, 3, byrow = TRUE)))
    expect_lt(max(abs(to_vertex[1:2] - atan(1 / phi) * 180 / pi)), 1e-7)
    to_face <- sort(row_angle(s$normals, matrix(C, 20, 3, byrow = TRUE)))
    expect_lt(abs(to_face[1] - to_face[2]), 1e-7)
  }
})

test_that("a face orientation puts the centre on the centre of a face", {
  for (centre in list(c(10, 46.5), c(-75, -10), c(0, 90), c(170, -89), c(-180, 0))) {
    g <- hex_grid(resolution = 5, orientation = "face", region = centre)
    s <- icosa_solid(orient_arg(g))
    C <- drop(unit_rows(centre[1], centre[2]))
    to_face <- sort(row_angle(s$normals, matrix(C, 20, 3, byrow = TRUE)))
    expect_lt(to_face[1], 1e-9)
  }
  expect_error(hex_grid(resolution = 4, orientation = "face"), "needs a region")
  expect_error(hex_grid(resolution = 4, orientation = "random", region = c(0, 0)),
               "region applies")
})

test_that("a region can be an sf object", {
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = 5, ymin = 44, xmax = 15, ymax = 49),
                                   crs = 4326))
  centre <- region_centre(box)
  expect_equal(hex_grid(resolution = 5, orientation = "region", region = box)@orientation,
               hex_grid(resolution = 5, orientation = "region", region = centre)@orientation)
  expect_equal(centre, region_centre(sf::st_bbox(box)))
  expect_error(region_centre(c(0, 100)), "region must")
})

test_that("a random orientation follows the seed", {
  set.seed(3)
  a <- hex_grid(resolution = 4, orientation = "random")@orientation
  set.seed(3)
  b <- hex_grid(resolution = 4, orientation = "random")@orientation
  expect_identical(a, b)
  expect_true(a[["vert0_lat"]] >= -90 && a[["vert0_lat"]] <= 90)
  expect_true(a[["azimuth"]] >= 0 && a[["azimuth"]] < 360)
})

test_that("a grid saved before grids carried an orientation reads the standard one", {
  g <- hex_grid(resolution = 5)
  attr(g, "orientation") <- NULL
  expect_false(.hasSlot(g, "orientation"))
  expect_identical(grid_orientation(g), ISEA_ORIENTATION)
  expect_identical(extract_grid(g)@orientation, ISEA_ORIENTATION)
  expect_identical(lonlat_to_cell(10, 45, g), lonlat_to_cell(10, 45, hex_grid(resolution = 5)))
})

test_that("cells and centres match DGGRID under other orientations", {
  # paper/bench/make_dggrid_orientation_fixture.R writes this from DGGRID
  fx <- read.csv(test_path("data", "dggrid_orientation.csv"),
                 colClasses = c(aperture = "character"))
  for (key in unique(paste(fx$placement, fx$aperture, fx$resolution))) {
    d <- fx[paste(fx$placement, fx$aperture, fx$resolution) == key, ]
    g <- if (d$placement[1] == "region") {
      hex_grid(resolution = d$resolution[1], aperture = d$aperture[1],
               orientation = "region", region = c(d$region_lon[1], d$region_lat[1]))
    } else {
      hex_grid(resolution = d$resolution[1], aperture = d$aperture[1],
               orientation = c(d$vert0_lon[1], d$vert0_lat[1], d$azimuth[1]))
    }
    expect_identical(lonlat_to_cell(d$lon, d$lat, g), as.numeric(d$seqnum), info = key)
    ctr <- cell_to_lonlat(d$seqnum, g)
    # DGGRID writes centres to 12 decimals
    gap <- row_angle(unit_rows(ctr[[1]], ctr[[2]]), unit_rows(d$centre_lon, d$centre_lat))
    expect_lt(max(gap), 1e-7)
  }
})
