test_that("hex_corners_to_sf builds a valid polygon", {
  skip_if_not_installed("sf")

  lon <- c(0, 1, 1, 0, -1, -1)
  lat <- c(0, 0.5, 1, 1, 0.5, 0)

  poly <- hex_corners_to_sf(lon, lat)

  # basic structure
  expect_s3_class(poly, "sf")
  expect_true(sf::st_is_valid(poly))
  expect_equal(nrow(poly), 1L)
  expect_identical(as.character(sf::st_geometry_type(poly)), "POLYGON")

  # coordinates layout: X, Y, L1, L2 (so at least 4 columns)
  coords <- sf::st_coordinates(poly)
  expect_true(all(c("X","Y") %in% colnames(coords)))

  # XY must match the provided points + closing vertex
  # XY must match the provided points + closing vertex
  xy_expected <- rbind(cbind(lon, lat), c(lon[1], lat[1]))

  actual_xy   <- unname(as.matrix(coords[, c("X","Y"), drop = FALSE]))
  expected_xy <- unname(as.matrix(xy_expected))

  expect_equal(actual_xy, expected_xy, tolerance = 0)
})
