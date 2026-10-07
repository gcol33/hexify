test_that("a DGGRS definition carries the schema's required fields", {
  for (g in list(hex_grid(resolution = 4, aperture = 7),
                 hex_grid(resolution = 4, aperture = 3, projection = "fuller"),
                 hex_grid(resolution = 5, aperture = "4/3"),
                 hex_grid(resolution = 3, aperture = 7, polyhedron = "octahedron"),
                 hex_grid(resolution = 3, type = "h3"))) {
    d <- dggrs_definition(g)
    expect_true(all(c("dggh", "zirs", "subZoneOrder") %in% names(d)))
    def <- d$dggh$definition
    expect_equal(def$spatialDimensions, 2L)
    expect_equal(def$temporalDimensions, 0L)
    expect_true(is.character(d$zirs$textZIRS$description))
    expect_true(is.character(d$subZoneOrder$description))
  }
})

test_that("refinement ratio, strategy, zone types and equal size follow the grid", {
  d7 <- dggrs_definition(hex_grid(resolution = 4, aperture = 7))
  expect_equal(d7$title, "ISEA7H")
  expect_equal(d7$dggh$definition$refinementRatio, 7L)
  expect_equal(d7$dggh$definition$refinementStrategy,
               c("centredChildCell", "nodeSharingChildCell"))
  expect_true(d7$dggh$definition$constraints$cellEqualSized)
  expect_equal(d7$dggh$parameters$orientation$longitude, 11.25)

  d43 <- dggrs_definition(hex_grid(resolution = 5, aperture = "4/3"))
  expect_equal(d43$title, "ISEA43H")
  expect_equal(d43$dggh$definition$refinementRatio, c(4L, 4L, 3L, 3L, 3L))
  expect_equal(d43$dggh$definition$refinementStrategy,
               c("centredChildCell", "edgeCentredChildCell", "nodeCentredChildCell"))

  df <- dggrs_definition(hex_grid(resolution = 3, aperture = 3, projection = "fuller"))
  expect_equal(df$title, "FULLER3H")
  expect_false(df$dggh$definition$constraints$cellEqualSized)

  do <- dggrs_definition(hex_grid(resolution = 3, aperture = 7, polyhedron = "octahedron"))
  expect_equal(do$dggh$definition$zoneTypes, c("hexagon", "square"))
})

test_that("sub-zones come in ascending order of cell ID, as the definition states", {
  for (g in list(hex_grid(resolution = 3, aperture = 7),
                 hex_grid(resolution = 4, aperture = 3),
                 hex_grid(resolution = 3, aperture = 4),
                 hex_grid(resolution = 3, aperture = "4/3"),
                 hex_grid(resolution = 3, type = "h3"))) {
    pts <- sphere_test_points(40)
    kids <- get_children(lonlat_to_cell(pts$lon, pts$lat, g), g, levels = 2)
    for (k in kids) {
      key <- if (is_h3_grid(g)) as.numeric(paste0("0x", k)) else k
      expect_false(is.unsorted(key))
    }
  }
})
