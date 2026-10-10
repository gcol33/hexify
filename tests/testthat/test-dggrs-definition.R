test_that("a DGGRS definition carries the schema's required fields", {
  for (g in list(hex_grid(resolution = 4, aperture = 7),
                 hex_grid(resolution = 4, aperture = 3, projection = "fuller"),
                 hex_grid(resolution = 5, aperture = "4/3"),
                 hex_grid(resolution = 3, aperture = 7, polyhedron = "octahedron"),
                 hex_grid(resolution = 3, aperture = 3, polyhedron = "tetrahedron"),
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

test_that("a grid on an ellipsoid states its earth model, and OGC's grids OGC's identifiers", {
  g <- hex_grid(resolution = 4, aperture = 3, ellipsoid = "WGS84", orientation = "ogc")
  d <- dggrs_definition(g)
  expect_identical(d$dggh$parameters$ellipsoid, "[EPSG:7030]")
  # OGC's ISEA3H states the vertex's geodetic latitude as 58.397145907431
  expect_equal(d$dggh$parameters$orientation$latitude, 58.397145907431, tolerance = 1e-12)
  expect_equal(d$dggh$parameters$orientation$longitude, 11.20)
  expect_identical(d$zirs$textZIRS$type, "levelRootFaceHexRowMajorSubZone")
  expect_identical(d$links[[1]]$href, "https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H")
  expect_null(d$uri)
  expect_match(d$description, "WGS84 ellipsoid")
  expect_identical(dggrs_definition(hex_grid(resolution = 3, aperture = 7, projection = "ivea",
                                             ellipsoid = "WGS84", orientation = "ogc"))$links[[1]]$href,
                   "https://www.opengis.net/def/dggrs/OGC/1.0/IVEA7H")
  # the standard orientation, or no ellipsoid, is not OGC's grid
  expect_null(dggrs_definition(hex_grid(resolution = 4, aperture = 3, ellipsoid = "WGS84"))$links)
  expect_null(dggrs_definition(hex_grid(resolution = 4, aperture = 3, orientation = "ogc"))$links)
  dm <- dggrs_definition(hex_grid(resolution = 3, aperture = 7, ellipsoid = "mars"))
  expect_equal(dm$dggh$parameters$ellipsoid$semiMajorAxis_km, 3396.19)
  expect_identical(dm$zirs$textZIRS$type, "hierarchicalConcatenation")
})
