# OGC's zone identifiers (textZIRS, uint64ZIRS) of ISEA3H, ISEA7H, IVEA3H
# and IVEA7H against DGGAL's. The fixture holds a sample of DGGAL's zones
# (paper/bench/make_dggal_zirs_fixture.R); every zone up to ISEA3H level 10
# and ISEA7H level 5 is compared in paper/bench/bench_dggal_zirs.R.

ogc_grid <- function(dggrs, level) {
  hex_grid(resolution = level, aperture = as.integer(substr(dggrs, 5, 5)),
           projection = tolower(substr(dggrs, 1, 4)), ellipsoid = "WGS84",
           orientation = "ogc")
}

test_that("cells carry DGGAL's identifiers and centroids on OGC's grids", {
  fx <- read.csv(test_path("data", "dggal_zirs.csv"), colClasses = c(uint64 = "character"))
  for (key in unique(paste(fx$dggrs, fx$level))) {
    z <- fx[paste(fx$dggrs, fx$level) == key, ]
    g <- ogc_grid(z$dggrs[1], z$level[1])
    cell <- lonlat_to_cell(z$lon, z$lat, g)
    expect_false(anyDuplicated(cell) > 0, label = key)
    # DGGAL's centroids are geodetic on WGS84: the conversion agrees
    ctr <- cell_to_lonlat(cell, g)
    dlon <- (ctr$lon_deg - z$lon + 180) %% 360 - 180
    expect_lt(max(abs(dlon) * cos(z$lat * pi / 180), abs(ctr$lat_deg - z$lat)), 1e-10,
              label = key)
    expect_identical(cell_to_index(cell, g, form = "textZIRS"), z$text, label = key)
    expect_identical(as.character(cell_to_index(cell, g, form = "uint64ZIRS")), z$uint64,
                     label = key)
    expect_identical(index_to_cell(z$text, g, form = "textZIRS"), cell, label = key)
    expect_identical(index_to_cell(z$uint64, g, form = "uint64ZIRS"), cell, label = key)
    expect_identical(index_to_cell(bit64::as.integer64(z$uint64), g, form = "uint64ZIRS"),
                     cell, label = key)
  }
})

test_that("ISEA3H zone E6-317-A lies on the Crimean peninsula", {
  g <- ogc_grid("ISEA3H", 8L)
  cell <- index_to_cell("E6-317-A", g, form = "textZIRS")
  expect_identical(lonlat_to_cell(34.78, 45.43, g), cell)
  expect_identical(cell_to_index(cell, g, form = "textZIRS"), "E6-317-A")
  expect_identical(as.character(cell_to_index(cell, g, form = "uint64ZIRS")),
                   "630503947831872604")
})

test_that("every cell of a grid has its own identifier, read back to it", {
  for (spec in list(list("ISEA3H", 5L), list("ISEA7H", 3L), list("IVEA7H", 2L))) {
    g <- ogc_grid(spec[[1]], spec[[2]])
    ids <- grid_cells(g)
    text <- cell_to_index(ids, g, form = "textZIRS")
    expect_false(anyDuplicated(text) > 0)
    expect_identical(index_to_cell(text, g, form = "textZIRS"), ids)
    u <- cell_to_index(ids, g, form = "uint64ZIRS")
    expect_identical(index_to_cell(u, g, form = "uint64ZIRS"), ids)
  }
  # the identifiers name positions, not places: any orientation has them
  g <- hex_grid(resolution = 4, aperture = 7)
  ids <- grid_cells(g)
  expect_identical(cell_to_index(ids, g, form = "textZIRS"),
                   cell_to_index(ids, ogc_grid("ISEA7H", 4L), form = "textZIRS"))
})

test_that("strings that name no zone read as NA, other levels as an error", {
  g <- ogc_grid("ISEA3H", 3L)
  bad <- c("B0-0-E", "BC-0-B", "B0-9-B", "BA-1-B", "nonsense", NA)
  expect_identical(index_to_cell(bad, g, form = "textZIRS"), as_cell_id(rep(NA, 6)))
  expect_error(index_to_cell("B0-0-A", g, form = "textZIRS"), "resolution-2")
  expect_error(index_to_cell("C0-0-A", g, form = "textZIRS"), "resolution-4")
  g7 <- ogc_grid("ISEA7H", 3L)
  expect_true(is.na(index_to_cell("DA-0-H", g7, form = "textZIRS")))
  expect_true(is.na(index_to_cell("D0-0-H", g7, form = "textZIRS")))
  expect_false(is.na(index_to_cell("D0-1-H", g7, form = "textZIRS")))
})

test_that("only aperture-3 and aperture-7 grids on the icosahedron have them", {
  expect_error(cell_to_index(1, hex_grid(resolution = 3, aperture = 4), form = "textZIRS"),
               "aperture-3 and aperture-7")
  expect_error(cell_to_index(1, hex_grid(resolution = 3, aperture = 3, polyhedron = "octahedron"),
                             form = "textZIRS"), "icosahedron")
  expect_error(cell_to_index("8001fffffffffff", hex_grid(resolution = 0, type = "h3"),
                             form = "uint64ZIRS"), "aperture-3 and aperture-7")
  expect_error(cell_to_index(1, hex_grid(resolution = 20, aperture = 7), form = "textZIRS"),
               "up to level 19")
})

test_that("index_to_cell() reads hexify's own index strings back", {
  for (g in list(hex_grid(resolution = 5, aperture = 3), hex_grid(resolution = 4, aperture = 4),
                 hex_grid(resolution = 3, aperture = 7), hex_grid(resolution = 3, aperture = "4/3"))) {
    ids <- grid_cells(g)[seq(1, n_cells(g), by = 7)]
    expect_identical(index_to_cell(cell_to_index(ids, g), g), ids)
  }
  h3 <- hex_grid(resolution = 3, type = "h3")
  cells <- lonlat_to_cell(c(0, 10), c(45, 50), h3)
  expect_identical(index_to_cell(cell_to_index(cells, h3), h3), cells)
  g <- hex_grid(resolution = 3, aperture = 7)
  expect_error(index_to_cell(cell_to_index(5, hex_grid(resolution = 2, aperture = 7)), g),
               "resolution-2")
})
