# IGEO7's integer forms of the Z7 index: the packed 64-bit index (as
# integer64 and as hexadecimal) and the monotonic ID.

library(testthat)

# The signed 64-bit integer whose bits 16 hexadecimal digits spell
hex_to_integer64 <- function(hex) {
  from_hex <- function(h) {
    v <- bit64::as.integer64(0)
    for (ch in strsplit(h, "")[[1]]) v <- v * 16L + strtoi(ch, 16L)
    v
  }
  do.call(c, lapply(hex, function(h) {
    if (strtoi(substr(h, 1, 1), 16L) < 8L) return(from_hex(h))
    complement <- paste(sprintf("%x", 15L - strtoi(strsplit(h, "")[[1]], 16L)), collapse = "")
    -(from_hex(complement) + 1L)
  }))
}

test_that("the integer forms of the slides' Table 3 rows are IGEO7's", {
  # Kmoch et al. (FOSS4G Europe 2026), Table 3: the children of "0330" at level 4
  tab <- data.frame(
    z7 = c("033000", "033001", "033002", "033003", "033004", "033005", "033006",
           "033010", "033011", "033012", "033013", "033014", "033015", "033016",
           "033020", "033065", "033066"),
    monotonic = c(8232:8246, 8279, 8280),
    hex = c("0x3600ffffffffffff", "0x3601ffffffffffff", "0x3602ffffffffffff",
            "0x3603ffffffffffff", "0x3604ffffffffffff", "0x3605ffffffffffff",
            "0x3606ffffffffffff", "0x3608ffffffffffff", "0x3609ffffffffffff",
            "0x360affffffffffff", "0x360bffffffffffff", "0x360cffffffffffff",
            "0x360dffffffffffff", "0x360effffffffffff", "0x3610ffffffffffff",
            "0x3635ffffffffffff", "0x3636ffffffffffff"),
    uint64 = c("3891391553024819199", "3891673028001529855", "3891954502978240511",
               "3892235977954951167", "3892517452931661823", "3892798927908372479",
               "3893080402885083135", "3893643352838504447", "3893924827815215103",
               "3894206302791925759", "3894487777768636415", "3894769252745347071",
               "3895050727722057727", "3895332202698768383", "3895895152652189695",
               "3906309726790483967", "3906591201767194623"),
    stringsAsFactors = FALSE)
  g <- hex_grid(resolution = 4, aperture = 7)
  cell <- index_to_cell(tab$z7, g)
  expect_false(anyNA(cell))

  expect_identical(cell_to_index(cell, g, "hex"), sub("^0x", "", tab$hex))
  expect_identical(as.character(cell_to_index(cell, g, "int")), tab$uint64)
  expect_identical(cell_to_index(cell, g, "monotonic"),
                   bit64::as.integer64(tab$monotonic))

  expect_identical(index_to_cell(tab$hex, g, "hex"), cell)
  expect_identical(index_to_cell(toupper(sub("^0x", "", tab$hex)), g, "hex"), cell)
  expect_identical(index_to_cell(bit64::as.integer64(tab$uint64), g, "int"), cell)
  expect_identical(index_to_cell(tab$monotonic, g, "monotonic"), cell)
})

test_that("the packed index is DGGRID's INT64 output for every cell", {
  # paper/bench/make_dggrid_z7_int_fixture.R writes this from DGGRID: every cell
  # at resolutions 0 to 3 and the cells of random points at 5, 10, 15 and 20
  fx <- read.csv(test_path("data", "dggrid_z7_int.csv"),
                 colClasses = c(seqnum = "character", z7 = "character",
                                z7_hex = "character"))
  for (res in unique(fx$resolution)) {
    at <- fx[fx$resolution == res, ]
    g <- hex_grid(resolution = res, aperture = 7)
    ids <- bit64::as.integer64(at$seqnum)
    info <- sprintf("resolution %d", res)

    expect_identical(cell_to_index(ids, g), at$z7, info = info)
    expect_identical(cell_to_index(ids, g, "hex"), at$z7_hex, info = info)
    expect_identical(index_to_cell(at$z7_hex, g, "hex"), ids, info = info)
    expect_identical(index_to_cell(at$z7, g), ids, info = info)

    packed <- cell_to_index(ids, g, "int")
    expect_identical(packed, hex_to_integer64(at$z7_hex), info = info)
    expect_identical(index_to_cell(packed, g, "int"), ids, info = info)

    mono <- cell_to_index(ids, g, "monotonic")
    expect_identical(index_to_cell(mono, g, "monotonic"), ids, info = info)
  }
})

test_that("the packed index of base cells 08-11 reads as a negative integer64", {
  g <- hex_grid(resolution = 2, aperture = 7)
  z <- c("0000", "0766", "0800", "1166")
  packed <- cell_to_index(index_to_cell(z, g), g, "int")
  expect_identical(as.vector(packed < 0), c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(cell_to_index(index_to_cell(z, g), g, "hex"),
                   c("003fffffffffffff", "7dbfffffffffffff",
                     "803fffffffffffff", "bdbfffffffffffff"))
  expect_identical(packed, hex_to_integer64(c("003fffffffffffff", "7dbfffffffffffff",
                                              "803fffffffffffff", "bdbfffffffffffff")))
  # signed order puts base cells 08-11 first; the monotonic ID keeps Z7 order
  expect_identical(order(packed), c(3L, 4L, 1L, 2L))
  expect_identical(order(cell_to_index(index_to_cell(z, g), g, "monotonic")), 1:4)
})

test_that("base cell 08's pentagon at resolution 20 packs to the bits of NA_integer64_", {
  g <- hex_grid(resolution = 20, aperture = 7)
  pent <- index_to_cell(paste0("08", strrep("0", 20)), g)
  expect_identical(cell_to_index(pent, g, "hex"), "8000000000000000")
  expect_true(is.na(cell_to_index(pent, g, "int")))
  expect_identical(index_to_cell("8000000000000000", g, "hex"), pent)

  # every other index of resolution 20, the finest, is exact as integer64
  far <- index_to_cell(c(paste0("08", strrep("0", 19), "1"), paste0("11", strrep("6", 20))), g)
  packed <- cell_to_index(far, g, "int")
  expect_false(anyNA(packed))
  expect_identical(index_to_cell(packed, g, "int"), far)
  expect_identical(cell_to_index(far, g, "hex"), c("8000000000000001", "bdb6db6db6db6db6"))
})

test_that("monotonic IDs count every digit string, leaving the pentagons' deleted ones unused", {
  g <- hex_grid(resolution = 1, aperture = 7)
  mono <- as.integer(cell_to_index(as_cell_id(seq_len(72)), g, "monotonic"))
  # base * 7 + digit; digit 2 is deleted under base cells 00-05, 5 under 06-11
  deleted <- c(0:5 * 7 + 2, 6:11 * 7 + 5)
  expect_identical(sort(mono), setdiff(0:83, deleted))

  g3 <- hex_grid(resolution = 3, aperture = 7)
  ids <- as_cell_id(seq_len(n_cells(g3)))
  expect_identical(order(cell_to_index(ids, g3, "monotonic")),
                   order(cell_to_index(ids, g3)))
})

test_that("monotonic IDs reach resolution 21, the packed index resolution 20", {
  g21 <- hex_grid(resolution = 21, aperture = 7)
  cell <- lonlat_to_cell(16.37, 48.21, g21)
  mono <- cell_to_index(cell, g21, "monotonic")
  expect_identical(index_to_cell(mono, g21, "monotonic"), cell)
  expect_error(cell_to_index(cell, g21, "int"), "resolutions 0 to 20")
  expect_error(cell_to_index(cell, g21, "hex"), "resolutions 0 to 20")
})

test_that("integer forms reject grids and indices they do not cover", {
  g <- hex_grid(resolution = 3, aperture = 7)
  for (other in list(hex_grid(resolution = 3, aperture = 3),
                     hex_grid(resolution = 3, aperture = 7, polyhedron = "octahedron"),
                     hex_grid(resolution = 3, type = "h3"))) {
    expect_error(cell_to_index(1, other, "int"), "IGEO7")
    expect_error(index_to_cell("0", other, "hex"), "IGEO7")
  }
  # a digit after the first 7
  expect_error(index_to_cell("0e3fffffffffffff", g, "hex"), "not an IGEO7 packed index")
  expect_error(index_to_cell("3600ffffffffffff", g, "hex"),
               "resolution-4 cell, not one of this resolution-3 grid")
  expect_error(index_to_cell("c000ffffffffffff", g, "hex"), "not an IGEO7 packed index")
  expect_error(index_to_cell("xyz", g, "hex"), "hexadecimal")
  expect_error(index_to_cell(bit64::as.integer64(12 * 7^3), g, "monotonic"),
               "monotonic ID")
  cell <- index_to_cell("03300", g)
  expect_identical(index_to_cell(c(NA, cell_to_index(cell, g, "hex")), g, "hex"),
                   c(bit64::NA_integer64_, cell))
  expect_identical(cell_to_index(c(bit64::NA_integer64_, cell), g, "hex"),
                   c(NA, "3607ffffffffffff"))
})

test_that("index_to_cell() inverts cell_to_index() on every kind of grid", {
  set.seed(20261009)
  grids <- list(hex_grid(resolution = 5, aperture = 3),
                hex_grid(resolution = 4, aperture = 4),
                hex_grid(resolution = 4, aperture = 7),
                hex_grid(resolution = 3, aperture = 7, projection = "fuller"),
                hex_grid(resolution = 3, aperture = 7, polyhedron = "octahedron"),
                hex_grid(resolution = 3, aperture = "4/3"),
                hex_grid(resolution = 5, type = "h3"))
  for (g in grids) {
    lon <- runif(20, -180, 180)
    lat <- asin(runif(20, -1, 1)) * 180 / pi
    cells <- lonlat_to_cell(lon, lat, g)
    expect_identical(index_to_cell(cell_to_index(cells, g), g), cells,
                     info = paste(class(g), g@aperture))
  }
  g <- hex_grid(resolution = 3, aperture = 7)
  expect_error(index_to_cell("0100", g), "resolution-2 cell; the grid's cells are at resolution 3")
})
