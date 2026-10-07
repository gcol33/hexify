# Tests for the aperture-7 Z7 index, IGEO7's on the icosahedron.

library(testthat)

test_that("Z7 strings are DGGRID's for every cell", {
  # paper/bench/make_dggrid_fixture.R writes this from DGGRID: every cell at
  # resolutions 0 to 3 and a sample at resolution 7
  fx <- read.csv(test_path("data", "dggrid_z7.csv"),
                 colClasses = c(seqnum = "numeric", z7 = "character"))
  for (res in unique(fx$resolution)) {
    at <- fx[fx$resolution == res, ]
    g <- hex_grid(resolution = res, aperture = 7)
    expect_identical(cell_to_index(at$seqnum, g), at$z7, info = sprintf("resolution %d", res))
    expect_equal(hexify:::isea_index_to_cells(at$z7, 7L, "z7", hexify:::icosa_arg(g)),
                 at$seqnum, info = sprintf("resolution %d", res))
  }
})

test_that("Z7: resolution 0 is the twelve base cells", {
  for (face in 0:11) {
    idx <- hexify_cell_to_index(face, 0L, 0L, 0L, 7L, "z7")
    expect_equal(idx, sprintf("%02d", face))
    result <- hexify_index_to_cell(idx, 7L, "z7")
    expect_equal(as.integer(result$face), face)
    expect_equal(c(result$i, result$j), c(0, 0))
    expect_equal(as.integer(result$resolution), 0L)
  }
})

test_that("Z7: a pentagon has six children, missing the deleted direction", {
  kids <- hexify_get_children(c("01", "11", "0100", "0110"), 7L, "z7")
  expect_identical(kids[[1]], paste0("01", c(0, 1, 3, 4, 5, 6)))
  expect_identical(kids[[2]], paste0("11", c(0, 1, 2, 3, 4, 6)))
  expect_identical(kids[[3]], paste0("0100", c(0, 1, 3, 4, 5, 6)))
  expect_identical(kids[[4]], paste0("0110", 0:6))

  g <- hex_grid(resolution = 1, aperture = 7)
  expect_length(get_children(hexify:::isea_index_to_cells("01", 7L, "z7",
                                                          hexify:::icosa_arg(g)),
                             hexify:::grid_at_resolution(g, 0))[[1]], 6L)
})

test_that("Z7: a string in the deleted subsequence canonicalises to its cell", {
  for (idx in c("012", "0102", "01026", "115", "11005")) {
    canon <- hexify_z7_canonical(idx)
    expect_false(identical(canon, idx), info = idx)
    expect_equal(hexify_z7_canonical(canon), canon, info = idx)
    expect_equal(nchar(canon), nchar(idx), info = idx)
  }
  for (idx in c("0103", "110001", "0944444", "1066666")) {
    expect_equal(hexify_z7_canonical(idx), idx)
  }
})

test_that("Z7: malformed indices are rejected", {
  expect_error(hexify_index_to_cell("1", 7L, "z7"))
  expect_error(hexify_index_to_cell("", 7L, "z7"))
  expect_error(hexify_index_to_cell("12", 7L, "z7"), "base cell")
  expect_error(hexify_index_to_cell("990", 7L, "z7"), "base cell")
})

test_that("Z7: the index carries the two-digit base cell and one digit per level", {
  for (res in 1:3) {
    cell <- hexify_lonlat_to_cell(10, 45, res, 7L)
    qij <- hexify_cell_to_quad_ij(cell, res, 7L)
    idx <- hexify_cell_to_index(qij$quad, qij$i, qij$j, res, 7L, "z7")
    expect_equal(nchar(idx), 2L + res)
  }
})

test_that("Z7: a coordinate outside its quad is rejected", {
  # (49, 49) at resolution 3 lies outside quad 5's substrate box, so no cell of
  # that quad carries it and the hierarchy walk leaves the quad's own base cell.
  expect_error(hexify_cell_to_index(5L, 49L, 49L, 3L, 7L, "z7"),
               "does not lie in the given quad")
})

test_that("Z7: every cell round-trips through its index (#53)", {
  set.seed(53)
  for (res in 1:5) {
    n_cells <- 10 * 7^res + 2
    ids <- if (n_cells <= 600) seq_len(n_cells) else sample.int(n_cells, 600)

    qij <- hexify_cell_to_quad_ij(ids, res, 7L)
    idx <- vapply(seq_along(ids), function(k) {
      hexify_cell_to_index(qij$quad[k], qij$i[k], qij$j[k], res, 7L, "z7")
    }, character(1))

    back <- cell_id_unlist(lapply(idx, function(s) {
      r <- hexify_index_to_cell(s, 7L, "z7")
      hexify_quad_ij_to_cell(r$face, r$i, r$j, r$resolution, 7L)
    }))

    expect_equal(back, as_cell_id(ids),
                 info = sprintf("resolution %d round-trip", res))
    expect_equal(length(unique(idx)), length(idx),
                 info = sprintf("resolution %d indices are distinct", res))
    expect_true(all(nchar(idx) == 2L + res),
                info = sprintf("resolution %d index length", res))
  }
})

test_that("Z7: the two cells #53 reported as sharing a string get DGGRID's distinct strings", {
  # DGGRID writes "0061540" and "0045310" for these at resolution 5
  cells <- hexify_lonlat_to_cell(c(5, -34.9), c(45, 60.2), 5, 7L)
  qij <- hexify_cell_to_quad_ij(cells, 5, 7L)
  idx <- vapply(1:2, function(k) {
    hexify_cell_to_index(qij$quad[k], qij$i[k], qij$j[k], 5, 7L, "z7")
  }, character(1))
  expect_identical(idx, c("0061540", "0045310"))
})
