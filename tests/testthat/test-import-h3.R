# tests/testthat/test-import-h3.R
# Tests for import_h3() (R/h3_compat.R). Previously had zero test coverage.

h3_ids <- c("8528342bfffffff", "85283473fffffff", "85283447fffffff")

test_that("import_h3 returns a HexGridInfo when data is NULL", {
  grid <- import_h3(h3_ids)
  expect_true(is_hex_grid(grid))
  expect_equal(grid@grid_type, "h3")
  expect_equal(grid@resolution, cpp_h3_getResolution(h3_ids)[1])
})

test_that("import_h3 returns a HexData when data is supplied", {
  df <- data.frame(species = c("oak", "pine", "birch"), count = c(10, 5, 3))
  hd <- import_h3(h3_ids, data = df)

  expect_true(is_hex_data(hd))
  expect_equal(nrow(as.data.frame(hd)), 3)
  expect_equal(as.character(hd@cell_id), h3_ids)
})

test_that("import_h3 validates cell_ids type", {
  expect_error(import_h3(c(1, 2, 3)), "character vector")
})

test_that("import_h3 errors on empty input", {
  expect_error(import_h3(character(0)), "must not be empty")
})

test_that("import_h3 errors on all-NA input", {
  expect_error(import_h3(c(NA_character_, NA_character_)), "only NA")
})

test_that("import_h3 drops NA entries but keeps valid ones for resolution inference", {
  ids_with_na <- c(h3_ids[1], NA, h3_ids[2])
  grid <- import_h3(ids_with_na)
  expect_true(is_hex_grid(grid))
})

test_that("import_h3 rejects invalid H3 cell IDs when validate = TRUE", {
  expect_error(import_h3(c(h3_ids[1], "not-a-real-cell")), "invalid H3 cell")
})

test_that("import_h3 skips validation when validate = FALSE", {
  # Even a malformed ID doesn't error at import time when validation is off;
  # only resolution inference runs.
  expect_no_error(import_h3(h3_ids, validate = FALSE))
})

test_that("import_h3 errors on mixed H3 resolutions", {
  # h3_ids[1:2] are resolution 5; get a resolution 6 child of a different cell
  child <- get_children(h3_ids[1], hex_grid(resolution = 5, type = "h3"))[[1]][1]
  expect_error(import_h3(c(h3_ids[1], child)), "same H3 resolution")
})

test_that("import_h3 validates data dimensions", {
  bad_df <- data.frame(x = 1:2)  # only 2 rows for 3 cell_ids
  expect_error(import_h3(h3_ids, data = bad_df), "must match")
})

test_that("import_h3 requires data to be a data.frame", {
  expect_error(import_h3(h3_ids, data = list(x = 1:3)), "data.frame")
})
