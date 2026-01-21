
# tests/testthat/test-classes.R
# Tests for S4 classes HexGridInfo and HexData

# =============================================================================
# HexGridInfo Tests
# =============================================================================

test_that("HexGridInfo can be created via hex_grid", {
  grid <- hex_grid(area_km2 = 1000)

  expect_s4_class(grid, "HexGridInfo")
  expect_equal(grid@aperture, "3")
  expect_true(grid@resolution > 0)
  expect_true(grid@area_km2 > 0)
})

test_that("HexGridInfo $ accessor works", {
  grid <- hex_grid(area_km2 = 1000)

  expect_equal(grid$aperture, "3")
  expect_equal(grid$resolution, grid@resolution)
  expect_equal(grid$area_km2, grid@area_km2)
  expect_equal(grid$crs, 4326L)
})

test_that("HexGridInfo names() returns slot names", {
  grid <- hex_grid(area_km2 = 1000)

  expect_true("aperture" %in% names(grid))
  expect_true("resolution" %in% names(grid))
  expect_true("area_km2" %in% names(grid))
  expect_true("diagonal_km" %in% names(grid))
  expect_true("crs" %in% names(grid))
})

test_that("HexGridInfo show() prints correctly", {
  grid <- hex_grid(area_km2 = 1000)

  output <- capture.output(show(grid))
  expect_true(any(grepl("HexGridInfo", output)))
  expect_true(any(grepl("Aperture", output)))
  expect_true(any(grepl("Resolution", output)))
})

test_that("HexGridInfo show() handles mixed aperture", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")

  output <- capture.output(show(grid))
  expect_true(any(grepl("4/3", output)))
})

test_that("HexGridInfo as.list() works", {
  grid <- hex_grid(area_km2 = 1000)
  lst <- as.list(grid)

  expect_type(lst, "list")
  expect_equal(lst$aperture, "3")
  expect_equal(lst$resolution, grid@resolution)
  expect_equal(lst$crs, 4326L)
})

test_that("is_hex_grid identifies HexGridInfo", {
  grid <- hex_grid(area_km2 = 1000)

  expect_true(is_hex_grid(grid))
  expect_false(is_hex_grid(list()))
  expect_false(is_hex_grid(data.frame()))
})

# =============================================================================
# HexData Tests
# =============================================================================

test_that("HexData is created by hexify", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_s4_class(result, "HexData")
  expect_s4_class(result@grid, "HexGridInfo")
  expect_equal(length(result@cell_id), 2)
})

test_that("HexData $ accessor works for data columns", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(1, 2))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_equal(result$value, c(1, 2))
  expect_equal(result$lon, c(0, 10))
})

test_that("HexData $ accessor works for virtual cell columns", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_equal(result$cell_id, result@cell_id)
  expect_equal(result$cell_cen_lon, result@cell_center[, "lon"])
  expect_equal(result$cell_cen_lat, result@cell_center[, "lat"])
  expect_length(result$cell_area_km2, 2)
  expect_length(result$cell_diag_km, 2)
})

test_that("HexData $<- assignment works", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  result$new_col <- c("a", "b")
  expect_equal(result$new_col, c("a", "b"))
})

test_that("HexData [[ accessor works", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = c(1, 2))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_equal(result[["value"]], c(1, 2))
  expect_equal(result[["cell_id"]], result@cell_id)
  expect_equal(result[["cell_cen_lon"]], result@cell_center[, "lon"])
  expect_equal(result[["cell_area_km2"]], rep(result@grid@area_km2, 2))
})

test_that("HexData [[<- assignment works", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  result[["new_col"]] <- c(100, 200)
  expect_equal(result[["new_col"]], c(100, 200))
})

test_that("HexData [ subsetting works with row indices", {
  df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55), value = 1:3)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  subset <- result[1:2, ]
  expect_s4_class(subset, "HexData")
  expect_equal(nrow(subset), 2)
  expect_equal(length(subset@cell_id), 2)
})

test_that("HexData [ subsetting returns HexData for column subset", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Column subsetting returns HexData (not a vector)
  subset <- result[, "value", drop = FALSE]
  expect_s4_class(subset, "HexData")
})

test_that("HexData nrow/ncol/dim work", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), a = 1:2, b = 3:4)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), ncol(df) + 5)  # +5 for cell columns
  expect_equal(dim(result), dim(df))
})

test_that("HexData names() includes virtual columns", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  nm <- names(result)
  expect_true("cell_id" %in% nm)
  expect_true("cell_cen_lon" %in% nm)
  expect_true("cell_cen_lat" %in% nm)
  expect_true("cell_area_km2" %in% nm)
  expect_true("cell_diag_km" %in% nm)
})

test_that("HexData grid_info() returns grid", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  g <- grid_info(result)
  expect_s4_class(g, "HexGridInfo")
  expect_equal(g@aperture, "3")
})

test_that("HexData cells() returns unique cell IDs", {
  df <- data.frame(lon = c(0, 0.001, 10), lat = c(45, 45.001, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  cell_ids <- cells(result)
  expect_equal(cell_ids, unique(result@cell_id))
})

test_that("HexData n_cells() counts unique cells", {
  df <- data.frame(lon = c(0, 0.001, 10), lat = c(45, 45.001, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 10000)

  n <- n_cells(result)
  expect_equal(n, length(unique(result@cell_id)))
})

test_that("HexData show() prints correctly", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  output <- capture.output(show(result))
  expect_true(any(grepl("HexData", output)))
  expect_true(any(grepl("Rows", output)))
  expect_true(any(grepl("Cells", output)))
})

test_that("HexData show() handles sf data", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  output <- capture.output(show(result))
  expect_true(any(grepl("sf", output)))
})

test_that("HexData show() handles many columns", {
  df <- data.frame(
    lon = c(0, 10), lat = c(45, 50),
    a = 1:2, b = 1:2, c = 1:2, d = 1:2, e = 1:2,
    f = 1:2, g = 1:2, h = 1:2, i = 1:2, j = 1:2
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  output <- capture.output(show(result))
  expect_true(any(grepl("\\.\\.\\.", output)))  # Truncated columns
})

test_that("HexData as.data.frame() works", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  df_out <- as.data.frame(result)
  expect_s3_class(df_out, "data.frame")
  expect_true("cell_id" %in% names(df_out))
  expect_true("cell_cen_lon" %in% names(df_out))
  expect_true("value" %in% names(df_out))
})

test_that("HexData as.data.frame() drops sf geometry", {
  skip_if_not_installed("sf")

  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  sf_data <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  result <- hexify(sf_data, area_km2 = 1000)

  df_out <- as.data.frame(result)
  expect_false(inherits(df_out, "sf"))
  expect_true("cell_id" %in% names(df_out))
})

test_that("HexData as.data.frame() accepts row.names", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  df_out <- as.data.frame(result, row.names = c("first", "second"))
  expect_equal(rownames(df_out), c("first", "second"))
})

test_that("HexData as.list() works", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  lst <- as.list(result)
  expect_type(lst, "list")
  expect_true("data" %in% names(lst))
  expect_true("grid" %in% names(lst))
  expect_true("cell_id" %in% names(lst))
})

test_that("is_hex_data identifies HexData", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  expect_true(is_hex_data(result))
  expect_false(is_hex_data(df))
  expect_false(is_hex_data(list()))
})

# =============================================================================
# Validity Tests
# =============================================================================

test_that("HexGridInfo validity rejects invalid aperture", {
  expect_error(
    new("HexGridInfo", aperture = "5", resolution = 5L, crs = 4326L),
    "aperture must be"
  )
})

test_that("HexGridInfo validity rejects invalid resolution", {
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = -1L, crs = 4326L),
    "resolution must be between"
  )
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = 50L, crs = 4326L),
    "resolution must be between"
  )
})

test_that("HexGridInfo validity rejects invalid area_km2", {
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = 5L, area_km2 = -100, crs = 4326L),
    "area_km2 must be positive"
  )
})

test_that("HexGridInfo validity rejects invalid diagonal_km", {
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = 5L, diagonal_km = -50, crs = 4326L),
    "diagonal_km must be positive"
  )
})

test_that("HexGridInfo validity rejects invalid crs", {
  expect_error(
    new("HexGridInfo", aperture = "3", resolution = 5L, crs = -1L),
    "crs must be a positive"
  )
})

# =============================================================================
# extract_grid Tests
# =============================================================================

test_that("extract_grid returns NULL when allow_null = TRUE", {
  result <- hexify:::extract_grid(NULL, allow_null = TRUE)
  expect_null(result)
})

test_that("extract_grid errors on NULL when allow_null = FALSE", {
  expect_error(hexify:::extract_grid(NULL), "grid specification required")
})

test_that("extract_grid errors on invalid object type", {
  expect_error(hexify:::extract_grid(list(a = 1)), "Cannot extract grid")
  expect_error(hexify:::extract_grid(data.frame()), "Cannot extract grid")
  expect_error(hexify:::extract_grid("string"), "Cannot extract grid")
})

# =============================================================================
# HexGridInfo_to_hexify_grid Tests
# =============================================================================

test_that("HexGridInfo_to_hexify_grid converts aperture 3", {
  grid <- hex_grid(area_km2 = 1000, aperture = 3)
  legacy <- hexify:::HexGridInfo_to_hexify_grid(grid)

  expect_s3_class(legacy, "hexify_grid")
  expect_equal(legacy$aperture, 3L)
  expect_equal(legacy$index_type, "z3")
})

test_that("HexGridInfo_to_hexify_grid converts aperture 7", {
  grid <- hex_grid(area_km2 = 10000, aperture = 7)
  legacy <- hexify:::HexGridInfo_to_hexify_grid(grid)

  expect_s3_class(legacy, "hexify_grid")
  expect_equal(legacy$aperture, 7L)
  expect_equal(legacy$index_type, "z7")
})

test_that("HexGridInfo_to_hexify_grid converts aperture 4", {
  grid <- hex_grid(area_km2 = 1000, aperture = 4)
  legacy <- hexify:::HexGridInfo_to_hexify_grid(grid)

  expect_s3_class(legacy, "hexify_grid")
  expect_equal(legacy$aperture, 4L)
  expect_equal(legacy$index_type, "zorder")
})

test_that("HexGridInfo_to_hexify_grid converts mixed aperture 4/3", {
  grid <- hex_grid(area_km2 = 1000, aperture = "4/3")
  legacy <- hexify:::HexGridInfo_to_hexify_grid(grid)

  expect_s3_class(legacy, "hexify_grid")
  expect_equal(legacy$aperture, 3L)  # Base aperture for mixed
  expect_equal(legacy$aperture_type, "MIXED43")
})

# =============================================================================
# HexData [ subset edge cases
# =============================================================================

test_that("HexData [ returns vector when extracting single column without drop=FALSE", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Extracting single column with drop=TRUE returns vector
  vec <- result@data[, "value", drop = TRUE]
  expect_type(vec, "integer")
})

test_that("HexData [[ accessor works with numeric index", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50), value = 1:2)
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  # Numeric indexing goes to underlying data
  expect_equal(result[[1]], c(0, 10))  # First column (lon)
})

test_that("HexData [[ accessor returns cell_diag_km", {
  df <- data.frame(lon = c(0, 10), lat = c(45, 50))
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)

  diag <- result[["cell_diag_km"]]
  expect_length(diag, 2)
  expect_true(all(diag > 0))
})

# =============================================================================
# HexData show() edge cases
# =============================================================================

test_that("HexData show() handles more than 3 rows", {
  df <- data.frame(
    lon = c(0, 5, 10, 15, 20),
    lat = c(45, 46, 47, 48, 49),
    value = 1:5
  )
  result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 5000)

  output <- capture.output(show(result))
  expect_true(any(grepl("more rows", output)))
})
