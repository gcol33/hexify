# test-h3-native.R
# Tests for native H3 C backend wrappers

# =============================================================================
# latLngToCell / cellToLatLng round-trip
# =============================================================================

test_that("cpp_h3_latLngToCell returns valid hex strings", {
  cells <- cpp_h3_latLngToCell(c(16.37, 2.35), c(48.21, 48.86), 5L)
  expect_true(is.character(cells))
  expect_equal(length(cells), 2)
  expect_true(all(nchar(cells) == 15))  # H3 hex strings are 15 chars
})

test_that("cpp_h3_latLngToCell handles NA inputs", {
  cells <- cpp_h3_latLngToCell(c(16.37, NA), c(48.21, 48.86), 5L)
  expect_false(is.na(cells[1]))
  expect_true(is.na(cells[2]))
})

test_that("latLngToCell → cellToLatLng round-trip is close", {
  lon <- c(16.37, -73.99, 139.69)
  lat <- c(48.21, 40.75, 35.69)
  cells <- cpp_h3_latLngToCell(lon, lat, 8L)
  centers <- cpp_h3_cellToLatLng(cells)

  for (i in seq_along(lon)) {
    expect_true(abs(centers$lon[i] - lon[i]) < 0.1)
    expect_true(abs(centers$lat[i] - lat[i]) < 0.1)
  }
})

# =============================================================================
# isValidCell
# =============================================================================

test_that("cpp_h3_isValidCell validates correctly", {
  cells <- cpp_h3_latLngToCell(c(0, 10), c(45, 50), 5L)
  valid <- cpp_h3_isValidCell(cells)
  expect_true(all(valid))

  invalid <- cpp_h3_isValidCell(c("garbage", "0000000000000000"))
  expect_false(any(invalid))
})

test_that("cpp_h3_isValidCell handles NA", {
  result <- cpp_h3_isValidCell(c("851e15b7fffffff", NA))
  expect_true(result[1])
  expect_true(is.na(result[2]))
})

# =============================================================================
# cellToParent / cellToChildren
# =============================================================================

test_that("cpp_h3_cellToParent returns coarser cell", {
  cell <- cpp_h3_latLngToCell(0, 45, 8L)
  parent <- cpp_h3_cellToParent(cell, 5L)
  expect_true(is.character(parent))
  expect_equal(length(parent), 1)
  expect_false(parent == cell)
})

test_that("cpp_h3_cellToChildren returns 7 children", {
  cell <- cpp_h3_latLngToCell(0, 45, 5L)
  children <- cpp_h3_cellToChildren(cell, 6L)
  expect_true(is.list(children))
  expect_equal(length(children[[1]]), 7)
  expect_true(all(cpp_h3_isValidCell(children[[1]])))
})

# =============================================================================
# cellToBoundary
# =============================================================================

test_that("cpp_h3_cellToBoundary returns closed rings", {
  cell <- cpp_h3_latLngToCell(0, 45, 5L)
  bndry <- cpp_h3_cellToBoundary(cell)
  ring <- bndry[[1]]
  expect_true(is.matrix(ring))
  expect_equal(ncol(ring), 2)
  expect_equal(ring[1, ], ring[nrow(ring), ])  # closed ring
  expect_true(nrow(ring) == 7)  # 6 hex verts + 1 closing
})

# =============================================================================
# polygonToCells
# =============================================================================

test_that("cpp_h3_polygonToCells fills a bounding box", {
  bbox <- matrix(c(16, 48, 17, 48, 17, 49, 16, 49, 16, 48),
                 ncol = 2, byrow = TRUE)
  cells <- cpp_h3_polygonToCells(bbox, 5L)
  expect_true(length(cells) > 0)
  expect_true(all(cpp_h3_isValidCell(cells)))
})

# =============================================================================
# cellAreaKm2
# =============================================================================

test_that("cpp_h3_cellAreaKm2 returns positive areas", {
  cells <- cpp_h3_latLngToCell(c(0, 0), c(0, 80), 5L)
  areas <- cpp_h3_cellAreaKm2(cells)
  expect_true(all(areas > 0))
  # Areas should differ (H3 is not equal-area)
  expect_false(areas[1] == areas[2])
})

test_that("cpp_h3_cellAreaKm2 magnitude is reasonable for res 5", {
  cell <- cpp_h3_latLngToCell(0, 45, 5L)
  area <- cpp_h3_cellAreaKm2(cell)
  # H3 res 5 avg area is ~252.9 km^2
  expect_true(area > 200 && area < 300)
})
