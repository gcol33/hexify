
# tests/testthat/test-assign.R
# Tests for hexify_assign function

setup_icosa <- function() {
  cpp_build_icosa()
}

# =============================================================================
# BASIC FUNCTIONALITY
# =============================================================================

test_that("hexify_assign returns correct structure", {
  setup_icosa()

  lon <- c(0, 45, -120)
  lat <- c(0, 30, -45)

  result <- hexify_assign(lon, lat, effective_res = 3)

  expect_s3_class(result, "data.frame")
  expect_named(result, c(
    "id", "face", "effective_res", "center_lon", "center_lat"
  ))
  expect_equal(nrow(result), 3)
})

test_that("hexify_assign produces valid cell IDs", {
  setup_icosa()

  lon <- c(16.37, 2.35, -3.70)
  lat <- c(48.21, 48.86, 40.42)

  result <- hexify_assign(lon, lat, effective_res = 5)

  # Z3 index: 2-digit quad followed by one digit per resolution level
  expect_true(all(grepl("^\\d{2}[0-2]{5}$", result$id)))

  # Face is the quad (0-11)
  expect_true(all(result$face >= 0))
  expect_true(all(result$face <= 11))
})

test_that("hexify_assign agrees with the grid pipeline", {
  setup_icosa()

  lon <- c(16.37, 2.35, -3.70, 151.21, -58.38)
  lat <- c(48.21, 48.86, 40.42, -33.87, -34.60)

  for (res in 1:8) {
    grid <- hex_grid(resolution = res, aperture = 3)
    cells <- lonlat_to_cell(lon, lat, grid)
    centers <- cell_to_lonlat(cells, grid)

    result <- hexify_assign(lon, lat, effective_res = res)

    expect_equal(result$center_lon, centers$lon_deg)
    expect_equal(result$center_lat, centers$lat_deg)
    expect_equal(result$id, unname(cell_to_index(cells, grid)))
  }
})

test_that("hexify_assign centers are valid coordinates", {
  setup_icosa()

  lon <- c(0, 90, -90, 180, -180)
  lat <- c(0, 45, -45, 60, -60)

  result <- hexify_assign(lon, lat, effective_res = 4)

  # Centers should be valid lon/lat
  expect_true(all(result$center_lon >= -180 & result$center_lon <= 180))
  expect_true(all(result$center_lat >= -90 & result$center_lat <= 90))
})

# =============================================================================
# RESOLUTION HANDLING
# =============================================================================

test_that("hexify_assign handles different resolutions", {
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  for (res in c(1, 3, 5, 7)) {
    result <- hexify_assign(lon, lat, effective_res = res)

    expect_equal(result$effective_res, res)
    expect_equal(nrow(result), 1)
  }
})

test_that("hexify_assign resolution affects cell ID length", {
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  result_low <- hexify_assign(lon, lat, effective_res = 2)
  result_high <- hexify_assign(lon, lat, effective_res = 6)

  # Higher resolution should produce longer digit sequences
  # Extract digit part of ID
  digits_low <- sub(".*:Z3:", "", result_low$id)
  digits_high <- sub(".*:Z3:", "", result_high$id)

  expect_lt(nchar(digits_low), nchar(digits_high))
})

# =============================================================================
# CENTER ACCURACY
# =============================================================================

test_that("hexify_assign centers lie inside the assigned cell", {
  setup_icosa()

  # Uniform points on the sphere
  set.seed(23)
  n <- 400
  lon <- 360 * runif(n) - 180
  lat <- asin(2 * runif(n) - 1) * 180 / pi

  gc_km <- function(lon1, lat1, lon2, lat2) {
    r <- pi / 180
    EARTH_RADIUS_KM * acos(pmin(1, pmax(-1,
      sin(lat1 * r) * sin(lat2 * r) +
        cos(lat1 * r) * cos(lat2 * r) * cos((lon2 - lon1) * r))))
  }

  for (res in 1:8) {
    result <- hexify_assign(lon, lat, effective_res = res)
    d <- gc_km(lon, lat, result$center_lon, result$center_lat)

    # Center spacing of a grid of 10 * 3^res + 2 cells over the sphere; a point
    # lies at most one circumradius, about 0.6 of the spacing, from its cell
    # center.
    n_cells <- 10 * 3^res + 2
    spacing <- sqrt(2 * EARTH_SURFACE_KM2 / (sqrt(3) * n_cells))
    expect_lt(max(d), 0.7 * spacing)
  }
})

# =============================================================================
# POLYGON GENERATION
# =============================================================================

test_that("hexify_assign with make_polygons=TRUE returns sf", {
  skip_if_not_installed("sf")
  setup_icosa()


  lon <- c(10, 20)
  lat <- c(50, 55)

  result <- hexify_assign(lon, lat, effective_res = 3, make_polygons = TRUE)

  expect_s3_class(result, "sf")
  expect_true("geometry" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("hexify_assign polygons have 6 corners", {
  skip_if_not_installed("sf")
  setup_icosa()

  lon <- c(10)
  lat <- c(50)

  result <- hexify_assign(lon, lat, effective_res = 3, make_polygons = TRUE)

  # Get coordinates of polygon
  coords <- sf::st_coordinates(result$geometry[[1]])

  # Polygon should have 7 points (6 corners + closing point)
  expect_equal(nrow(coords), 7)
})

# =============================================================================
# INPUT VALIDATION
# =============================================================================

test_that("hexify_assign validates input lengths", {
  setup_icosa()

  expect_error(
    hexify_assign(c(0, 1), c(0), effective_res = 3),
    "length"
  )
})

test_that("hexify_assign validates resolution", {
  setup_icosa()

  expect_error(
    hexify_assign(0, 0, effective_res = 0),
    "effective_res"
  )

  expect_error(
    hexify_assign(0, 0, effective_res = c(1, 2)),
    "effective_res"
  )
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("hexify_assign handles poles", {
  setup_icosa()

  # North pole
  result_north <- hexify_assign(0, 90, effective_res = 3)
  expect_equal(nrow(result_north), 1)
  expect_true(is.finite(result_north$center_lat))

  # South pole
  result_south <- hexify_assign(0, -90, effective_res = 3)
  expect_equal(nrow(result_south), 1)
  expect_true(is.finite(result_south$center_lat))
})

test_that("hexify_assign handles date line", {
  setup_icosa()

  # Points near date line
  lon <- c(179, -179, 180, -180)
  lat <- c(0, 0, 45, -45)

  result <- hexify_assign(lon, lat, effective_res = 3)
  expect_equal(nrow(result), 4)
})
