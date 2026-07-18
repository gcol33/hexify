# tests/testthat/test-globe.R
# Tests for globe_centers, plot_globe(), and the internal resolve_center()
# helper (R/plot_globe.R). Previously had zero test coverage.

# =============================================================================
# globe_centers
# =============================================================================

test_that("globe_centers is a named list of lon/lat presets", {
  expect_type(globe_centers, "list")
  expect_true(all(c("europe", "north_america", "south_america", "africa",
                     "asia", "oceania", "middle_east", "south_asia",
                     "pacific", "caribbean", "arctic", "antarctic")
                   %in% names(globe_centers)))

  for (nm in names(globe_centers)) {
    preset <- globe_centers[[nm]]
    expect_length(preset, 2)
    expect_equal(names(preset), c("lon", "lat"), info = nm)
    expect_true(preset["lon"] >= -180 && preset["lon"] <= 180, info = nm)
    expect_true(preset["lat"] >= -90 && preset["lat"] <= 90, info = nm)
  }
})

test_that("globe_centers$europe matches documented value", {
  expect_equal(unname(globe_centers$europe), c(10, 50))
})

# =============================================================================
# resolve_center() (internal helper used by plot_globe())
# =============================================================================

test_that("resolve_center accepts a valid preset name", {
  expect_equal(hexify:::resolve_center("europe"), globe_centers$europe)
})

test_that("resolve_center errors on an unknown preset name", {
  expect_error(hexify:::resolve_center("nowhere"), "Unknown center preset")
})

test_that("resolve_center accepts numeric c(lon, lat), named or not", {
  expect_equal(hexify:::resolve_center(c(lon = 5, lat = 40)), c(lon = 5, lat = 40))
  expect_equal(hexify:::resolve_center(c(5, 40)), c(lon = 5, lat = 40))
})

test_that("resolve_center errors on invalid input", {
  expect_error(hexify:::resolve_center(c(1, 2, 3)), "preset name or numeric")
  expect_error(hexify:::resolve_center(TRUE), "preset name or numeric")
})

# =============================================================================
# plot_globe()
# =============================================================================

test_that("plot_globe(return_data = TRUE) returns valid, back-side-filtered geometry", {
  skip_on_cran()  # Slow sf/orthographic-projection operations
  skip_if_not_installed("sf")

  data <- plot_globe(area = 500000, center = "europe", return_data = TRUE)

  expect_type(data, "list")
  expect_true("hexagons" %in% names(data))
  expect_true("ocean_circle" %in% names(data))
  expect_s3_class(data$hexagons, "sf")

  # Back-side hexagons are filtered out: fewer cells than a full global grid
  full_grid <- hex_grid(area_km2 = 500000)
  full_count <- 10 * as.integer(full_grid@aperture)^full_grid@resolution + 2
  expect_true(nrow(data$hexagons) < full_count)

  # Invalid geometries are repaired (st_buffer(0)) before being returned
  expect_true(all(sf::st_is_valid(data$hexagons)))
})

test_that("plot_globe(return_data = TRUE) works with a custom numeric center", {
  skip_on_cran()  # Slow sf/orthographic-projection operations
  skip_if_not_installed("sf")

  data <- plot_globe(area = 500000, center = c(0, 0), return_data = TRUE)
  expect_s3_class(data$hexagons, "sf")
  expect_true(nrow(data$hexagons) > 0)
})

test_that("plot_globe errors on an unknown center preset", {
  skip_if_not_installed("sf")
  expect_error(plot_globe(area = 500000, center = "nowhere"), "Unknown center preset")
})
