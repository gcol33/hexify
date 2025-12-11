# test-integer_seqnum.R
# Tests for integer sequential numbering functions

test_that("lonlat_to_seqnum returns positive integers for aperture 3", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)
  seqnums <- hexify_lonlat_to_seqnum(lons, lats, resolution = 10, aperture = 3)

  expect_true(all(seqnums > 0))
  expect_true(is.numeric(seqnums))
  expect_true(all(is.finite(seqnums)))
})

test_that("lonlat_to_seqnum returns positive integers for aperture 4", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)
  seqnums <- hexify_lonlat_to_seqnum(lons, lats, resolution = 10, aperture = 4)

  expect_true(all(seqnums > 0))
  expect_true(is.numeric(seqnums))
  expect_true(all(is.finite(seqnums)))
})

test_that("lonlat_to_seqnum returns positive integers for aperture 7", {
  lons <- c(0, 10, 20, -30, 45)
  lats <- c(0, 45, -30, 15, -60)
  seqnums <- hexify_lonlat_to_seqnum(lons, lats, resolution = 5, aperture = 7)

  expect_true(all(seqnums > 0))
  expect_true(is.numeric(seqnums))
  expect_true(all(is.finite(seqnums)))
})

test_that("round-trip conversion works for aperture 3", {
  test_lons <- c(0, 5, 10, -30, 45, 90, -90)
  test_lats <- c(0, 45, -30, 15, -60, 30, -15)
  resolution <- 10
  aperture <- 3
  tolerance <- 1.0  # degrees

  seqnums <- hexify_lonlat_to_seqnum(test_lons, test_lats, resolution, aperture)
  coords <- hexify_seqnum_to_lonlat(seqnums, resolution, aperture)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])

    # Account for longitude wrapping
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

test_that("round-trip conversion works for aperture 4", {
  test_lons <- c(0, 5, 10, -30, 45)
  test_lats <- c(0, 45, -30, 15, -60)
  resolution <- 10
  aperture <- 4
  tolerance <- 1.0

  seqnums <- hexify_lonlat_to_seqnum(test_lons, test_lats, resolution, aperture)
  coords <- hexify_seqnum_to_lonlat(seqnums, resolution, aperture)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

test_that("round-trip conversion works for aperture 7", {
  test_lons <- c(0, 5, 10, -30, 45)
  test_lats <- c(0, 45, -30, 15, -60)
  resolution <- 5
  aperture <- 7
  # Aperture 7 uses surrogate-substrate conversion with rounding,
  # which can introduce up to ~4 degrees of error at res 5
  tolerance <- 4.0

  seqnums <- hexify_lonlat_to_seqnum(test_lons, test_lats, resolution, aperture)
  coords <- hexify_seqnum_to_lonlat(seqnums, resolution, aperture)

  for (i in seq_along(test_lons)) {
    lon_error <- abs(coords$lon_deg[i] - test_lons[i])
    lat_error <- abs(coords$lat_deg[i] - test_lats[i])
    if (lon_error > 180) lon_error <- 360 - lon_error

    expect_lt(lon_error, tolerance)
    expect_lt(lat_error, tolerance)
  }
})

test_that("seqnum_to_cell returns valid cell info", {
  seqnums <- c(1702, 1954, 100)
  info <- hexify_seqnum_to_cell(seqnums, resolution = 5, aperture = 3)

  expect_true(is.data.frame(info))
  expect_true("face" %in% names(info))
  expect_true("i" %in% names(info))
  expect_true("j" %in% names(info))

  # Faces should be 0-19
  expect_true(all(info$face >= 0 & info$face <= 19))

  # i, j should be non-negative
  expect_true(all(info$i >= 0))
  expect_true(all(info$j >= 0))
})

test_that("resolution 0 returns face numbers (1-20)", {
  # At resolution 0, there are only 20 cells (one per face)
  lons <- c(0, 60, 120, -60, -120, 0)
  lats <- c(45, 45, 45, 45, 45, -45)

  seqnums <- hexify_lonlat_to_seqnum(lons, lats, resolution = 0, aperture = 3)

  # Seqnums should be 1-20
  expect_true(all(seqnums >= 1 & seqnums <= 20))
})

test_that("extreme coordinates work correctly", {
  extreme_coords <- data.frame(
    lon = c(-180, 180, 0, 0),
    lat = c(0, 0, -89, 89)
  )

  for (i in 1:nrow(extreme_coords)) {
    seqnum <- hexify_lonlat_to_seqnum(
      extreme_coords$lon[i],
      extreme_coords$lat[i],
      resolution = 10,
      aperture = 3
    )

    expect_true(is.finite(seqnum))
    expect_true(seqnum > 0)
  }
})

test_that("batch processing handles large datasets", {
  set.seed(123)
  n <- 1000
  lons <- runif(n, -180, 180)
  lats <- runif(n, -90, 90)

  seqnums <- hexify_lonlat_to_seqnum(lons, lats, resolution = 10, aperture = 3)

  expect_length(seqnums, n)
  expect_true(all(seqnums > 0))

  # Check reasonable uniqueness (should be high for random points)
  uniqueness <- length(unique(seqnums)) / n
  expect_gt(uniqueness, 0.8)
})
