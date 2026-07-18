
# tests/testthat/test-icosahedron.R
# Tests for icosahedron construction and face detection
#
# Functions tested:
# - hexify_build_icosa()
# - hexify_face_centers()
# - hexify_which_face()

# =============================================================================
# ICOSAHEDRON CONSTRUCTION
# =============================================================================

test_that("icosahedron builds successfully with default orientation", {
  expect_no_error(hexify_build_icosa())
})

test_that("icosahedron builds with custom orientation", {
  # ISEA3H default orientation
  expect_no_error(hexify_build_icosa(11.25, 58.28252559, 0))

  # Alternative orientations
  expect_no_error(hexify_build_icosa(0, 58.28252559, 0))
  expect_no_error(hexify_build_icosa(11.25, 58.28252559, 45))
})

# =============================================================================
# FACE CENTERS
# =============================================================================

test_that("face_centers returns 20 faces with valid coordinates", {
  hexify_build_icosa()

  centers <- hexify_face_centers()

  expect_equal(nrow(centers), 20)
  expect_true("lon" %in% names(centers))
  expect_true("lat" %in% names(centers))

  # All coordinates should be finite
  expect_true(all(is.finite(centers$lon)))
  expect_true(all(is.finite(centers$lat)))

  # hexify_face_centers() returns radians; valid range is [-pi, pi] / [-pi/2, pi/2]
  expect_true(all(centers$lon >= -pi & centers$lon <= pi))
  expect_true(all(centers$lat >= -pi / 2 & centers$lat <= pi / 2))
})

test_that("face_centers returns the exact standard ISEA orientation", {
  # The 20 icosahedron face centers are a fixed, closed-form geometric
  # structure for the default ISEA3H orientation (pole_lon=11.25,
  # pole_lat=58.28252559, azimuth=0) - not just "some finite value".
  hexify_build_icosa()
  centers <- hexify_face_centers()

  expected_lon <- c(
    -1.3744467860, -0.5890486226, 0.1963495408, 0.9817477043, 1.7671458677,
    -2.1598449437, -1.0095829578, 0.1963495479, 1.4022820395, 2.5525440367,
    -1.7393106141, -0.5890486168, 0.9817477099, 2.1320096958, -2.9452431057,
    -2.1598449424, -1.3744467859, 1.7671458676, 2.5525440380, -2.9452431041
  )
  expected_lat <- c(
    1.2059324925, 0.6154797087, 0.3648638281, 0.6154797087, 1.2059325048,
    0.6154797060, 0.0000000000, -0.3648638281, 0.0000000000, 0.6154797113,
    -0.0000000000, -0.6154797113, -0.6154797060, -0.0000000000, 0.3648638281,
    -0.6154797054, -1.2059324973, -1.2059325001, -0.6154797120, -0.3648638281
  )

  expect_equal(centers$lon, expected_lon, tolerance = 1e-8)
  expect_equal(centers$lat, expected_lat, tolerance = 1e-8)
})

test_that("face centers span the globe", {
  hexify_build_icosa()
  centers <- hexify_face_centers()

  # Should have faces in both hemispheres
  expect_true(any(centers$lat > 0))
  expect_true(any(centers$lat < 0))

  # Should have faces across longitude range
  # Note: coordinates may be in radians (range ~[-pi, pi]) or degrees
  lon_range <- max(centers$lon) - min(centers$lon)
  expect_true(lon_range > 2.0)  # Works for both radians (~5.5 rad) and degrees (>100°)
})

# =============================================================================
# FACE DETECTION (which_face)
# =============================================================================

test_that("which_face returns valid face numbers", {
  hexify_build_icosa()

  test_points <- data.frame(
    lon = c(0, 90, -90, 180, 0, 0),
    lat = c(0, 0, 0, 0, 45, -45)
  )

  for (i in seq_len(nrow(test_points))) {
    face <- hexify_which_face(test_points$lon[i], test_points$lat[i])
    expect_true(face >= 0 && face <= 19,
                info = sprintf("lon=%.2f, lat=%.2f",
                              test_points$lon[i], test_points$lat[i]))
  }
})

test_that("which_face is deterministic", {
  hexify_build_icosa()

  lon <- 16.37
  lat <- 48.21

  face1 <- hexify_which_face(lon, lat)
  face2 <- hexify_which_face(lon, lat)

  expect_equal(face1, face2)
})

test_that("which_face handles extreme coordinates", {
  hexify_build_icosa()

  # Near poles
  expect_true(hexify_which_face(0, 89.9999) %in% 0:19)
  expect_true(hexify_which_face(0, -89.9999) %in% 0:19)

  # Near antimeridian
  expect_true(hexify_which_face(179.9999, 0) %in% 0:19)
  expect_true(hexify_which_face(-179.9999, 0) %in% 0:19)
})

test_that("which_face matches reference data", {
  hexify_build_icosa()

  truth <- read.csv(testthat::test_path("data/truth_faces.csv"))
  expect_true(all(c("lon", "lat", "tnum") %in% names(truth)))

  mismatches <- 0L
  msgs <- character()

  for (i in seq_len(nrow(truth))) {
    got <- hexify_which_face(truth$lon[i], truth$lat[i])
    if (got != truth$tnum[i]) {
      mismatches <- mismatches + 1L
      if (mismatches <= 5L) {
        msgs <- c(msgs, sprintf(
          "row %d: got %d vs expected %d (lon=%.9f lat=%.9f)",
          i, got, truth$tnum[i], truth$lon[i], truth$lat[i]
        ))
      }
    }
  }

  if (mismatches > 0L) {
    message(paste(msgs, collapse = "\n"))
  }
  expect_equal(mismatches, 0L)
})

# =============================================================================
# FACE CENTER CONSISTENCY
# =============================================================================

test_that("face center detection returns the originating face", {
  # hexify_face_centers() returns radians; hexify_which_face() expects
  # degrees (see its @param docs) - convert before round-tripping.
  hexify_build_icosa()
  centers <- hexify_face_centers()

  for (face in 0:19) {
    lon_deg <- centers$lon[face + 1] * 180 / pi
    lat_deg <- centers$lat[face + 1] * 180 / pi
    detected <- hexify_which_face(lon_deg, lat_deg)
    expect_equal(detected, face,
                info = sprintf("Face center %d was detected as face %d", face, detected))
  }
})
