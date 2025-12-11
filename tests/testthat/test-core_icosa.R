test_that("core icosa builds and which_face matches reference", {
  hexify_build_icosa()

  centers <- hexify_face_centers()
  expect_equal(nrow(centers), 20)
  for (i in seq_len(nrow(centers))) {
    expect_true(is.finite(centers$lon[i]), info = paste("center lon face", i-1))
    expect_true(is.finite(centers$lat[i]), info = paste("center lat face", i-1))
  }

  truth <- read.csv(testthat::test_path("data/truth_faces.csv"))
  expect_true(all(c("lon","lat","tnum") %in% names(truth)))

  bad <- 0L
  msgs <- character()
  for (i in seq_len(nrow(truth))) {
    got <- hexify_which_face(truth$lon[i], truth$lat[i])
    if (got != truth$tnum[i]) {
      bad <- bad + 1L
      if (bad <= 10L) {
        msgs <- c(msgs, sprintf("row %d: got %d vs %d (lon=%.9f lat=%.9f)",
                                i, got, truth$tnum[i], truth$lon[i], truth$lat[i]))
      }
    }
  }
  if (bad > 0L) message(paste(msgs, collapse = "\n"))
  expect_identical(bad, 0L)
})
