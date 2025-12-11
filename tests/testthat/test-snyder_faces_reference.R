test_that("which_face matches truth_faces.csv", {
  hexify_build_icosa(11.25, 58.28252559, 0)
  truth <- read.csv(testthat::test_path("data/truth_faces.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("lon","lat","tnum") %in% names(truth)))
  for (i in seq_len(nrow(truth))) {
    f <- hexify_which_face(truth$lon[i], truth$lat[i])
    expect_identical(as.integer(f), as.integer(truth$tnum[i]))
  }
})
