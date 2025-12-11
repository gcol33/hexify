test_that("Snyder forward matches truth_txty.csv", {
  # Initialize to ISEA3H orientation
  hexify_build_icosa(11.25, 58.28252559, 0)
  hexify_set_precision("ultra")

  truth <- read.csv(testthat::test_path("data/truth_txty.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("lon","lat","tnum","tx","ty") %in% names(truth)))

  atol <- 2e-8
  bad  <- 0L
  msgs <- character()

  for (i in seq_len(nrow(truth))) {
    out  <- hexify_forward(truth$lon[i], truth$lat[i])
    face <- as.integer(out[["face"]])
    px   <- as.numeric(out[["tx"]])
    py   <- as.numeric(out[["ty"]])

    ok <- (face == as.integer(truth$tnum[i])) &&
      (abs(px - truth$tx[i]) <= atol) &&
      (abs(py - truth$ty[i]) <= atol)

    if (!ok) {
      bad <- bad + 1L
      if (bad <= 10L) {
        msgs <- c(msgs, sprintf(
          "Mismatch lon=%.10f lat=%.10f: face %d vs %d, tx %.12g vs %.12g, ty %.12g vs %.12g",
          truth$lon[i], truth$lat[i], face, truth$tnum[i], px, truth$tx[i], py, truth$ty[i]
        ))
      }
    }
  }

  if (bad > 0L) message(paste(msgs, collapse = "\n"))
  expect_identical(bad, 0L)
})
