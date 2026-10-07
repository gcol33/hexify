# DGGAL's IVEA zones for the test suite (tests/testthat/data/dggal_ivea.csv),
# which runs without DGGAL: every zone of IVEA3H at level 3 and a random
# subset of IVEA7H at level 3, each with its centroid and vertices, latitudes
# on the authalic sphere (dggal_zones.py).
#
# Usage: Rscript paper/bench/make_dggal_fixture.R
# Set DGGAL_PYTHON to a Python interpreter with the dggal package.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                             value = TRUE))))
DGGAL_PYTHON <- Sys.getenv("DGGAL_PYTHON", "python")
N_AP7 <- 120L
SEED <- 20261008L

zones <- function(name, level) {
  out <- tempfile(fileext = ".csv")
  status <- system2(DGGAL_PYTHON, c(shQuote(file.path(script_dir, "dggal_zones.py")),
                                    name, level, shQuote(out)))
  if (!identical(status, 0L)) stop("dggal_zones.py failed for ", name, " level ", level)
  d <- read.csv(out)
  unlink(out)
  cbind(dggrs = name, level = level, d)
}

ap3 <- zones("IVEA3H", 3L)
ap7 <- zones("IVEA7H", 3L)
set.seed(SEED)
keep <- sample(unique(ap7$zone), N_AP7)
ap7 <- ap7[ap7$zone %in% keep, ]
out <- rbind(ap3, ap7)
out$lon <- sprintf("%.11f", out$lon)
out$lat <- sprintf("%.11f", out$lat)
write.csv(out, file.path(script_dir, "..", "..", "tests", "testthat", "data", "dggal_ivea.csv"),
          row.names = FALSE, quote = FALSE)
