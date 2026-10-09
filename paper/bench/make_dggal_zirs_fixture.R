# DGGAL's zone identifiers for the test suite (tests/testthat/data/dggal_zirs.csv),
# which runs without DGGAL: every zone of ISEA3H at levels 0 to 3 and of
# ISEA7H at levels 0 to 2, random samples at ISEA3H level 8 (with zone
# E6-317-A), ISEA7H level 3, IVEA3H level 5 and IVEA7H level 3, each with its
# textZIRS and uint64ZIRS identifiers and its centroid in WGS84 geodetic
# latitude (dggal_zone_ids.py).
#
# Usage: Rscript paper/bench/make_dggal_zirs_fixture.R
# Set DGGAL_PYTHON to a Python interpreter with the dggal package.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                             value = TRUE))))
DGGAL_PYTHON <- Sys.getenv("DGGAL_PYTHON", "python")
SEED <- 20261009L

zones <- function(name, level, n = NULL, keep = character(0)) {
  out <- tempfile(fileext = ".csv")
  status <- system2(DGGAL_PYTHON, c(shQuote(file.path(script_dir, "dggal_zone_ids.py")),
                                    name, level, shQuote(out)))
  if (!identical(status, 0L)) stop("dggal_zone_ids.py failed for ", name, " level ", level)
  d <- read.csv(out, colClasses = c(uint64 = "character"))
  unlink(out)
  if (!is.null(n)) d <- d[sort(union(sample(nrow(d), n), which(d$text %in% keep))), ]
  cbind(dggrs = name, level = level, d)
}

set.seed(SEED)
out <- rbind(
  do.call(rbind, lapply(0:3, function(l) zones("ISEA3H", l))),
  zones("ISEA3H", 8L, 60L, keep = "E6-317-A"),
  do.call(rbind, lapply(0:2, function(l) zones("ISEA7H", l))),
  zones("ISEA7H", 3L, 200L),
  zones("IVEA3H", 5L, 100L),
  zones("IVEA7H", 3L, 100L)
)
out$lon <- sprintf("%.13f", out$lon)
out$lat <- sprintf("%.13f", out$lat)
write.csv(out, file.path(script_dir, "..", "..", "tests", "testthat", "data", "dggal_zirs.csv"),
          row.names = FALSE, quote = FALSE)
