# DGGRID's cells and centres under non-standard orientations, for the test
# suite (tests/testthat/data/dggrid_orientation.csv), which runs without
# DGGRID. Random points are assigned to cells by DGGRID for icosahedra placed
# at a given vertex 0 and azimuth and about region centres, and each cell's
# centre is read back.
#
# Usage: Rscript paper/bench/make_dggrid_orientation_fixture.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

N_POINTS <- 150L
SEED <- 20261005L

# placement, aperture as hex_grid() takes it, CONFIGS entry, resolution,
# orientation or NULL, region centre or NULL
CASES <- list(
  list("vert0", "3", CONFIGS[[1]], 6, c(-40, 20, 33), NULL),
  list("vert0", "7", CONFIGS[[3]], 4, c(0, 90, 0), NULL),
  list("vert0", "4/3", CONFIGS[[4]], 6, c(150, -60, 200), NULL),
  list("region", "4", CONFIGS[[2]], 5, NULL, c(10, 46.5)),
  list("region", "3", CONFIGS[[1]], 5, NULL, c(0, 90))
)

set.seed(SEED)
out <- do.call(rbind, lapply(CASES, function(cs) {
  lon <- runif(N_POINTS, -180, 180)
  lat <- asin(runif(N_POINTS, -1, 1)) * 180 / pi
  orient <- orient_lines(cs[[5]], cs[[6]])
  seqnum <- dggrid_seqnum(cs[[3]], cs[[4]], lon, lat, orient)
  ctr <- dggrid_centres(cs[[3]], cs[[4]], seqnum, orient)
  data.frame(
    placement = cs[[1]], aperture = cs[[2]], resolution = cs[[4]],
    vert0_lon = if (is.null(cs[[5]])) NA else cs[[5]][1],
    vert0_lat = if (is.null(cs[[5]])) NA else cs[[5]][2],
    azimuth = if (is.null(cs[[5]])) NA else cs[[5]][3],
    region_lon = if (is.null(cs[[6]])) NA else cs[[6]][1],
    region_lat = if (is.null(cs[[6]])) NA else cs[[6]][2],
    lon = sprintf("%.12f", lon), lat = sprintf("%.12f", lat),
    seqnum = format(seqnum, scientific = FALSE, trim = TRUE),
    centre_lon = sprintf("%.12f", ctr$lon), centre_lat = sprintf("%.12f", ctr$lat)
  )
}))

dest <- file.path(dirname(dirname(dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))))),
  "tests", "testthat", "data", "dggrid_orientation.csv")
write.csv(out, dest, row.names = FALSE, quote = FALSE)
message("wrote ", nrow(out), " rows to ", dest)
