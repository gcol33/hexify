# DGGRID's cells and centres under non-standard orientations and on the
# Fuller projection, for the test suite (tests/testthat/data/dggrid_reference.csv),
# which runs without DGGRID. Random points are assigned to cells by DGGRID for
# icosahedra placed at a given vertex 0 and azimuth and about region centres,
# on the ISEA and FULLER projections, and each cell's centre is read back.
# DGGRID's aperture-7 Z7 strings go to tests/testthat/data/dggrid_z7.csv.
#
# Usage: Rscript paper/bench/make_dggrid_fixture.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

N_POINTS <- 150L
SEED <- 20261005L

# placement, aperture as hex_grid() takes it, CONFIGS entry, resolution,
# orientation or NULL, region centre or NULL, projection
CASES <- list(
  list("vert0", "3", CONFIGS[[1]], 6, c(-40, 20, 33), NULL, "isea"),
  list("vert0", "7", CONFIGS[[3]], 4, c(0, 90, 0), NULL, "isea"),
  list("vert0", "4/3", CONFIGS[[4]], 6, c(150, -60, 200), NULL, "isea"),
  list("region", "4", CONFIGS[[2]], 5, NULL, c(10, 46.5), "isea"),
  list("region", "3", CONFIGS[[1]], 5, NULL, c(0, 90), "isea"),
  list("standard", "3", CONFIGS[[1]], 6, NULL, NULL, "fuller"),
  list("standard", "4", CONFIGS[[2]], 5, NULL, NULL, "fuller"),
  list("standard", "7", CONFIGS[[3]], 4, NULL, NULL, "fuller"),
  list("standard", "7", CONFIGS[[3]], 5, NULL, NULL, "fuller"),
  list("standard", "4/3", CONFIGS[[4]], 6, NULL, NULL, "fuller"),
  list("vert0", "4", CONFIGS[[2]], 5, c(-40, 20, 33), NULL, "fuller")
)

set.seed(SEED)
out <- do.call(rbind, lapply(CASES, function(cs) {
  lon <- runif(N_POINTS, -180, 180)
  lat <- asin(runif(N_POINTS, -1, 1)) * 180 / pi
  icosa <- c(if (cs[[1]] != "standard") orient_lines(cs[[5]], cs[[6]]),
             proj_lines(cs[[7]]))
  seqnum <- dggrid_seqnum(cs[[3]], cs[[4]], lon, lat, icosa)
  ctr <- dggrid_centres(cs[[3]], cs[[4]], seqnum, icosa)
  data.frame(
    placement = cs[[1]], projection = cs[[7]], aperture = cs[[2]], resolution = cs[[4]],
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
  "tests", "testthat", "data", "dggrid_reference.csv")
write.csv(out, dest, row.names = FALSE, quote = FALSE)
message("wrote ", nrow(out), " rows to ", dest)

# DGGRID's Z7 strings of every aperture-7 cell at resolutions 0 to 3, and of
# Z7_SAMPLE random cells at resolution 7 (tests/testthat/data/dggrid_z7.csv).
Z7_SAMPLE <- 2000L
z7_cells <- c(lapply(0:3, function(res) list(res = res, seqnum = seq_len(10 * 7^res + 2))),
              list(list(res = 7L, seqnum = sort(sample.int(10 * 7^7 + 2, Z7_SAMPLE)))))
z7 <- do.call(rbind, lapply(z7_cells, function(cs) {
  data.frame(resolution = cs$res, seqnum = format(cs$seqnum, scientific = FALSE, trim = TRUE),
             z7 = dggrid_z7(CONFIGS[[3]], cs$res, seqnum = cs$seqnum))
}))
dest_z7 <- file.path(dirname(dest), "dggrid_z7.csv")
write.csv(z7, dest_z7, row.names = FALSE, quote = FALSE)
message("wrote ", nrow(z7), " rows to ", dest_z7)
