# DGGRID's IGEO7 Z7 indices in both of its forms, the digit string and the
# packed 64-bit index (DGGRID's INT64 output, 16 hexadecimal digits), for the
# test suite (tests/testthat/data/dggrid_z7_int.csv): every cell at
# resolutions 0 to 3, and the cells of N_POINTS random points at resolutions
# 5, 10, 15 and 20. For the points, DGGRID writes the SEQNUM and both Z7 forms
# from the point itself; the SEQNUMs are kept as the digit strings DGGRID
# writes, so the test reads them exactly at every resolution. (DGGRID built on
# Windows misreads a SEQNUM above 2^32 given as input, so the points do not
# pass through SEQNUM input.)
#
# Usage: Rscript paper/bench/make_dggrid_z7_int_fixture.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

N_POINTS <- 200L
SEED <- 20261009L
AP7 <- CONFIGS[[3]]

# DGGRID's Z7 index of the cell holding each point, in hier_ndx_form `form`
dggrid_z7_of_points <- function(res, lon, lat, form) {
  inp <- file.path(work, "z7_pts.txt"); outp <- file.path(work, "z7_pts_out.txt")
  write.table(data.frame(sprintf("%.12f", lon), sprintf("%.12f", lat)), inp,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  run_dggrid(c("dggrid_operation TRANSFORM_POINTS", dggs_lines(AP7, res),
               paste("input_file_name", inp), "input_address_type GEO",
               "input_delimiter \" \"", paste("output_file_name", outp),
               "output_address_type HIERNDX", "output_hier_ndx_system Z7",
               paste("output_hier_ndx_form", form), "output_delimiter \" \""))
  readLines(outp)
}

all_cells <- do.call(rbind, lapply(0:3, function(res) {
  seqnum <- format(seq_len(10 * 7^res + 2), scientific = FALSE, trim = TRUE)
  data.frame(resolution = res, seqnum = seqnum,
             z7 = dggrid_z7(AP7, res, seqnum = seqnum),
             z7_hex = dggrid_z7(AP7, res, seqnum = seqnum, form = "INT64"))
}))

set.seed(SEED)
point_cells <- do.call(rbind, lapply(c(5L, 10L, 15L, 20L), function(res) {
  lon <- runif(N_POINTS, -180, 180)
  lat <- asin(runif(N_POINTS, -1, 1)) * 180 / pi
  unique(data.frame(resolution = res,
                    seqnum = dggrid_seqnum(AP7, res, lon, lat, exact = TRUE),
                    z7 = dggrid_z7_of_points(res, lon, lat, "DIGIT_STRING"),
                    z7_hex = dggrid_z7_of_points(res, lon, lat, "INT64")))
}))

out <- rbind(all_cells, point_cells)

dest <- file.path(dirname(dirname(dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))))),
  "tests", "testthat", "data", "dggrid_z7_int.csv")
write.csv(out, dest, row.names = FALSE, quote = FALSE)
message("wrote ", nrow(out), " rows to ", dest)
