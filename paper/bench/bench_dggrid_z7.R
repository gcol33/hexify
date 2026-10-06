# hexify's aperture-7 Z7 index against DGGRID's, the reference implementation
# of IGEO7's Z7 index (Kmoch et al. 2025). For every cell of each grid both
# programs write the cell's Z7 string; the script counts the cells whose
# strings agree, the cells each program reads back from its own string, and
# the cells whose hexify parent is the cell DGGRID's string names once its last
# digit is dropped. Grids: ISEA at resolutions 0 to MAX_RES in the standard
# orientation, and resolutions 0 to 5 under a given icosahedron orientation and
# on the FULLER projection.
#
# Usage: Rscript paper/bench/bench_dggrid_z7.R
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_dggrid_common.R"))

MAX_RES <- 7L
AP7 <- CONFIGS[[3]]
ORIENTATION <- c(-40, 20, 33)

# placement, projection, orientation or NULL, resolutions
GRIDS <- list(
  list("standard", "isea", NULL, 0:MAX_RES),
  list("vert0", "isea", ORIENTATION, 0:5),
  list("standard", "fuller", NULL, 0:5)
)

rows <- list()
for (gr in GRIDS) for (res in gr[[4]]) {
  icosa <- c(if (!is.null(gr[[3]])) orient_lines(gr[[3]]), proj_lines(gr[[2]]))
  g <- hex_grid(resolution = res, aperture = 7, projection = gr[[2]],
                orientation = if (is.null(gr[[3]])) "standard" else gr[[3]])
  ids <- seq_len(n_cells(g))
  h <- cell_to_index(ids, g)
  d <- dggrid_z7(AP7, res, seqnum = ids, icosa = icosa)
  d_back <- dggrid_z7(AP7, res, z7 = d, icosa = icosa)
  h_back <- hexify:::isea_index_to_cells(h, 7L, "z7", hexify:::icosa_arg(g))
  n_parent <- NA_integer_
  if (res > 0) {
    d_parent <- hexify:::isea_index_to_cells(substr(d, 1, nchar(d) - 1L), 7L, "z7",
                                             hexify:::icosa_arg(g))
    n_parent <- sum(get_parent(ids, g) == d_parent)
  }
  message(sprintf("%-8s %-6s res %d: %d cells, %d strings agree, %d parents agree",
                  gr[[1]], gr[[2]], res, length(ids), sum(h == d), n_parent))
  rows[[length(rows) + 1]] <- data.frame(
    placement = gr[[1]], projection = gr[[2]], resolution = res, n_cells = length(ids),
    n_identical = sum(h == d), n_distinct_hexify = length(unique(h)),
    n_distinct_dggrid = length(unique(d)),
    n_roundtrip_hexify = sum(h_back == ids), n_roundtrip_dggrid = sum(d_back == ids),
    n_parent_agree = n_parent)
}

dg_version <- tryCatch(system2("git", c("-C", shQuote(dirname(dirname(dirname(dirname(dirname(DGGRID_EXE)))))),
                                        "log", "-1", "--format=%h"), stdout = TRUE),
                       error = function(e) NA)
write_result(do.call(rbind, rows), "dggrid_z7",
             extra = c(paste("orientation:", paste(ORIENTATION, collapse = " ")),
                       paste("DGGRID executable:", DGGRID_EXE),
                       paste("DGGRID commit:", dg_version)))
