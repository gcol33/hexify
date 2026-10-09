# Aperture-7 neighbours by IGEO7 digit arithmetic against the quad-lattice
# path that sends every quad-leaving step through lon/lat (the inverse and the
# forward projection). A cell whose six lattice steps stay in its quad takes
# them in both paths; only cells on a quad edge differ. For each even
# resolution from 2 to 20 the script times both paths on N_EDGE cells along the
# edges of the ten diamond quads and on N_POINTS cells of random points, and
# counts the cells whose neighbours the two paths agree on. At
# resolutions 5 and 6 it also times and compares every cell of the grid.
#
# Usage: Rscript paper/bench/bench_z7_neighbors.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

N_EDGE <- 20000L
N_POINTS <- 20000L
REPS <- 5L
MIN_SECONDS <- 0.25
SEED <- 20261009L

neighbours <- function(ids, res, lonlat) {
  f <- if (lonlat) hexify:::cpp_get_neighbors_isea_lonlat else hexify:::cpp_get_neighbors_isea
  f(numeric(0), ids, res, 7L, integer(0))
}

# Median seconds per cell over REPS runs, each repeating the call until it
# has run MIN_SECONDS, well above the timer's resolution
per_cell <- function(ids, res, lonlat) {
  t <- vapply(seq_len(REPS), function(k) {
    calls <- 1L
    repeat {
      elapsed <- system.time(for (m in seq_len(calls)) neighbours(ids, res, lonlat))[["elapsed"]]
      if (elapsed >= MIN_SECONDS) return(elapsed / calls)
      calls <- calls * 2L
    }
  }, numeric(1))
  stats::median(t) / length(ids)
}

agree <- function(ids, res) {
  a <- neighbours(ids, res, FALSE)
  b <- neighbours(ids, res, TRUE)
  sum(vapply(seq_along(a), function(k) identical(a[[k]], b[[k]]), logical(1)))
}

# Cells along the four edges of the diamond quads: at even resolutions the
# stored (i, j) of a cell is its substrate coordinate, and the quad's cells
# fill [0, S) x [0, S), S the quad edge in substrate steps.
edge_cells <- function(res, n) {
  S <- hexify:::cpp_quad_edge_dim(res, 7L, integer(0))
  n_diamonds <- 10L
  quad <- sample.int(n_diamonds, n, replace = TRUE)
  along <- floor(stats::runif(n, 0, S))
  side <- sample.int(4L, n, replace = TRUE)
  i <- ifelse(side == 1L, 0, ifelse(side == 2L, S - 1, along))
  j <- ifelse(side == 3L, 0, ifelse(side == 4L, S - 1, along))
  unique(hexify:::cpp_quad_ij_to_cell(numeric(0), quad, i, j, res, 7L, integer(0)))
}

set.seed(SEED)
rows <- list()
for (res in seq(2L, 20L, by = 2L)) {
  g <- hex_grid(resolution = res, aperture = 7)
  edge <- edge_cells(res, N_EDGE)
  pts <- sphere_points(N_POINTS, SEED + res)
  interior <- unique(lonlat_to_cell(pts$lon, pts$lat, g))
  rows[[length(rows) + 1]] <- data.frame(
    resolution = res, cells = "quad edge", n = length(edge),
    us_per_cell_z7 = 1e6 * per_cell(edge, res, FALSE),
    us_per_cell_lonlat = 1e6 * per_cell(edge, res, TRUE),
    n_agree = agree(edge, res))
  rows[[length(rows) + 1]] <- data.frame(
    resolution = res, cells = "random points", n = length(interior),
    us_per_cell_z7 = 1e6 * per_cell(interior, res, FALSE),
    us_per_cell_lonlat = 1e6 * per_cell(interior, res, TRUE),
    n_agree = agree(interior, res))
  message(sprintf("res %2d: edge %.3f vs %.3f us, agree %d/%d", res,
                  rows[[length(rows) - 1]]$us_per_cell_z7,
                  rows[[length(rows) - 1]]$us_per_cell_lonlat,
                  rows[[length(rows) - 1]]$n_agree, length(edge)))
}
for (res in 5:6) {
  ids <- bit64::as.integer64(seq_len(10 * 7^res + 2))
  rows[[length(rows) + 1]] <- data.frame(
    resolution = res, cells = "all", n = length(ids),
    us_per_cell_z7 = 1e6 * per_cell(ids, res, FALSE),
    us_per_cell_lonlat = 1e6 * per_cell(ids, res, TRUE),
    n_agree = agree(ids, res))
}

write_result(do.call(rbind, rows), "z7_neighbors",
             extra = c(paste("repetitions:", REPS), paste("seconds per timing:", MIN_SECONDS),
                       paste("seed:", SEED)))
