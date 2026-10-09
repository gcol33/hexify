# Does binning commute with the hierarchy? Griffin (2026, Sec. 12) bins 10,000
# points at a fine resolution and takes each cell's ancestor k levels up, bins
# the same points directly at the coarse resolution, and counts how often the
# two cells differ. A nested hierarchy gives 0; hexagons do not nest, so a
# point in a cell that reaches out of its parent can land in a neighbour of
# that parent.
#
# One level up, the lattice gives the rate exactly. Aperture 3: two thirds of
# the area lies in cells centred on parent corners, each a third in each of
# three parents, so 2/3 * 2/3 = 4/9. Aperture 4: three quarters in cells
# centred on parent edges, half in each of two parents, 3/4 * 1/2 = 3/8.
# Aperture 7: six sevenths in cells off the parent's centre, 1/12 of each in a
# neighbour of the parent, 6/7 * 1/12 = 1/14. Hex9 (aperture 9, Griffin's
# rule: the parent holds the cell's mode-0 half): three of nine cells
# straddle their parent, half in a neighbour, 3/9 * 1/2 = 1/6; Griffin
# predicts (1/6) 3^(1 - k) for k levels, as straddling cells k levels down
# sit in a band along the ancestor's edge that thins by 3 a level. With
# LIBHEX9_DUMP pointing at the hex9_libhex9_dump tool, libhex9's own rate
# (hex9_cell_parent, k times) is measured on the same points. Snyder's projection is
# equal-area, so the rates on the sphere are those of the plane, up to the
# cells beside the twelve vertices. Fuller's projection and H3 are not
# equal-area, so their rates follow the lattice only approximately.
#
# Usage: Rscript paper/bench/bench_hierarchy_commutation.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

N_POINTS <- 10000
pts <- sphere_points(N_POINTS, seed = 2026)

GRIDS <- list(
  list(name = "ISEA3H", fine = 12, make = function(r) hex_grid(resolution = r, aperture = 3), one = 4 / 9),
  list(name = "ISEA4H", fine = 10, make = function(r) hex_grid(resolution = r, aperture = 4), one = 3 / 8),
  list(name = "ISEA7H", fine = 8, make = function(r) hex_grid(resolution = r, aperture = 7), one = 1 / 14),
  list(name = "FULLER7H", fine = 8, make = function(r) hex_grid(resolution = r, aperture = 7, projection = "fuller"), one = 1 / 14),
  # The family refines by 4 for the first floor(r / 2) levels, so resolution 12
  # holds one more aperture-4 level than resolution 11: one level up is a step
  # of aperture 4
  list(name = "ISEA43H", fine = 12, make = function(r) hex_grid(resolution = r, aperture = "4/3"), one = 3 / 8),
  list(name = "ISEA[4,7]H", fine = 8, make = function(r) hex_grid(resolution = r, aperture = c(4, 7)), one = 1 / 14),
  list(name = "H3", fine = 8, make = function(r) hex_grid(resolution = r, type = "h3"), one = 1 / 14),
  list(name = "Hex9", fine = 10, make = function(r) hex_grid(resolution = r, aperture = 9), one = 1 / 6)
)

rows <- list()
for (spec in GRIDS) {
  g <- spec$make(spec$fine)
  fine <- lonlat_to_cell(pts$lon, pts$lat, g)
  for (k in 1:6) {
    coarse <- spec$make(spec$fine - k)
    up <- get_parent(fine, g, levels = k)
    rate <- mean(up != lonlat_to_cell(pts$lon, pts$lat, coarse))
    rows[[length(rows) + 1]] <- data.frame(
      grid = spec$name, fine_resolution = spec$fine, levels_up = k,
      disagree = rate, se = sqrt(rate * (1 - rate) / N_POINTS),
      lattice_one_level = if (k == 1) spec$one else NA)
    print(rows[[length(rows)]], digits = 4)
  }
}

dump <- Sys.getenv("LIBHEX9_DUMP")
if (nzchar(dump)) {
  input <- tempfile(fileext = ".txt")
  writeLines(sprintf("%.17g %.17g", pts$lon, pts$lat), input)
  ref <- read.csv(text = system2(dump, c("commute", 4, 6), stdin = input, stdout = TRUE),
                  colClasses = "character")
  for (k in 1:6) {
    rate <- mean(ref[[paste0("k", k)]] != ref$direct)
    rows[[length(rows) + 1]] <- data.frame(
      grid = "Hex9 (libhex9)", fine_resolution = 4 + k, levels_up = k,
      disagree = rate, se = sqrt(rate * (1 - rate) / N_POINTS),
      lattice_one_level = if (k == 1) 1 / 6 else NA)
    print(rows[[length(rows)]], digits = 4)
  }
}

write_result(do.call(rbind, rows), "hierarchy_commutation",
             extra = sprintf("points: %d uniform on the sphere, seed 2026", N_POINTS))
