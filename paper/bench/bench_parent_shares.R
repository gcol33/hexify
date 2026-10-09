# How a cell divides among the coarser cells it overlaps. On the plane lattice
# of one refinement step the shares are 1/3 (aperture 3, a cell centred on a
# parent corner), 1/2 (aperture 4, on a parent edge) and 11/12 and 1/12
# (aperture 7); hex_aggregate() uses them on every equal-area grid. This
# script measures every share of every cell on the sphere the way
# hex_aggregate() does on Fuller's projection and H3: the cell and each coarser
# cell get_parent(overlapping = TRUE) returns are clipped against each other
# where both are straight (an ISEA-family cell on each face plane, an H3 cell
# on the gnomonic plane) and the pieces measured on the sphere. On the
# equal-area grids the measured shares should equal the lattice's to within
# the measuring error; on Fuller's projection and on H3 they depart from it.
#
# As an independent check, the shares are measured again with s2: the cells
# as spherical polygons from cell_to_sf() (ISEA-family walls densified to
# S2_DENSIFY, H3 cells from their corners, exact), intersected by s2.
#
# For each grid it also aggregates a uniform density (each cell's own area)
# one level up under both rules and reports how far each parent's total falls
# from the parent's area.
#
# Usage: Rscript paper/bench/bench_parent_shares.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

GRIDS <- list(
  list(name = "ISEA3H", res = 3, make = function(r) hex_grid(resolution = r, aperture = 3)),
  list(name = "ISEA3H", res = 5, make = function(r) hex_grid(resolution = r, aperture = 3)),
  list(name = "ISEA4H", res = 2, make = function(r) hex_grid(resolution = r, aperture = 4)),
  list(name = "ISEA4H", res = 4, make = function(r) hex_grid(resolution = r, aperture = 4)),
  list(name = "ISEA7H", res = 2, make = function(r) hex_grid(resolution = r, aperture = 7)),
  list(name = "ISEA7H", res = 3, make = function(r) hex_grid(resolution = r, aperture = 7)),
  list(name = "ISEA43H", res = 5, make = function(r) hex_grid(resolution = r, aperture = "4/3")),
  list(name = "ISEA[4,7,3,7]H", res = 4, make = function(r) hex_grid(resolution = r, aperture = c(4, 7, 3, 7)[seq_len(r)])),
  list(name = "IVEA7H", res = 3, make = function(r) hex_grid(resolution = r, aperture = 7, projection = "ivea")),
  list(name = "ISEA3H octahedron", res = 4, make = function(r) hex_grid(resolution = r, aperture = 3, polyhedron = "octahedron")),
  list(name = "FULLER3H", res = 4, make = function(r) hex_grid(resolution = r, aperture = 3, projection = "fuller")),
  list(name = "FULLER4H", res = 3, make = function(r) hex_grid(resolution = r, aperture = 4, projection = "fuller")),
  list(name = "FULLER7H", res = 3, make = function(r) hex_grid(resolution = r, aperture = 7, projection = "fuller")),
  list(name = "H3", res = 2, make = function(r) hex_grid(resolution = r, type = "h3")),
  list(name = "H3", res = 3, make = function(r) hex_grid(resolution = r, type = "h3"))
)

all_cells <- function(g) hexify:::grid_cells(g)

S2_DENSIFY <- 1e-5

s2_shares <- function(ids, ov, g, pg) {
  flat <- if (g@grid_type == "h3") unlist(ov, use.names = FALSE) else hexify:::cell_id_unlist(ov)
  parents <- unique(flat)
  child <- rep(seq_along(ids), lengths(ov))
  geog <- function(x, grid) {
    s2::as_s2_geography(sf::st_geometry(cell_to_sf(
      x, grid, wrap_dateline = FALSE,
      densify = if (grid@grid_type == "h3") NULL else S2_DENSIFY)))
  }
  piece <- s2::s2_area(s2::s2_intersection(geog(ids, g)[child],
                                           geog(parents, pg)[match(flat, parents)]),
                       radius = 1)
  piece / as.numeric(rowsum(piece, child, reorder = TRUE))[child]
}

rows <- list()
for (spec in GRIDS) {
  g <- spec$make(spec$res)
  pg <- spec$make(spec$res - 1)
  ids <- all_cells(g)
  t0 <- Sys.time()
  ov <- get_parent(ids, g, overlapping = TRUE)
  n_par <- lengths(ov)
  measured <- hexify:::sphere_shares(ids, ov, g)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  lattice <- if (g@grid_type == "h3") rep(NA_real_, length(measured)) else
    hexify:::lattice_shares(n_par, hexify:::step_aperture(g))
  split_share <- rep(n_par, n_par) > 1
  dev <- abs(measured - lattice)[split_share]
  s2_dev <- abs(measured - s2_shares(ids, ov, g, pg))

  density <- unname(cell_area(ids, g))
  rel <- function(rule) {
    up <- hex_aggregate(ids, density, g, rule = rule)
    stopifnot(nrow(up) == as.numeric(n_cells(pg)))
    max(abs(up$value / unname(cell_area(up$cell_id, pg)) - 1))
  }
  rows[[length(rows) + 1]] <- data.frame(
    grid = spec$name, resolution = spec$res, cells = length(ids),
    cells_split = sum(n_par > 1),
    parents_per_cell = paste(sprintf("%d:%d", sort(unique(n_par)),
                                     as.integer(table(n_par))), collapse = " "),
    max_share_dev = if (all(is.na(dev))) NA else max(dev),
    median_share_dev = if (all(is.na(dev))) NA else median(dev),
    min_split_share = min(measured[split_share]),
    max_s2_dev = max(s2_dev),
    uniform_area_rule_max_rel = rel("area"),
    uniform_centre_rule_max_rel = rel("centre"),
    measure_seconds = secs)
  print(rows[[length(rows)]], digits = 4)
}

write_result(do.call(rbind, rows), "parent_shares",
             extra = sprintf("s2 check: cell_to_sf densify %g", S2_DENSIFY))
