# hexify grids in the area and shape comparison of Kmoch et al. (2022).
#
# Kmoch, Vasilyev, Virro and Uuemaa (2022, Big Earth Data 6: 256-275) measured
# normalized area and the isoperimetric quotient 4 pi A / p^2 of every cell of
# ten open-source DGGS. Their pipeline (p. 262): cell polygons from the
# library in WGS84, cells crossing the antimeridian removed, each cell
# projected into a Lambert azimuthal equal-area plane on the WGS84 ellipsoid
# centred on the cell, area and perimeter measured in that plane. This script
# runs the same pipeline on hexify's cells, drawn as DGGRID draws them (corners
# joined by straight lines in the plane), and sets the result beside their
# per-resolution tables (data supplement, Zenodo record 6634479). For ISEA7H,
# FULLER7H, ISEA3H, FULLER3H and H3 it also compares every cell with their
# per-cell file, matched by cell ID; aperture-7 cells at odd resolutions whose
# DGGRID corners differ from hexify's (the corner bug in CLAUDE.md) are counted
# apart.
#
# A pentagon's walls bend where they cross face edges, so its corners joined by
# straight lines miss part of it. Kmoch et al. dropped cells with invalid
# geometry, and their smallest normalized areas show the pentagons kept at
# some resolutions (about 0.85) and dropped at others (above 0.97), so the hex_
# columns give the comparison without pentagons as well.
#
# Beside the pipeline's numbers it gives hexify's exact ones: area and
# perimeter on the sphere with the walls followed as the grid draws them
# (cell_metrics()).
#
# Usage: Rscript paper/bench/bench_kmoch_distortion.R
# Downloads about 220 MB to paper/bench/data/kmoch2022/ on first run.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

ZENODO <- "https://zenodo.org/records/6634479/files/"
DATA_DIR <- file.path(dirname(results_dir()), "data", "kmoch2022")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

fetch <- function(file) {
  path <- file.path(DATA_DIR, file)
  if (!file.exists(path)) {
    utils::download.file(paste0(ZENODO, file, "?download=1"), path, mode = "wb", quiet = TRUE)
  }
  path
}

# WGS84 (EPSG:7030)
A_WGS84 <- 6378137
F_WGS84 <- 1 / 298.257223563
E2 <- F_WGS84 * (2 - F_WGS84)
E <- sqrt(E2)

# Snyder (1987, eqs. 3-12, 14-15, 24-27): q, authalic latitude and the oblique
# ellipsoidal Lambert azimuthal equal-area projection centred on (lon0, lat0),
# vectorised over points that each carry their own centre.
q_of <- function(phi) {
  s <- sin(phi)
  (1 - E2) * (s / (1 - E2 * s^2) - log((1 - E * s) / (1 + E * s)) / (2 * E))
}
QP <- q_of(pi / 2)
RQ <- A_WGS84 * sqrt(QP / 2)

laea_xy <- function(lon, lat, lon0, lat0) {
  d2r <- pi / 180
  phi <- lat * d2r; phi1 <- lat0 * d2r; dl <- (lon - lon0) * d2r
  beta <- asin(pmax(-1, pmin(1, q_of(phi) / QP)))
  beta1 <- asin(pmax(-1, pmin(1, q_of(phi1) / QP)))
  m1 <- cos(phi1) / sqrt(1 - E2 * sin(phi1)^2)
  D <- A_WGS84 * m1 / (RQ * cos(beta1))
  B <- RQ * sqrt(2 / (1 + sin(beta1) * sin(beta) + cos(beta1) * cos(beta) * cos(dl)))
  cbind(x = B * D * cos(beta) * sin(dl),
        y = (B / D) * (cos(beta1) * sin(beta) - sin(beta1) * cos(beta) * cos(dl)))
}

# The closed form against PROJ on a few points, before anything is measured
check_laea <- function() {
  set.seed(7)
  lon0 <- 16.37; lat0 <- 48.21
  lon <- lon0 + runif(20, -3, 3); lat <- lat0 + runif(20, -3, 3)
  ours <- laea_xy(lon, lat, lon0, lat0)
  proj <- sf::sf_project("EPSG:4326",
                         sprintf("+proj=laea +lon_0=%s +lat_0=%s +ellps=WGS84 +units=m", lon0, lat0),
                         cbind(lon, lat))
  err <- max(abs(ours - proj))
  if (err > 1e-3) stop(sprintf("closed-form LAEA differs from PROJ by %.3g m", err))
  err
}

# Per-cell planar area and perimeter of corner polygons, each in its own LAEA
corner_metrics <- function(ids, g) {
  rings <- hexify:::cell_corner_rings(ids, g)
  n <- vapply(rings, nrow, integer(1))
  owner <- rep(seq_along(ids), n)
  lon <- unlist(lapply(rings, function(r) r[, 1]), use.names = FALSE)
  lat <- unlist(lapply(rings, function(r) r[, 2]), use.names = FALSE)
  span <- tapply(lon, owner, function(v) diff(range(v)))
  crossed <- as.vector(span > 180)

  ctr <- cell_to_lonlat(ids, g)
  xy <- laea_xy(lon, lat, ctr$lon_deg[owner], ctr$lat_deg[owner])
  nxt <- unlist(lapply(split(seq_along(owner), owner), function(k) c(k[-1], k[1])),
                use.names = FALSE)
  cross <- xy[, 1] * xy[nxt, 2] - xy[nxt, 1] * xy[, 2]
  seg <- sqrt((xy[nxt, 1] - xy[, 1])^2 + (xy[nxt, 2] - xy[, 2])^2)
  area <- abs(as.vector(tapply(cross, owner, sum))) / 2
  perimeter <- as.vector(tapply(seg, owner, sum))
  data.frame(cell_id = as.character(ids), crossed = crossed, area = area,
             perimeter = perimeter, ipq = 4 * pi * area / perimeter^2)
}

GRIDS <- list(
  list(name = "ISEA3H", kmoch = "DGGRID_ISEA3H", res = 2:6, make = function(r) hex_grid(resolution = r, aperture = 3)),
  list(name = "FULLER3H", kmoch = "DGGRID_FULLER3H", res = 2:6, make = function(r) hex_grid(resolution = r, aperture = 3, projection = "fuller")),
  list(name = "ISEA4H", kmoch = NA, res = 2:7, make = function(r) hex_grid(resolution = r, aperture = 4)),
  list(name = "FULLER4H", kmoch = NA, res = 2:7, make = function(r) hex_grid(resolution = r, aperture = 4, projection = "fuller")),
  list(name = "ISEA43H", kmoch = NA, res = 2:9, make = function(r) hex_grid(resolution = r, aperture = "4/3")),
  list(name = "ISEA7H", kmoch = "DGGRID_ISEA7H", res = 2:5, make = function(r) hex_grid(resolution = r, aperture = 7)),
  list(name = "FULLER7H", kmoch = "DGGRID_FULLER7H", res = 2:5, make = function(r) hex_grid(resolution = r, aperture = 7, projection = "fuller")),
  list(name = "H3", kmoch = "h3", res = 2:4, make = function(r) hex_grid(resolution = r, type = "h3"))
)

# Their per-cell file: a cell whose DGGRID corners lie more than 1 m from
# hexify's is one of DGGRID's misplaced corners, not a difference in method.
compare_cells <- function(ours, g, spec, res) {
  theirs <- sf::st_read(fetch(sprintf("%s_%d_LAEA_step2.fgb", spec$kmoch, res)), quiet = TRUE)
  theirs <- theirs[!theirs$crossed, ]
  key <- as.character(theirs$cell_id)
  at <- match(key, ours$cell_id)
  ok <- !is.na(at) & !ours$crossed[at]
  theirs <- theirs[ok, ]; at <- at[ok]

  open_ring <- function(r) if (all(r[1, ] == r[nrow(r), ])) r[-nrow(r), , drop = FALSE] else r
  rings <- lapply(hexify:::cell_corner_rings(
    if (hexify:::is_h3_grid(g)) key[ok] else hexify:::as_cell_id(key[ok]), g), open_ring)
  their_xy <- lapply(sf::st_geometry(theirs), function(p) open_ring(p[[1]]))
  corner_gap <- mapply(function(a, b) {
    if (nrow(a) != nrow(b)) return(Inf)
    max(vapply(seq_len(nrow(a)), function(i) {
      min(gc_km(a[i, 1], a[i, 2], b[, 1], b[, 2]))
    }, numeric(1))) * 1000
  }, rings, their_xy)
  same <- corner_gap <= 1
  rel_area <- abs(theirs$area / ours$area[at] - 1)
  d_ipq <- abs(theirs$ipq - ours$ipq[at])
  data.frame(grid = spec$name, resolution = res, cells_compared = length(at),
             cells_corners_differ = sum(!same),
             max_rel_area_diff_same_corners = max(rel_area[same]),
             max_ipq_diff_same_corners = max(d_ipq[same]),
             max_rel_area_diff_corners_differ = if (any(!same)) max(rel_area[!same]) else NA,
             max_ipq_diff_corners_differ = if (any(!same)) max(d_ipq[!same]) else NA)
}

laea_error_m <- check_laea()
message(sprintf("closed-form LAEA agrees with PROJ to %.2g m", laea_error_m))

rows <- list(); cells <- list()
for (spec in GRIDS) {
  for (res in spec$res) {
    g <- spec$make(res)
    ids <- hexify:::grid_cells(g)
    t0 <- proc.time()[["elapsed"]]
    m <- corner_metrics(ids, g)
    keep <- !m$crossed
    na <- m$area[keep] / mean(m$area[keep])
    hex <- !is_pentagon(ids, g)[keep]
    exact <- cell_metrics(ids, g)
    row <- data.frame(
      grid = spec$name, resolution = res, n_cells = length(ids),
      n_crossed = sum(m$crossed),
      norm_area_sd = sd(na), norm_area_min = min(na), norm_area_max = max(na),
      ipq_mean = mean(m$ipq[keep]), ipq_sd = sd(m$ipq[keep]),
      ipq_min = min(m$ipq[keep]), ipq_max = max(m$ipq[keep]),
      hex_norm_area_sd = sd(na[hex]), hex_norm_area_min = min(na[hex]),
      hex_norm_area_max = max(na[hex]), hex_ipq_mean = mean(m$ipq[keep][hex]),
      exact_norm_area_sd = sd(exact$normalized_area),
      exact_norm_area_min = min(exact$normalized_area),
      exact_norm_area_max = max(exact$normalized_area),
      exact_ipq_mean = mean(exact$ipq), exact_ipq_sd = sd(exact$ipq),
      kmoch_n_cells = NA, kmoch_norm_area_sd = NA, kmoch_norm_area_min = NA,
      kmoch_norm_area_max = NA, kmoch_ipq_mean = NA, kmoch_ipq_sd = NA,
      seconds = proc.time()[["elapsed"]] - t0)
    if (!is.na(spec$kmoch)) {
      k <- utils::read.csv(fetch(sprintf("%s_%d_stats.csv", spec$kmoch, res)))
      row$kmoch_n_cells <- k$num_cells
      row$kmoch_norm_area_sd <- k$norm_area_std
      row$kmoch_norm_area_min <- k$norm_area_min
      row$kmoch_norm_area_max <- k$norm_area_max
      row$kmoch_ipq_mean <- k$ipq_mean
      row$kmoch_ipq_sd <- k$ipq_std
      cells[[length(cells) + 1]] <- compare_cells(m, g, spec, res)
      print(cells[[length(cells)]], digits = 3)
    }
    rows[[length(rows) + 1]] <- row
    print(row[, 1:11], digits = 4)
  }
}

out <- do.call(rbind, rows)
write_result(out, "kmoch_distortion",
             extra = sprintf("closed-form LAEA vs PROJ, max |dx|, |dy|: %.3g m", laea_error_m))
write_result(do.call(rbind, cells), "kmoch_percell")
