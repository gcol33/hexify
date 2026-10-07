# Agreement of hexify's grids with DGGAL (Jacovella-St-Louis et al. 2025), the
# reference implementation of the IVEA grids: van Leeuwen and Strebe's (2006)
# vertex-oriented equal-area projection on the icosahedron. ISEA3H and ISEA7H
# run alongside as the control, on the same orientation and pipeline.
#
# DGGAL places the icosahedron 0.05 degrees west of the standard ISEA
# orientation, vertex 0 at 11.20E, 58.2825N, azimuth 0, and reads coordinates
# as WGS84 geodetic, mapping them to the authalic sphere before projecting.
# dggal_zones.py lists every zone of a grid with its centroid and vertices,
# latitudes taken back to the authalic sphere, and the hexify grid is built on
# that orientation. For each zone, its centroid is assigned to a hexify cell;
# the script records whether the zones map one-to-one onto hexify's cells, the
# largest and median distance between DGGAL's centroid and hexify's centre of
# that cell, and the largest distance from a DGGAL vertex to the nearest
# hexify corner of the same cell, with how many zones have a different number
# of vertices than their hexify cell has corners.
#
# Usage: Rscript paper/bench/bench_dggal_agreement.R
# Set DGGAL_PYTHON to a Python interpreter with the dggal package (pip install
# dggal; wheels exist for Python 3.6 to 3.13) if it is not the default.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

DGGAL_PYTHON <- Sys.getenv("DGGAL_PYTHON", "python")
ZONES_PY <- file.path(dirname(sys_script_path()), "dggal_zones.py")
DGGAL_ORIENTATION <- c(11.20, 58.282525588538995, 0)
# DGGAL name, hexify projection, aperture, levels
CONFIGS <- list(
  list("ISEA3H", "isea", 3L, 1:8),
  list("IVEA3H", "ivea", 3L, 1:8),
  list("ISEA7H", "isea", 7L, 1:5),
  list("IVEA7H", "ivea", 7L, 1:5)
)

dggal_zones <- function(name, level) {
  out <- tempfile(fileext = ".csv")
  status <- system2(DGGAL_PYTHON, c(shQuote(ZONES_PY), name, level, shQuote(out)))
  if (!identical(status, 0L)) stop("dggal_zones.py failed for ", name, " level ", level)
  d <- read.csv(out)
  unlink(out)
  d
}

rows <- list()
for (cfg in CONFIGS) for (res in cfg[[4]]) {
  d <- dggal_zones(cfg[[1]], res)
  g <- hex_grid(resolution = res, aperture = cfg[[3]], projection = cfg[[2]],
                orientation = DGGAL_ORIENTATION)
  ce <- d[d$kind == "centroid", ]
  vx <- d[d$kind == "vertex", ]
  id <- lonlat_to_cell(ce$lon, ce$lat, g)
  hc <- cell_to_lonlat(id, g)
  centre_gap_m <- gc_km(ce$lon, ce$lat, hc[[1]], hc[[2]]) * 1000

  rings <- hexify:::isea_cell_rings(id, res, g@aperture, hexify:::icosa_arg(g), 0)
  vsplit <- split(vx[c("lon", "lat")], factor(vx$zone, levels = unique(ce$zone)))
  corner_gap_m <- numeric(nrow(ce))
  count_differs <- logical(nrow(ce))
  for (i in seq_len(nrow(ce))) {
    h <- unique(rings[[i]][, 1:2, drop = FALSE])
    v <- vsplit[[ce$zone[i]]]
    corner_gap_m[i] <- max(vapply(seq_len(nrow(v)), function(j)
      min(gc_km(v$lon[j], v$lat[j], h[, 1], h[, 2])), numeric(1))) * 1000
    count_differs[i] <- nrow(v) != nrow(h)
  }

  message(sprintf(paste("%s level %d: %d zones onto %d cells of %s, centre gap max %.3g m,",
                        "corner gap max %.3g m, %d vertex counts differ"),
                  cfg[[1]], res, nrow(ce), length(unique(id)), n_cells(g),
                  max(centre_gap_m), max(corner_gap_m), sum(count_differs)))
  rows[[length(rows) + 1]] <- data.frame(
    dggrs = cfg[[1]], projection = cfg[[2]], aperture = cfg[[3]], resolution = res,
    n_zones = nrow(ce), n_hexify_cells = n_cells(g), n_cells_hit = length(unique(id)),
    max_centre_gap_m = max(centre_gap_m), median_centre_gap_m = median(centre_gap_m),
    max_corner_gap_m = max(corner_gap_m), median_corner_gap_m = median(corner_gap_m),
    n_vertex_count_differs = sum(count_differs))
}

dggal_version <- tryCatch(system2(DGGAL_PYTHON, c("-c", shQuote(
  "import importlib.metadata as m; print(m.version('dggal'))")), stdout = TRUE),
  error = function(e) NA_character_)
write_result(do.call(rbind, rows), "dggal_agreement",
             extra = c(paste("DGGAL python:", DGGAL_PYTHON),
                       paste("dggal version:", dggal_version),
                       paste("orientation:", paste(DGGAL_ORIENTATION, collapse = ", "))))
