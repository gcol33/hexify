# OGC zone identifiers and the WGS84 earth model against DGGAL
# (Jacovella-St-Louis et al. 2025), the implementation behind OGC's ISEA3H,
# ISEA7H, IVEA3H and IVEA7H DGGRS definitions.
#
# For every zone of each grid and level, dggal_zone_ids.py lists DGGAL's
# textZIRS and uint64ZIRS identifiers and its centroid in WGS84 geodetic
# latitude. The hexify grid is OGC's: hex_grid(ellipsoid = "WGS84",
# orientation = "ogc"). Each centroid is assigned to a hexify cell; the
# script records whether the zones map one-to-one onto the cells, the largest
# distance between DGGAL's centroid and the cell's centre, and how many cells
# cell_to_index() names as DGGAL does and index_to_cell() reads back to the
# same cell, in both identifier forms.
#
# Usage: Rscript paper/bench/bench_dggal_zirs.R
# Set DGGAL_PYTHON to a Python interpreter with the dggal package.

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

DGGAL_PYTHON <- Sys.getenv("DGGAL_PYTHON", "python")
IDS_PY <- file.path(dirname(sys_script_path()), "dggal_zone_ids.py")
# DGGAL name, hexify projection, aperture, levels
CONFIGS <- list(
  list("ISEA3H", "isea", 3L, 0:10),
  list("IVEA3H", "ivea", 3L, 0:10),
  list("ISEA7H", "isea", 7L, 0:5),
  list("IVEA7H", "ivea", 7L, 0:5)
)

dggal_ids <- function(name, level) {
  out <- tempfile(fileext = ".csv")
  status <- system2(DGGAL_PYTHON, c(shQuote(IDS_PY), name, level, shQuote(out)))
  if (!identical(status, 0L)) stop("dggal_zone_ids.py failed for ", name, " level ", level)
  d <- read.csv(out, colClasses = c(uint64 = "character"))
  unlink(out)
  d
}

rows <- list()
for (cfg in CONFIGS) for (res in cfg[[4]]) {
  d <- dggal_ids(cfg[[1]], res)
  g <- hex_grid(resolution = res, aperture = cfg[[3]], projection = cfg[[2]],
                ellipsoid = "WGS84", orientation = "ogc")
  cell <- lonlat_to_cell(d$lon, d$lat, g)
  ctr <- cell_to_lonlat(cell, g)
  text <- cell_to_index(cell, g, form = "textZIRS")
  u64 <- as.character(cell_to_index(cell, g, form = "uint64ZIRS"))
  back_text <- index_to_cell(d$text, g, form = "textZIRS")
  back_u64 <- index_to_cell(d$uint64, g, form = "uint64ZIRS")
  row <- data.frame(
    dggrs = cfg[[1]], level = res, n_zones = nrow(d), n_hexify_cells = n_cells(g),
    n_cells_hit = length(unique(cell)),
    max_centre_gap_m = max(gc_km(d$lon, d$lat, ctr$lon_deg, ctr$lat_deg)) * 1000,
    text_equal = sum(text == d$text), uint64_equal = sum(u64 == d$uint64),
    text_read_back = sum(back_text == cell, na.rm = TRUE),
    uint64_read_back = sum(back_u64 == cell, na.rm = TRUE))
  print(row)
  rows[[length(rows) + 1]] <- row
}

dggal_version <- tryCatch(system2(DGGAL_PYTHON, c("-c", shQuote(
  "import importlib.metadata as m; print(m.version('dggal'))")), stdout = TRUE),
  error = function(e) NA_character_)
write_result(do.call(rbind, rows), "dggal_zirs",
             extra = c(paste("DGGAL python:", DGGAL_PYTHON),
                       paste("dggal version:", dggal_version)))
