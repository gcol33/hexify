# DGGRID helpers for the agreement script. DGGRID runs
# as its command-line program, built from source (github.com/sahrk/DGGRID).
# Set DGGRID_EXE to the dggrid executable if it is not at the default path.

DGGRID_EXE <- Sys.getenv("DGGRID_EXE",
  "C:/Users/Gilles Colling/Documents/dev/upstream/DGGRID/build/src/apps/dggrid/dggrid.exe")
stopifnot(file.exists(DGGRID_EXE))

# label, hexify aperture argument, DGGRID aperture lines, resolutions. DGGRID
# reads two sequence digits beyond the resolution, so the sequence is longer
# than the deepest resolution tested; hexify takes its first `res` entries.
CONFIGS <- list(
  list("3", 3, c("dggs_aperture_type PURE", "dggs_aperture 3"), c(3, 7, 11, 15)),
  list("4", 4, c("dggs_aperture_type PURE", "dggs_aperture 4"), c(2, 5, 8, 11)),
  list("7", 7, c("dggs_aperture_type PURE", "dggs_aperture 7"), c(2, 4, 6, 8)),
  list("4/3", "4/3", NULL, c(4, 8, 12)),
  list("4,3,4,7,4,7", c(4, 3, 4, 7, 4, 7, 4, 3),
       c("dggs_aperture_type SEQUENCE", "dggs_aperture_sequence 43474743"), c(3, 6))
)

work <- file.path(tempdir(), "dggrid_agreement")
dir.create(work, showWarnings = FALSE)

run_dggrid <- function(lines) {
  meta <- file.path(work, "run.meta")
  writeLines(c(lines, "verbosity 0"), meta)
  out <- system2(DGGRID_EXE, shQuote(meta), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("dggrid failed:\n", paste(out, collapse = "\n"))
  invisible(out)
}

# `orient` is further DGGRID lines placing the icosahedron (see
# orient_lines()); NULL keeps DGGRID's standard orientation.
dggs_lines <- function(cfg, res, orient = NULL) {
  ap_lines <- if (identical(cfg[[1]], "4/3")) {
    c("dggs_aperture_type MIXED43", sprintf("dggs_num_aperture_4_res %d", res %/% 2))
  } else cfg[[3]]
  c("dggs_type CUSTOM", "dggs_topology HEXAGON", "dggs_proj ISEA", ap_lines,
    sprintf("dggs_res_spec %d", res), orient)
}

# DGGRID lines for an orientation c(vert0_lon, vert0_lat, azimuth), or for
# DGGRID's own REGION_CENTER placement about c(lon, lat).
orient_lines <- function(orientation = NULL, region = NULL) {
  if (!is.null(region)) {
    return(c("dggs_orient_specify_type REGION_CENTER",
             sprintf("region_center_lon %.15f", region[1]),
             sprintf("region_center_lat %.15f", region[2])))
  }
  c("dggs_orient_specify_type SPECIFIED",
    sprintf("dggs_vert0_lon %.15f", orientation[1]),
    sprintf("dggs_vert0_lat %.15f", orientation[2]),
    sprintf("dggs_vert0_azimuth %.15f", orientation[3]))
}

dggrid_seqnum <- function(cfg, res, lon, lat, orient = NULL) {
  inp <- file.path(work, "pts.txt"); outp <- file.path(work, "seq.txt")
  write.table(data.frame(sprintf("%.12f", lon), sprintf("%.12f", lat)), inp,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  run_dggrid(c("dggrid_operation TRANSFORM_POINTS", dggs_lines(cfg, res, orient),
               paste("input_file_name", inp), "input_address_type GEO",
               "input_delimiter \" \"", paste("output_file_name", outp),
               "output_address_type SEQNUM", "output_delimiter \" \""))
  as.numeric(readLines(outp))
}

dggrid_generate <- function(cfg, res, seqnum, cells = FALSE, orient = NULL) {
  inp <- file.path(work, "clip.txt")
  cstem <- file.path(work, "cells"); pstem <- file.path(work, "points")
  writeLines(format(unique(seqnum), scientific = FALSE, trim = TRUE), inp)
  unlink(c(paste0(cstem, ".gen"), paste0(pstem, ".txt")))
  run_dggrid(c("dggrid_operation GENERATE_GRID", dggs_lines(cfg, res, orient),
               "clip_subset_type ADDRESS_FILES", "input_address_type SEQNUM",
               paste("clip_region_files", inp),
               if (cells) c("cell_output_type AIGEN", paste("cell_output_file_name", cstem))
               else "cell_output_type NONE",
               "point_output_type TEXT", paste("point_output_file_name", pstem),
               "densification 0", "precision 12"))
  p <- read.csv(paste0(pstem, ".txt"), header = FALSE)
  ctr <- data.frame(seqnum = as.numeric(p[[1]]), lon = p[[2]], lat = p[[3]])
  list(centres = ctr, cells = if (cells) read_aigen(paste0(cstem, ".gen")))
}

# AIGEN: per cell a line "seqnum lon lat" (the centre), the corner lines
# "lon lat", and "END"; the file ends with a further "END".
read_aigen <- function(path) {
  x <- trimws(readLines(path))
  ends <- which(x == "END")
  starts <- c(1L, head(ends, -1) + 1L)
  keep <- starts < ends
  starts <- starts[keep]; ends <- ends[keep]
  ids <- as.numeric(sub("\\s.*$", "", x[starts]))
  corners <- lapply(seq_along(starts), function(k) {
    v <- x[(starts[k] + 1L):(ends[k] - 1L)]
    m <- do.call(rbind, lapply(strsplit(v, "[ ,]+"), as.numeric))
    m[, 1:2, drop = FALSE]
  })
  list(seqnum = ids, corners = corners)
}

dggrid_centres <- function(cfg, res, seqnum, orient = NULL) {
  ctr <- dggrid_generate(cfg, res, seqnum, orient = orient)$centres
  ctr[match(seqnum, ctr$seqnum), c("lon", "lat")]
}

dggrid_polygons <- function(cfg, res, seqnum, orient = NULL) {
  dggrid_generate(cfg, res, seqnum, cells = TRUE, orient = orient)$cells
}

# Coordinates in DGGRID's PLANE system (the unfolded icosahedron), for points
# given as GEO lon/lat or for cell centres given as SEQNUMs.
dggrid_plane <- function(cfg, res, lon = NULL, lat = NULL, seqnum = NULL) {
  inp <- file.path(work, "plane_in.txt"); outp <- file.path(work, "plane_out.txt")
  if (is.null(seqnum)) {
    write.table(data.frame(sprintf("%.12f", lon), sprintf("%.12f", lat)), inp,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    in_type <- "GEO"
  } else {
    writeLines(format(seqnum, scientific = FALSE, trim = TRUE), inp)
    in_type <- "SEQNUM"
  }
  run_dggrid(c("dggrid_operation TRANSFORM_POINTS", dggs_lines(cfg, res),
               paste("input_file_name", inp), paste("input_address_type", in_type),
               "input_delimiter \" \"", paste("output_file_name", outp),
               "output_address_type PLANE", "output_delimiter \" \"", "precision 15"))
  m <- read.table(outp)
  data.frame(x = m[[1]], y = m[[2]])
}
