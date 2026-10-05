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
  list("7", 7, c("dggs_aperture_type PURE", "dggs_aperture 7"), 2:8),
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

# `icosa` is further DGGRID lines placing the icosahedron (see
# orient_lines()) and projecting its faces (see proj_lines()); NULL keeps
# DGGRID's standard orientation, and the ISEA projection unless a dggs_proj
# line is given.
dggs_lines <- function(cfg, res, icosa = NULL) {
  ap_lines <- if (identical(cfg[[1]], "4/3")) {
    c("dggs_aperture_type MIXED43", sprintf("dggs_num_aperture_4_res %d", res %/% 2))
  } else cfg[[3]]
  proj <- if (!any(startsWith(icosa, "dggs_proj"))) proj_lines("isea")
  c("dggs_type CUSTOM", "dggs_topology HEXAGON", proj, ap_lines,
    sprintf("dggs_res_spec %d", res), icosa)
}

# The DGGRID line for a face projection as hex_grid() names it.
proj_lines <- function(projection) paste("dggs_proj", toupper(projection))

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

dggrid_seqnum <- function(cfg, res, lon, lat, icosa = NULL) {
  inp <- file.path(work, "pts.txt"); outp <- file.path(work, "seq.txt")
  write.table(data.frame(sprintf("%.12f", lon), sprintf("%.12f", lat)), inp,
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  run_dggrid(c("dggrid_operation TRANSFORM_POINTS", dggs_lines(cfg, res, icosa),
               paste("input_file_name", inp), "input_address_type GEO",
               "input_delimiter \" \"", paste("output_file_name", outp),
               "output_address_type SEQNUM", "output_delimiter \" \""))
  as.numeric(readLines(outp))
}

dggrid_generate <- function(cfg, res, seqnum, cells = FALSE, icosa = NULL) {
  inp <- file.path(work, "clip.txt")
  cstem <- file.path(work, "cells"); pstem <- file.path(work, "points")
  writeLines(format(unique(seqnum), scientific = FALSE, trim = TRUE), inp)
  unlink(c(paste0(cstem, ".gen"), paste0(pstem, ".txt")))
  run_dggrid(c("dggrid_operation GENERATE_GRID", dggs_lines(cfg, res, icosa),
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

dggrid_centres <- function(cfg, res, seqnum, icosa = NULL) {
  ctr <- dggrid_generate(cfg, res, seqnum, icosa = icosa)$centres
  ctr[match(seqnum, ctr$seqnum), c("lon", "lat")]
}

dggrid_polygons <- function(cfg, res, seqnum, icosa = NULL) {
  dggrid_generate(cfg, res, seqnum, cells = TRUE, icosa = icosa)$cells
}

# A corner further than this from every corner of the other program's cell
# counts as a corner gap; corners otherwise agree to within a millimetre.
CORNER_GAP_KM <- 1e-3

# Largest distance from a corner of `cd` to the nearest corner of `ch`, km.
corner_gap_km <- function(cd, ch) {
  max(vapply(seq_len(nrow(cd)), function(i)
    min(gc_km(cd[i, 1], cd[i, 2], ch[, 1], ch[, 2])), numeric(1)))
}

# How many of `rings` have a corner within CORNER_GAP_KM of each row of `p`.
corner_owners <- function(p, rings) {
  vapply(seq_len(nrow(p)), function(k)
    sum(vapply(rings, function(r)
      any(gc_km(p[k, 1], p[k, 2], r[, 1], r[, 2]) < CORNER_GAP_KM), logical(1))),
    integer(1))
}

# Cell corners of both programs on `cells`. Each DGGRID corner is matched to
# the nearest hexify corner of the same cell; the largest distance is the
# cell's corner gap. Cells in a tiling meet three to a corner, so for every
# cell with a gap the corners are checked against the cell's neighbours in
# each program: a DGGRID corner no neighbour carries leaves DGGRID's polygons
# open there, and hexify's corners of that cell should each meet three cells.
corner_check <- function(cfg, res, g, cells, icosa = NULL) {
  ocpp <- hexify:::icosa_arg(g)
  rings_h <- function(ids) {
    lapply(hexify:::isea_cell_rings(ids, g@resolution, g@aperture, ocpp, 0),
           function(m) m[, 1:2, drop = FALSE])
  }
  dp <- dggrid_polygons(cfg, res, cells, icosa)
  hr <- rings_h(dp$seqnum)
  gap <- vapply(seq_along(dp$seqnum), function(i)
    corner_gap_km(dp$corners[[i]], hr[[i]]), numeric(1))

  gap_idx <- which(gap > CORNER_GAP_KM)
  dggrid_open <- logical(length(gap_idx))
  hexify_closed <- logical(length(gap_idx))
  if (length(gap_idx)) {
    nb <- hexify:::grid_neighbors_isea(dp$seqnum[gap_idx], g)
    nb_ids <- unique(unlist(nb))
    dn <- dggrid_polygons(cfg, res, nb_ids, icosa)
    dn_rings <- dn$corners[match(nb_ids, dn$seqnum)]
    hn_rings <- rings_h(nb_ids)
    for (k in seq_along(gap_idx)) {
      i <- gap_idx[k]
      own <- match(nb[[k]], nb_ids)
      d <- dp$corners[[i]]
      far <- vapply(seq_len(nrow(d)), function(j)
        min(gc_km(d[j, 1], d[j, 2], hr[[i]][, 1], hr[[i]][, 2])) > CORNER_GAP_KM, logical(1))
      dggrid_open[k] <- all(corner_owners(d[far, , drop = FALSE],
                                          c(dp$corners[i], dn_rings[own])) < 3L)
      hexify_closed[k] <- all(corner_owners(hr[[i]], c(hr[i], hn_rings[own])) == 3L)
    }
  }
  explained <- seq_along(gap) %in% gap_idx[dggrid_open & hexify_closed]
  list(n_polygons = length(gap), max_corner_gap_m = max(gap) * 1000,
       n_corner_gap = length(gap_idx),
       n_gap_dggrid_open = sum(dggrid_open), n_gap_hexify_closed = sum(hexify_closed),
       max_corner_gap_other_m = max(c(0, gap[!explained])) * 1000)
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
