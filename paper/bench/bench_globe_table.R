# The globe widget's cells' table against the per-cell upload of a base commit:
# R build time, page size, bytes the grid's cells take on the graphics card
# and the time of a frame, for grids with values on every cell, values on a
# sample of cells, and no values.
#
#   Rscript paper/bench/bench_globe_table.R [base commit, default 7996043]
#
# Run from the repository root. The base commit (the per-cell upload: values
# in a storage buffer read by cell ID, or by binary search among sorted cell
# IDs) and the working tree are each installed into a temporary library, and
# each is measured in its own R process, since one session holds one hexify.
# Frames are timed in headless Chrome on the graphics adapter it picks.

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------------------------------------------------------
# One version, in its own process: Rscript bench_globe_table.R --measure <lib> <label> <out>
# ---------------------------------------------------------------------------

if (length(args) >= 1 && args[1] == "--measure") {
  .libPaths(c(args[2], .libPaths()))
  suppressPackageStartupMessages(library(hexify, lib.loc = args[2]))
  stopifnot(normalizePath(dirname(find.package("hexify")), winslash = "/") ==
              normalizePath(args[2], winslash = "/"))
  label <- args[3]

  page_mb <- function(widget) {
    dir <- tempfile()
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE))
    html <- file.path(dir, "globe.html")
    htmlwidgets::saveWidget(widget, html, selfcontained = FALSE, libdir = "lib")
    file.size(html) / 1e6
  }

  # Bytes of the grid's cells on the graphics card: the table's textures, or
  # the per-cell values and sorted keys
  gpu_mb <- function(widget) {
    g <- widget$x$grid
    b64 <- function(s) if (is.null(s) || !nzchar(s)) 0 else floor(nchar(s) * 3 / 4)
    if (!is.null(g$textures)) {
      texel <- c(r32uint = 4, rg32uint = 8)
      bytes <- sum(vapply(g$textures, function(t) prod(t$size) * texel[[t$format]], numeric(1)))
    } else {
      bytes <- b64(g$values) + b64(g$keys)
    }
    bytes / 1e6
  }

  # Time of a frame on the graphics card: 40 frames of the view encoded into
  # one submission, which is timed from submission until the card is done,
  # ten times, so that the time per frame is the card's and not the round
  # trip of one submission; the median over the ten, per frame. And the
  # adapter's description.
  frames <- function(widget, size) {
    hexify:::globe_in_chrome(widget, size, size, 1, 300, function(session, read) {
      ms <- read(paste0(
        "(async () => {",
        "  const g = document.querySelector('.hexify-globe').hexGlobe;",
        "  const target = g.device.createTexture({ size: [g.canvas.width, g.canvas.height],",
        "    format: g.format, usage: GPUTextureUsage.RENDER_ATTACHMENT });",
        "  const batch = async () => {",
        "    g.writeCamera();",
        "    const enc = g.device.createCommandEncoder();",
        "    for (let k = 0; k < 40; k++) g.encodeFrame(enc, target);",
        "    const t0 = performance.now();",
        "    g.device.queue.submit([enc.finish()]);",
        "    await g.device.queue.onSubmittedWorkDone();",
        "    return (performance.now() - t0) / 40;",
        "  };",
        "  for (let k = 0; k < 3; k++) await batch();",
        "  const t = [];",
        "  for (let k = 0; k < 10; k++) { g.state.lon += 7; t.push(await batch()); }",
        "  t.sort((a, b) => a - b);",
        "  return (t[4] + t[5]) / 2;",
        "})()"))
      adapter <- read(paste0(
        "(async () => { const a = await navigator.gpu.requestAdapter();",
        "  const i = a.info || {}; return [i.vendor, i.architecture, i.description].join(' '); })()"))
      list(ms = ms, adapter = adapter)
    })
  }

  set.seed(93)
  sample_cells <- function(grid, n) {
    ids <- unique(lonlat_to_cell(runif(2 * n, -180, 180), asin(runif(2 * n, -1, 1)) * 180 / pi,
                                 grid))
    ids[seq_len(n)]
  }
  g10 <- hex_grid(resolution = 10, aperture = 3)
  g12 <- hex_grid(resolution = 12, aperture = 3)
  g14 <- hex_grid(resolution = 14, aperture = 4)
  lat <- function(grid) cell_to_lonlat(seq_len(n_cells(grid)), grid)$lat_deg
  s5 <- sample_cells(g14, 1e5)
  s6 <- sample_cells(g14, 1e6)
  cases <- list(
    list(case = "aperture 3 resolution 10, no values", grid = g10, args = list()),
    list(case = "aperture 3 resolution 10, values on every cell", grid = g10,
         args = list(values = lat(g10))),
    list(case = "aperture 3 resolution 12, values on every cell", grid = g12,
         args = list(values = lat(g12))),
    list(case = "aperture 4 resolution 14, values on 100,000 cells", grid = g14,
         args = list(cells = s5, values = as.numeric(s5) %% 1000)),
    list(case = "aperture 4 resolution 14, values on 1,000,000 cells", grid = g14,
         args = list(cells = s6, values = as.numeric(s6) %% 1000))
  )
  if ("smooth" %in% names(formals(hex_globe))) {
    cases <- c(cases, list(list(case = "aperture 3 resolution 10, values on every cell, smooth",
                                grid = g10, args = list(values = lat(g10), smooth = TRUE))))
  }

  rows <- lapply(cases, function(cs) {
    gc()
    build <- system.time(w <- do.call(hex_globe, c(list(cs$grid, land = FALSE), cs$args)))
    f800 <- frames(w, 800)
    f1600 <- frames(w, 1600)
    data.frame(version = label, case = cs$case, cells = as.numeric(n_cells(cs$grid)),
               build_s = build[["elapsed"]], page_mb = page_mb(w), gpu_mb = gpu_mb(w),
               frame_ms_800 = f800$ms, frame_ms_1600 = f1600$ms, adapter = f800$adapter)
  })
  write.csv(do.call(rbind, rows), args[4], row.names = FALSE)
  quit(save = "no")
}

# ---------------------------------------------------------------------------
# Both versions
# ---------------------------------------------------------------------------

source("paper/bench/bench_common.R")
base <- if (length(args) >= 1) args[1] else "7996043"
rscript <- file.path(R.home("bin"), "Rscript")
r <- file.path(R.home("bin"), "R")
work <- tempfile("bench_globe_")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

install <- function(src, lib) {
  dir.create(lib)
  status <- system2(r, c("CMD", "INSTALL", "--preclean", "--no-test-load",
                         paste0("--library=", shQuote(lib)), shQuote(src)))
  if (status != 0) stop("installing ", src, " failed")
}

base_src <- file.path(work, "base")
dir.create(base_src)
tar <- file.path(work, "base.tar")
if (system2("git", c("archive", "--format=tar", "-o", shQuote(tar), base)) != 0) {
  stop("git archive of ", base, " failed")
}
utils::untar(tar, exdir = base_src)
install(base_src, file.path(work, "lib_base"))
install(".", file.path(work, "lib_head"))

parts <- c(base = file.path(work, "base.csv"), head = file.path(work, "head.csv"))
labels <- c(base = paste0("per-cell upload (", base, ")"), head = "cells' table")
for (v in names(parts)) {
  status <- system2(rscript, c(shQuote(sys_script_path()), "--measure",
                               shQuote(file.path(work, paste0("lib_", v))),
                               shQuote(labels[[v]]), shQuote(parts[[v]])))
  if (status != 0) stop("measuring ", v, " failed")
}
out <- rbind(read.csv(parts[["base"]]), read.csv(parts[["head"]]))
print(out[, setdiff(names(out), "adapter")])
write_result(out, "globe_table")
