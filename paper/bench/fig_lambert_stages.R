# Still frames of Lambert's construction of Snyder's projection on one face
# (#98): (a) the Lambert step, the face on the plane touching the sphere at
# its centre, an equal-area triangle with curved edges; (b) Snyder's step,
# the plane triangle, edges straight and Tissot's ellipses stretched.
#
# Both panels show face 0 of the aperture-3 resolution-3 ISEA grid (the grid
# of the paper's icosahedron figure) at one scale, in units of the sphere's
# radius, with Tissot's ellipses of circles 0.025 radii across at the points
# of a lattice 7 steps to a face edge. The images carry no text: the paper
# sets the panel letters. Each is a square of half the text width less the
# gap, (5.38 in - 0.25 in) / 2, at 600 dpi.
#
# Writes results/fig-lambert-curved.png and results/fig-lambert-straight.png.

suppressPackageStartupMessages(library(hexify))

out_dir <- file.path("paper", "bench", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- hex_grid(resolution = 3, aperture = 3)
face <- 0L
stages <- c(curved = "lambert", straight = "snyder")
cell_col <- "#0072B2"
tissot_col <- "#CC79A7"
side_in <- (5.38 - 0.25) / 2
dpi <- 600

outline <- lapply(stages, function(s) hexify:::stage_outline(grid, face, s))
box <- apply(do.call(rbind, outline), 2, range, na.rm = TRUE)
centre <- colMeans(box)
half <- 0.53 * max(box[2, ] - box[1, ])

for (k in names(stages)) {
  s <- stages[[k]]
  file <- file.path(out_dir, paste0("fig-lambert-", k, ".png"))
  png(file, width = side_in * dpi, height = side_in * dpi, res = dpi, bg = "white")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(centre[1] + c(-half, half), centre[2] + c(-half, half), asp = 1,
              xaxs = "i", yaxs = "i")
  lines(hexify:::stage_walls(grid, face, s), col = cell_col, lwd = 1)
  lines(hexify:::stage_ellipses(grid, face, s), col = tissot_col, lwd = 1)
  polygon(outline[[k]], border = "black", lwd = 1.6)
  dev.off()
  message("wrote ", file)
}
