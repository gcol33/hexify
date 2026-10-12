# Build the animated hex logo man/figures/logo.svg.
#
# The hex is tessellated by the cells of an ISEA3H grid as they lie on
# Snyder's icosahedron unfolded into its plane net (net_layout(g, "plane")):
# the lattice of the net's cells is continued over the hex, and the cells
# whose centres fall on land (hexify_world) are filled. Towards the bottom
# the tessellation bends away like the lower half of a globe, and ends in a
# curve above the name, which sits on plain paper. The land cells pulse
# slowly, in a wave running from west to east. Under prefers-reduced-motion
# they stay still.
#
# Run from the package root:  Rscript data-raw/logo.R
# Needs sf only here; it is not a package dependency.

suppressMessages(library(sf))
sf_use_s2(FALSE)
pkgload::load_all(quiet = TRUE)

out <- "man/figures/logo.svg"

col <- c(border = "#2B3A4A", paper = "#FFFFFF", sea = "#CDD3DA", land = "#2E6F9E",
         name = "#2E6F9E")

# Hex geometry: pointy top, circumradius R, centred in a 520 x 600 viewBox.
W <- 520; H <- 600; R <- 300; cx <- W / 2; cy <- H / 2
hex_pts <- function(r) {
  a <- pi / 2 + (0:5) * pi / 3
  cbind(cx + r * cos(a), cy - r * sin(a))
}
pts_attr <- function(m) paste(sprintf("%.1f,%.1f", m[, 1], m[, 2]), collapse = " ")

# Net placement: net coordinates are in face edges, y up. (net_x, net_y)
# lands on the hex's centre, at `edge` viewBox units per face edge.
edge <- 230
net_x <- 2.9
net_y <- 1.45
grid <- hex_grid(resolution = 4, aperture = 3)
g <- extract_grid(grid)
net <- net_layout(grid, "plane")
resolve_layout(net, g)
to_view <- function(m) cbind(cx + edge * (m[, 1] - net_x), cy - edge * (m[, 2] - net_y))

# Every cell's parts on the net, with the part's centroid in viewBox units.
ids <- grid_cells(g, NULL)
paths <- grid_surface_paths(g, ids, step = 1)
ll <- cell_to_lonlat(ids, grid)
land <- st_union(st_geometry(hexify_world))
on_land <- lengths(suppressMessages(st_intersects(
  st_as_sf(data.frame(lon = ll$lon_deg, lat = ll$lat_deg), coords = c("lon", "lat"), crs = 4326),
  land))) > 0

rows <- split(seq_len(nrow(paths)), paths[, "cell"])
parts <- do.call(rbind, lapply(names(rows), function(k) {
  p <- cell_net_parts(paths[rows[[k]], , drop = FALSE], net)
  if (is.null(p)) return(NULL)
  xy <- to_view(cbind(p$x, p$y))
  data.frame(cell = as.integer(k), land = on_land[as.integer(k)], area = p$area,
             nv = vapply(p$polygon, nrow, 1L), x = xy[, 1], y = xy[, 2],
             poly = I(lapply(p$polygon, to_view)))
}))

# The lattice: a whole hexagonal part gives the cell's circumradius and the
# direction of its first corner; centres step by sqrt(3) circumradii at 30
# degrees past the corners.
whole <- parts[parts$nv == 6 & parts$area > 0.99 * max(parts$area[parts$nv == 6]), ][1, ]
corner <- whole$poly[[1]]
c0 <- c(whole$x, whole$y)
r <- mean(sqrt((corner[, 1] - c0[1])^2 + (corner[, 2] - c0[2])^2))
phi <- atan2(corner[1, 2] - c0[2], corner[1, 1] - c0[1])
a1 <- sqrt(3) * r * c(cos(phi + pi / 6), sin(phi + pi / 6))
a2 <- sqrt(3) * r * c(cos(phi + pi / 2), sin(phi + pi / 2))
n <- ceiling(1.5 * max(W, H) / (sqrt(3) * r))
ij <- expand.grid(i = -n:n, j = -n:n)
centres <- cbind(c0[1] + ij$i * a1[1] + ij$j * a2[1], c0[2] + ij$i * a1[2] + ij$j * a2[2])

# The bend. Above bend_from the plane is drawn as it is. Below it, with t
# running from 0 there to 1 at bend_to, the plane is foreshortened as if it
# turned away round a globe by up to `turn` radians (depth sin(t turn) /
# turn, slope 1 at t = 0 so the fold has no crease), drawn in towards the
# axis, and lifted at the sides, so its lower edge becomes a curve. Cells
# with centres below bend_to are left out.
bend_from <- 200
bend_to <- 430
turn <- 1
lift <- 45
squeeze <- 0.08
warp <- function(m) {
  t <- pmin(pmax((m[, 2] - bend_from) / (bend_to - bend_from), 0), 1)
  u <- (m[, 1] - cx) / (W / 2)
  y <- ifelse(m[, 2] > bend_from,
              bend_from + (bend_to - bend_from) * sin(t * turn) / turn, m[, 2])
  cbind(cx + (m[, 1] - cx) * (1 - squeeze * t^2), y - lift * u^2 * t^2)
}
centres <- centres[centres[, 1] > -r & centres[, 1] < W + r &
                     centres[, 2] > -r & centres[, 2] < bend_to, , drop = FALSE]

# A lattice cell is land where a land cell's part has its centroid inside it.
land_xy <- as.matrix(parts[parts$land, c("x", "y")])
is_land <- vapply(seq_len(nrow(centres)), function(i) {
  any((land_xy[, 1] - centres[i, 1])^2 + (land_xy[, 2] - centres[i, 2])^2 < (0.6 * r)^2)
}, logical(1))

# A cell's outline with a point midway along each edge, so the bend curves
# its edges, then warped.
ang <- phi + seq(0, 11) * pi / 6
rad_k <- rep(c(1, sqrt(3) / 2), 6)
hex_d <- function(p, s = 1) {
  m <- warp(cbind(p[1] + s * r * rad_k * cos(ang), p[2] + s * r * rad_k * sin(ang)))
  sprintf("M%sZ", paste(sprintf("%.1f,%.1f", m[, 1], m[, 2]), collapse = "L"))
}
sea_d <- vapply(which(!is_land), function(i) hex_d(centres[i, ]), "")

# Land cells pulse in a wave from west to east: each cell's delay is its
# place across the hex, so the wave runs over the land once per period.
period <- 4
land_svg <- vapply(which(is_land), function(i) {
  sprintf('<path class="land" d="%s" style="animation-delay:-%.2fs"/>',
          hex_d(centres[i, ], 0.9), period * (1 - pmin(pmax(centres[i, 1] / W, 0), 1)))
}, "")

css <- sprintf('
.sea { fill: none; stroke: %1$s; stroke-width: 2; stroke-linejoin: round; }
.land { fill: %2$s; animation: pulse %3$gs ease-in-out infinite; }
@keyframes pulse { 0%%, 100%% { opacity: 1; } 50%% { opacity: 0.45; } }
@media (prefers-reduced-motion: reduce) { .land { animation: none; } }
', col[["sea"]], col[["land"]], period)

svg <- c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" role="img" aria-label="hexify logo">', W, H, W, H),
  "<title>hexify logo</title>",
  sprintf("<style>%s</style>", css),
  "<defs>",
  sprintf('<clipPath id="inner"><polygon points="%s"/></clipPath>', pts_attr(hex_pts(R - 12))),
  "</defs>",
  sprintf('<polygon points="%s" fill="%s"/>', pts_attr(hex_pts(R - 12)), col[["paper"]]),
  '<g clip-path="url(#inner)">',
  sprintf('<path class="sea" d="%s"/>', paste(sea_d, collapse = "")),
  land_svg,
  sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-family="Georgia, \'Times New Roman\', serif" font-style="italic" font-weight="bold" font-size="78" fill="%s" stroke="%s" stroke-width="10" stroke-linejoin="round" paint-order="stroke">hexify</text>',
          cx, 472, col[["name"]], col[["paper"]]),
  "</g>",
  sprintf('<polygon points="%s" fill="none" stroke="%s" stroke-width="18" stroke-linejoin="round"/>', pts_attr(hex_pts(R - 9)), col[["border"]]),
  "</svg>"
)

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(svg, out, useBytes = TRUE)
message(sprintf("%s: %.0f KB, cell circumradius %.1f, %d land cells", out,
                file.size(out) / 1024, r, sum(is_land)))
