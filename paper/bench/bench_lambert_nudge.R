# How far Snyder's adjustment moves the points of a face away from Lambert's
# azimuthal equal-area projection about the face centre (#98).
#
# Lambert's construction of Snyder's projection (projection_stages()) takes a
# point S of a face to the Lambert point L on the plane touching the sphere
# at the face centre, moves it within that plane to the nudged point N, and
# scales the plane about the sphere's centre by R' onto the face plane, where
# it lands on Snyder's point P. This script measures, over one face:
#   - the swing |SL| (in space), over the face edge on the sphere's chord;
#   - the nudge |LN| within the tangent plane, over the edge of the nudged
#     triangle (the plane triangle's edge over R');
#   - |LP| between the Lambert point and Snyder's point, both as plane
#     coordinates from the face centre (the planes are parallel), over the
#     plane triangle's edge, and in space;
#   - the turn |Az' - Az| of the azimuth, in degrees.
# The points are a triangular lattice on the plane triangle, 200 steps along
# an edge (20,301 points with the corners and edges); Snyder's projection is
# equal-area, so they are spread evenly over the face's area and the plain
# mean is the mean over the face's area.

source("paper/bench/bench_common.R")

nudge_stats <- function(polyhedron, steps = 200L) {
  g <- hex_grid(resolution = 1, aperture = 4, polyhedron = polyhedron)
  ic <- hexify:::icosa_arg(g)
  ij <- expand.grid(i = 0:steps, j = 0:steps)
  ij <- ij[ij$i + ij$j <= steps, ]
  tx <- (ij$i + 0.5 * ij$j) / steps
  ty <- ij$j * sqrt(3) / 2 / steps
  ll <- t(vapply(seq_along(tx), function(k) {
    hexify:::cpp_face_xy_to_ll(ic, tx[k], ty[k], 0L)
  }, numeric(2)))
  s <- projection_stages(g, ll[, 1], ll[, 2], face = 0L)
  r1 <- attr(s, "r1")
  n_faces <- nrow(hexify:::icosa_solid(ic)$faces)
  edge <- sqrt(16 * pi / (sqrt(3) * n_faces))
  info <- hexify:::solid_info(polyhedron)
  chord <- 2 * sin(info$edge_arc_deg * pi / 360)
  stopifnot(max(abs(s$tx - tx)) < 1e-9, max(abs(s$ty - ty)) < 1e-9)

  z <- s$arc * pi / 180
  swing <- sqrt((s$lambert_u - s$sphere_u)^2 + (s$lambert_v - s$sphere_v)^2 +
                  (1 - cos(z))^2)
  nudge <- sqrt((s$nudge_u - s$lambert_u)^2 + (s$nudge_v - s$lambert_v)^2)
  lp_plane <- sqrt((s$snyder_u - s$lambert_u)^2 + (s$snyder_v - s$lambert_v)^2)
  lp_space <- sqrt(lp_plane^2 + (1 - r1)^2)
  turn <- abs((s$az_prime - s$az + 180) %% 360 - 180)
  row <- function(quantity, d, per, unit) {
    k <- which.max(d)
    data.frame(polyhedron = polyhedron, quantity = quantity, per = unit,
               max = max(d) / per, mean = mean(d) / per,
               max_at_arc_deg = s$arc[k], max_at_az_deg = s$az[k],
               max_at_tx = s$tx[k], max_at_ty = s$ty[k], n_points = nrow(s))
  }
  rbind(
    row("swing |SL|", swing, chord, "face edge chord"),
    row("nudge |LN|", nudge, edge / r1, "nudged triangle edge"),
    row("|LP| in the plane", lp_plane, edge, "plane triangle edge"),
    row("|LP| in space", lp_space, edge, "plane triangle edge"),
    row("turn |Az' - Az|", turn, 1, "degrees")
  )
}

out <- rbind(nudge_stats("icosahedron"), nudge_stats("octahedron"))
print(out, digits = 4)
write_result(out, "lambert_nudge")
