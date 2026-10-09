# The truncated icosahedron read through the icosahedral grid's face frame.
#
# Snyder's (1992) equal-area projection of the truncated icosahedron cuts each
# spherical icosahedron face into the spherical hexagon (twelve right
# triangles about the face centre C) and a fifth of each of the three pentagons
# (two right triangles about the vertex V each), and maps every right triangle
# onto a plane triangle by a slice from its polygon centre. A slice map
# followed by an affine map of its plane triangle is again a slice map, so the
# only way to carry this projection into the plane triangle of the icosahedral
# cell index, keeping its pieces straight-edged, its three-fold symmetry and
# its area exactly, is to slice the same eighteen spherical triangles onto
# eighteen plane triangles of the face frame with areas in proportion. Their
# corners are then fixed: the truncated icosahedron's vertices sit at t of the
# frame edge from each corner, and the midpoint of the hexagon-pentagon edge on
# the frame's bisector at s from the corner.
#
# For that composite this script measures the area shares, t and s, the
# offset between the two sides of each hexagon-pentagon edge, and the
# angular deformation over a face, beside Snyder's ISEA, IVEA and Fuller's
# projection from hexify, and beside the lower bound 2 asin(1/11) that holds
# near each vertex for every equal-area projection onto the icosahedron's
# faces (vignettes/theory.Rmd, "The truncated icosahedron").
#
# Usage: Rscript paper/bench/bench_trunc_icosa.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

EARTH_KM <- 6371.007180918475
N_FACE <- 100000L      # points on one face for the composite
N_GLOBE <- 2000000L    # points on the sphere for hexify's projections (100,000 per face)
CAP_DEG <- 0.5         # radius of the cap about a vertex
N_CAP <- 20000L

# ---- spherical geometry ------------------------------------------------------
nrm <- function(v) v / sqrt(sum(v^2))
cross <- function(a, b) c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                          a[1] * b[2] - a[2] * b[1])
arc <- function(a, b) atan2(sqrt(sum(cross(a, b)^2)), sum(a * b))
# Area of the spherical triangle abc on the unit sphere (Eriksson's formula)
excess <- function(a, b, c) {
  2 * atan2(abs(sum(a * cross(b, c))), 1 + sum(a * b) + sum(b * c) + sum(c * a))
}
unit_xyz <- function(lon, lat) {
  r <- pi / 180
  cbind(cos(lat * r) * cos(lon * r), cos(lat * r) * sin(lon * r), sin(lat * r))
}

# One icosahedron face: a vertex at the pole and two at latitude atan(1/2)
V <- list(c(0, 0, 1), c(cos(atan(0.5)), 0, sin(atan(0.5))),
          c(cos(atan(0.5)) * cos(2 * pi / 5), cos(atan(0.5)) * sin(2 * pi / 5), sin(atan(0.5))))
C <- nrm(V[[1]] + V[[2]] + V[[3]])
S_FACE <- 4 * pi / 20

# Truncated icosahedron: its vertex on edge ij at a third of the edge from V_i,
# the midpoint of the hexagon-pentagon edge at corner i, the midpoint of edge ij
ti_vertex <- function(i, j) nrm(2 / 3 * V[[i]] + 1 / 3 * V[[j]])
cut_mid <- function(i) {
  o <- setdiff(1:3, i)
  nrm(ti_vertex(i, o[1]) + ti_vertex(i, o[2]))
}
edge_mid <- function(i, j) nrm(V[[i]] + V[[j]])

hex_tri <- excess(C, ti_vertex(1, 2), edge_mid(1, 2))
hex_tri_cut <- excess(C, cut_mid(1), ti_vertex(1, 2))
pent_tri <- excess(V[[1]], ti_vertex(1, 2), cut_mid(1))
stopifnot(abs(hex_tri - hex_tri_cut) < 1e-14,
          abs(12 * hex_tri + 6 * pent_tri - S_FACE) < 1e-14)

# ---- the face frame ----------------------------------------------------------
# Unit-edge plane triangle; its area stands for the face's 4 pi / 20.
Vp <- list(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2))
Cp <- c(0.5, sqrt(3) / 6)
T_AREA <- sqrt(3) / 4
# Triangle C' A' M'_edge has area (1/2 - t) sqrt(3)/12; triangle V' A' M'_cut
# has area t s / 4 (angle 30 degrees at V').
t_ti <- 0.5 - (hex_tri / S_FACE * T_AREA) / (sqrt(3) / 12)
s_ti <- 4 * (pent_tri / S_FACE * T_AREA) / t_ti
ti_vertex_p <- function(i, j) Vp[[i]] + t_ti * (Vp[[j]] - Vp[[i]])
cut_mid_p <- function(i) Vp[[i]] + s_ti * (Cp - Vp[[i]]) / sqrt(sum((Cp - Vp[[i]])^2))
edge_mid_p <- function(i, j) (Vp[[i]] + Vp[[j]]) / 2

pieces <- list()
add_piece <- function(P, Q, R, Pp, Qp, Rp, kind) {
  pieces[[length(pieces) + 1L]] <<- list(P = P, Q = Q, R = R, Pp = Pp, Qp = Qp, Rp = Rp,
                                         kind = kind)
}
for (i in 1:3) for (j in setdiff(1:3, i)) {
  add_piece(V[[i]], ti_vertex(i, j), cut_mid(i), Vp[[i]], ti_vertex_p(i, j), cut_mid_p(i),
            "pentagon")
  add_piece(C, cut_mid(i), ti_vertex(i, j), Cp, cut_mid_p(i), ti_vertex_p(i, j), "hexagon")
  add_piece(C, ti_vertex(i, j), edge_mid(i, j), Cp, ti_vertex_p(i, j), edge_mid_p(i, j),
            "hexagon")
}
plane_area <- function(a, b, c) abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
for (p in pieces) {
  stopifnot(abs(plane_area(p$Pp, p$Qp, p$Rp) / T_AREA - excess(p$P, p$Q, p$R) / S_FACE) < 1e-14)
}

in_piece <- function(X, p) {
  s <- c(sum(cross(p$P, p$Q) * X), sum(cross(p$Q, p$R) * X), sum(cross(p$R, p$P) * X))
  all(s >= 0) || all(s <= 0)
}
which_piece <- function(X) {
  for (k in seq_along(pieces)) if (in_piece(X, pieces[[k]])) return(k)
  NA_integer_
}

# Slice of the spherical triangle PQR onto P'Q'R' from the apex P: the arc from
# P through X meets QR at D, which goes to the point of Q'R' that cuts off the
# same share of area; X goes to the point of P'D' at sin(PX/2) / sin(PD/2).
slice <- function(X, p) {
  D <- nrm(cross(cross(p$P, X), cross(p$Q, p$R)))
  if (sum(D * (p$Q + p$R)) < 0) D <- -D
  u <- excess(p$P, p$Q, D) / excess(p$P, p$Q, p$R)
  Dp <- p$Qp + u * (p$Rp - p$Qp)
  p$Pp + sin(arc(p$P, X) / 2) / sin(arc(p$P, D) / 2) * (Dp - p$Pp)
}

# Tissot scale factors of the slice holding X, by central differences in a
# tangent frame at X, with face-frame lengths in units of the unit sphere
L_FRAME <- sqrt(S_FACE / T_AREA)
tissot <- function(X, k, h) {
  e1 <- nrm(cross(c(0.3, 0.7, 0.2), X))
  e2 <- cross(X, e1)
  f <- function(Y) slice(nrm(Y), pieces[[k]])
  J <- L_FRAME * cbind(f(X + h * e1) - f(X - h * e1), f(X + h * e2) - f(X - h * e2)) / (2 * h)
  d <- svd(J)$d
  c(a = d[1], b = d[2])
}
omega_deg <- function(a, b) 2 * asin((a - b) / (a + b)) * 180 / pi

in_face <- function(Z) {
  fn <- cbind(cross(V[[1]], V[[2]]), cross(V[[2]], V[[3]]), cross(V[[3]], V[[1]]))
  rowSums((Z %*% fn) >= 0) == 3
}

# ---- offset along the hexagon-pentagon edges ---------------------------------
# A point D of the arc M_cut A goes to M'_cut + u_hex (A' - M'_cut) from the
# hexagon's slice and to A' + u_pent (M'_cut - A') from the pentagon's.
M_s <- cut_mid(1)
A_s <- ti_vertex(1, 2)
half_edge <- arc(M_s, A_s)
q <- seq(0, 1, length.out = 20001)[-c(1, 20001)]
offset <- vapply(q, function(qq) {
  D <- nrm(sin((1 - qq) * half_edge) * M_s + sin(qq * half_edge) * A_s)
  excess(C, M_s, D) / excess(C, M_s, A_s) - (1 - excess(V[[1]], A_s, D) / excess(V[[1]], A_s, M_s))
}, numeric(1))
offset_frac <- max(abs(offset))
seam_len_frame <- sqrt(sum((ti_vertex_p(1, 2) - cut_mid_p(1))^2))
offset_m <- offset_frac * seam_len_frame * L_FRAME * EARTH_KM * 1000

# ---- angular deformation of the composite ------------------------------------
set.seed(95)
Z <- NULL
while (is.null(Z) || nrow(Z) < N_FACE) {
  W <- matrix(rnorm(3 * 400000), ncol = 3)
  W <- W / sqrt(rowSums(W^2))
  Z <- rbind(Z, W[in_face(W), , drop = FALSE])
}
Z <- Z[seq_len(N_FACE), ]
ti <- t(apply(Z, 1, function(X) tissot(X, which_piece(X), 1e-6)))
ti_omega <- omega_deg(ti[, "a"], ti[, "b"])

cap_points <- function(centre, n, radius_deg) {
  w1 <- nrm(cross(c(0.3, 0.7, 0.2), centre))
  w2 <- cross(centre, w1)
  d <- acos(1 - runif(n) * (1 - cos(radius_deg * pi / 180)))
  az <- runif(n, 0, 2 * pi)
  cos(d) * matrix(centre, n, 3, byrow = TRUE) +
    sin(d) * (cos(az) %o% w1 + sin(az) %o% w2)
}
Zc <- cap_points(V[[1]], 5L * N_CAP, CAP_DEG)
Zc <- Zc[in_face(Zc), , drop = FALSE]
tic <- t(apply(Zc, 1, function(X) tissot(X, which_piece(X), 1e-8)))
ti_cap <- omega_deg(tic[, "a"], tic[, "b"])

summ <- function(name, equal_area, omega, areal, cap, seam_m) {
  data.frame(projection = name, equal_area = equal_area, continuous = seam_m == 0,
             omega_mean = mean(omega), omega_median = median(omega),
             omega_p95 = unname(quantile(omega, 0.95)), omega_max = max(omega),
             omega_min_near_vertex = min(cap), omega_max_near_vertex = max(cap),
             areal_min = min(areal), areal_max = max(areal), seam_offset_m = seam_m)
}
rows <- list(summ("truncated icosahedron in the face frame", TRUE, ti_omega,
                  ti[, "a"] * ti[, "b"], ti_cap, offset_m))

# ---- hexify's projections ----------------------------------------------------
pts <- sphere_points(N_GLOBE, 95)
for (pr in c("isea", "ivea", "fuller")) {
  g <- hex_grid(resolution = 3, aperture = 3, projection = pr)
  d <- projection_distortion(g, pts$lon, pts$lat)
  v0 <- unit_xyz(g@orientation[1], g@orientation[2])[1, ]
  cp <- cap_points(v0, N_CAP, CAP_DEG)
  dc <- projection_distortion(g, atan2(cp[, 2], cp[, 1]) * 180 / pi, asin(cp[, 3]) * 180 / pi)
  rows[[length(rows) + 1L]] <- summ(pr, pr != "fuller", d$angular, d$areal, dc$angular, 0)
}
out <- do.call(rbind, rows)
out$omega_bound_near_vertex <- ifelse(out$equal_area, 2 * asin(1 / 11) * 180 / pi, NA_real_)
print(out, digits = 6)

write_result(out, "trunc_icosa", extra = c(
  sprintf("spherical hexagon share of a face: %.10f (plane hexagon cut at 1/3: %.10f)",
          12 * hex_tri / S_FACE, 2 / 3),
  sprintf("spherical pentagon fifth share of a face: %.10f (plane corner cut at 1/3: %.10f)",
          2 * pent_tri / S_FACE, 1 / 9),
  sprintf("frame corners: t = %.10f of the edge, s = %.10f of the edge along the bisector",
          t_ti, s_ti),
  sprintf("hexagon-pentagon offset: %.6e of the half-edge, at %.4f of it from the midpoint",
          offset_frac, q[which.max(abs(offset))]),
  sprintf("points: %d on one face (composite), %d on the sphere (hexify), caps of %.1f deg",
          N_FACE, N_GLOBE, CAP_DEG)))
