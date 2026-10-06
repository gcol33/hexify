# Angles on the sphere and face-plane samples, shared by the test files.

# Angle between matching rows of two sets of unit vectors, in degrees
row_angle <- function(A, B) {
  atan2(sqrt(rowSums(hexify:::cross3(A, B)^2)), rowSums(A * B)) * 180 / pi
}

# Great-circle angle between points given in degrees, in radians
arc_between <- function(lon1, lat1, lon2, lat2) {
  d2r <- pi / 180
  u <- cbind(cos(lat1 * d2r) * cos(lon1 * d2r), cos(lat1 * d2r) * sin(lon1 * d2r),
             sin(lat1 * d2r))
  v <- cbind(cos(lat2 * d2r) * cos(lon2 * d2r), cos(lat2 * d2r) * sin(lon2 * d2r),
             sin(lat2 * d2r))
  cr <- cbind(u[, 2] * v[, 3] - u[, 3] * v[, 2], u[, 3] * v[, 1] - u[, 1] * v[, 3],
              u[, 1] * v[, 2] - u[, 2] * v[, 1])
  atan2(sqrt(rowSums(cr^2)), rowSums(u * v))
}

# Face-plane points: a barycentric grid over the triangle, edges included,
# and points on the three radii from the centre to the vertices, where
# Snyder's projection has its cusps, and to the edge midpoints, which with
# them make the face's symmetry lines.
face_plane_samples <- function(n = 30) {
  vx <- c(0, 1, 0.5)
  vy <- c(0, 0, sqrt(3) / 2)
  g <- expand.grid(i = 0:n, j = 0:n)
  g <- g[g$i + g$j <= n, ]
  b1 <- g$i / n
  b2 <- g$j / n
  b0 <- 1 - b1 - b2
  t <- seq(0.01, 1, length.out = 40)
  mx <- (vx + vx[c(2, 3, 1)]) / 2
  my <- (vy + vy[c(2, 3, 1)]) / 2
  cx <- mean(vx)
  cy <- mean(vy)
  data.frame(
    x = c(b0 * vx[1] + b1 * vx[2] + b2 * vx[3], cx + outer(t, vx - cx),
          cx + outer(t, mx - cx)),
    y = c(b0 * vy[1] + b1 * vy[2] + b2 * vy[3], cy + outer(t, vy - cy),
          cy + outer(t, my - cy))
  )
}
