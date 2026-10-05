# Unit vectors and angles on the sphere, shared by the test files.

unit_rows <- function(lon, lat) {
  r <- pi / 180
  cbind(cos(lat * r) * cos(lon * r), cos(lat * r) * sin(lon * r), sin(lat * r))
}

# Angle between matching rows of two sets of unit vectors, in degrees
row_angle <- function(A, B) {
  cr <- cbind(A[, 2] * B[, 3] - A[, 3] * B[, 2], A[, 3] * B[, 1] - A[, 1] * B[, 3],
              A[, 1] * B[, 2] - A[, 2] * B[, 1])
  atan2(sqrt(rowSums(cr^2)), rowSums(A * B)) * 180 / pi
}

rows_lonlat <- function(P) {
  cbind(atan2(P[, 2], P[, 1]), asin(pmax(-1, pmin(1, P[, 3])))) * 180 / pi
}
