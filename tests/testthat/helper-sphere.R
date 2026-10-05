# Angles on the sphere, shared by the test files.

# Angle between matching rows of two sets of unit vectors, in degrees
row_angle <- function(A, B) {
  atan2(sqrt(rowSums(hexify:::cross3(A, B)^2)), rowSums(A * B)) * 180 / pi
}
