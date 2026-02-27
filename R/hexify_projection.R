# hexify_projection.R
# ISEA Snyder Equal-Area Projection
#
# This module handles all projection operations between geographic coordinates
# (longitude/latitude) and the planar icosahedral faces.
#
# Note: The C++ layer handles lazy initialization of the icosahedron via ico().
# No R-level state management is needed.
#
# @name hexify-projection

# =============================================================================
# ICOSAHEDRON INITIALIZATION
# =============================================================================

#' Initialize icosahedron geometry
#'
#' Sets up the icosahedron state for ISEA projection. Uses standard ISEA3H
#' orientation by default (vertex 0 at 11.25E, 58.28N).
#'
#' @param vert0_lon Vertex 0 longitude in degrees (default ISEA_VERT0_LON_DEG)
#' @param vert0_lat Vertex 0 latitude in degrees (default ISEA_VERT0_LAT_DEG)
#' @param azimuth Azimuth rotation in degrees (default ISEA_AZIMUTH_DEG)
#'
#' @return Invisible NULL. Called for side effect.
#'
#' @details
#' The icosahedron is initialized lazily at the C++ level when first needed.
#' Manual call is only required for non-standard orientations.
#'
#' @family projection
#' @export
#' @examples
#' # Use standard ISEA3H orientation
#' hexify_build_icosa()
#'
#' # Custom orientation
#' hexify_build_icosa(vert0_lon = 0, vert0_lat = 90, azimuth = 0)
hexify_build_icosa <- function(vert0_lon = ISEA_VERT0_LON_DEG,
                                vert0_lat = ISEA_VERT0_LAT_DEG,
                                azimuth = ISEA_AZIMUTH_DEG) {
  cpp_build_icosa(
    as.numeric(vert0_lon),
    as.numeric(vert0_lat),
    as.numeric(azimuth)
  )
  invisible(NULL)
}

#' Get icosahedron face centers
#'
#' Returns the center coordinates of all 20 icosahedral faces.
#'
#' @return Data frame with 20 rows and columns lon, lat (degrees)
#'
#' @family projection
#' @export
#' @examples
#' centers <- hexify_face_centers()
#' plot(centers$lon, centers$lat)
hexify_face_centers <- function() {
  cpp_face_centers()
}

# =============================================================================
# FORWARD PROJECTION (lon/lat -> face coordinates)
# =============================================================================

#' Determine which face contains a point
#'
#' Returns the icosahedral face index (0-19) containing the given coordinates.
#'
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#'
#' @return Integer face index (0-19)
#'
#' @family projection
#' @export
#' @examples
#' face <- hexify_which_face(16.37, 48.21)
hexify_which_face <- function(lon, lat) {
  cpp_which_face(as.numeric(lon), as.numeric(lat))
}

#' Forward Snyder projection
#'
#' Projects geographic coordinates onto the icosahedron, returning
#' face index and planar coordinates (tx, ty).
#'
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#'
#' @return Named numeric vector: c(face, tx, ty)
#'
#' @details
#' tx and ty are normalized coordinates within the triangular face,
#' typically in range \[0, 1\].
#'
#' @family projection
#' @export
#' @examples
#' result <- hexify_forward(16.37, 48.21)
#' # result["face"], result["icosa_triangle_x"], result["icosa_triangle_y"]
hexify_forward <- function(lon, lat) {
  cpp_snyder_forward(as.numeric(lon), as.numeric(lat))
}

#' Forward projection to specific face
#'
#' Projects to a known face (skips face detection).
#'
#' @param face Face index (0-19)
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#'
#' @return Named numeric vector: c(icosa_triangle_x, icosa_triangle_y)
#'
#' @family projection
#' @export
hexify_forward_to_face <- function(face, lon, lat) {
  cpp_project_to_icosa_triangle(as.integer(face), as.numeric(lon), as.numeric(lat))
}

# =============================================================================
# INVERSE PROJECTION (face coordinates -> lon/lat)
# =============================================================================

#' Inverse Snyder projection
#'
#' Converts face plane coordinates back to geographic coordinates.
#'
#' @param x X coordinate on face plane
#' @param y Y coordinate on face plane
#' @param face Face index (0-19)
#' @param tol Convergence tolerance (NULL for default)
#' @param max_iters Maximum iterations (NULL for default)
#'
#' @return Named numeric vector: c(lon_deg, lat_deg)
#'
#' @family projection
#' @export
#' @examples
#' coords <- hexify_inverse(0.5, 0.3, face = 2)
hexify_inverse <- function(x, y, face, tol = NULL, max_iters = NULL) {
  stopifnot(length(x) == 1L, length(y) == 1L, length(face) == 1L)
  cpp_face_xy_to_ll(as.numeric(x), as.numeric(y), as.integer(face), tol, max_iters)
}

# =============================================================================
# PRECISION CONTROL
# =============================================================================

#' Set inverse projection precision
#'
#' Controls the accuracy/speed tradeoff for inverse Snyder projection.
#'
#' @param mode Preset mode: "fast", "default", "high", or "ultra"
#' @param tol Custom tolerance (overrides mode if provided)
#' @param max_iters Custom max iterations (overrides mode if provided)
#'
#' @return Invisible NULL
#'
#' @family projection
#' @export
#' @examples
#' hexify_set_precision("high")
#' hexify_set_precision(tol = 1e-12, max_iters = 100)
hexify_set_precision <- function(mode = c("fast", "default", "high", "ultra"),
                                  tol = NULL, max_iters = NULL) {
  mode <- match.arg(mode)
  cpp_snyder_inv_set_precision(mode, tol, max_iters)
  invisible(NULL)
}

#' Get current precision settings
#'
#' @return List with tol and max_iters
#'
#' @family projection
#' @export
hexify_get_precision <- function() {
  cpp_snyder_inv_get_precision()
}

#' Set verbose mode for inverse projection
#'
#' When enabled, prints convergence information.
#'
#' @param verbose Logical
#'
#' @return Invisible NULL
#'
#' @family projection
#' @export
hexify_set_verbose <- function(verbose = TRUE) {
  cpp_snyder_inv_set_verbose(isTRUE(verbose))
  invisible(NULL)
}

#' Get inverse projection statistics
#'
#' Returns and optionally resets convergence statistics.
#'
#' @param reset Whether to reset statistics after retrieval (default TRUE)
#'
#' @return List with statistics (iterations, convergence info, etc.)
#'
#' @family projection
#' @export
hexify_projection_stats <- function(reset = TRUE) {
  cpp_snyder_inv_get_stats_and_reset()
}

