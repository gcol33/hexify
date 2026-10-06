# hexify_projection.R
# Snyder's equal-area projection, and Fuller's, of a solid's faces
#
# This module handles all projection operations between geographic coordinates
# (longitude/latitude) and the planar faces of the icosahedron, the octahedron
# or the tetrahedron (Snyder 1992 defines the projection on each).
#
# The C++ layer builds each solid in each orientation once, when first read.
# Functions that take no grid read the solid's default orientation; the
# icosahedron's is the one hexify_build_icosa() sets.
#
# @name hexify-projection

# =============================================================================
# ICOSAHEDRON INITIALIZATION
# =============================================================================

#' Set the default icosahedron orientation
#'
#' Sets the orientation read by the functions that take no grid: the
#' projection functions in this family and the low-level conversions that take
#' a resolution and an aperture. The standard ISEA orientation (vertex 0 at
#' 11.25E, 58.28N, azimuth 0) is the default until this is called. A grid
#' carries its own orientation (see the \code{orientation} argument of
#' \code{\link{hex_grid}}), which this does not change.
#'
#' @param vert0_lon Vertex 0 longitude in degrees (default ISEA_VERT0_LON_DEG)
#' @param vert0_lat Vertex 0 latitude in degrees (default ISEA_VERT0_LAT_DEG)
#' @param azimuth Azimuth rotation in degrees (default ISEA_AZIMUTH_DEG)
#'
#' @return Invisible NULL. Called for side effect.
#'
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

#' Get face centers of a solid
#'
#' Returns the center coordinates of every face of the solid: 20 on the
#' icosahedron, 8 on the octahedron, 4 on the tetrahedron.
#'
#' @param polyhedron The solid: "icosahedron" (default), "octahedron" or
#'   "tetrahedron", in its default orientation
#'
#' @return Data frame with one row per face and columns lon, lat (radians)
#'
#' @family projection
#' @export
#' @examples
#' centers <- hexify_face_centers()
#' plot(centers$lon, centers$lat)
#' nrow(hexify_face_centers("octahedron"))
hexify_face_centers <- function(polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {
  cpp_face_centers(projection_icosa("isea", match.arg(polyhedron)))
}

# =============================================================================
# FORWARD PROJECTION (lon/lat -> face coordinates)
# =============================================================================

#' Determine which face contains a point
#'
#' Returns the index of the face of the solid containing the given
#' coordinates (0-19 on the icosahedron).
#'
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#' @inheritParams hexify_face_centers
#'
#' @return Integer face index
#'
#' @family projection
#' @export
#' @examples
#' face <- hexify_which_face(16.37, 48.21)
#' hexify_which_face(16.37, 48.21, polyhedron = "tetrahedron")
hexify_which_face <- function(lon, lat, polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {
  cpp_which_face(projection_icosa("isea", match.arg(polyhedron)),
                 as.numeric(lon), as.numeric(lat))
}

#' Forward face projection
#'
#' Projects geographic coordinates onto a face of the solid, returning
#' face index and planar coordinates (tx, ty).
#'
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#' @param projection Face projection: \code{"isea"} (Snyder's equal-area
#'   projection) or \code{"fuller"} (Fuller's projection, defined on the
#'   icosahedron only)
#' @param polyhedron The solid: "icosahedron" (default), "octahedron" or
#'   "tetrahedron", in its default orientation. Snyder (1992) gives his
#'   equal-area projection for each; its angular distortion grows with the
#'   faces' size.
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
#' hexify_forward(16.37, 48.21, polyhedron = "tetrahedron")
hexify_forward <- function(lon, lat, projection = c("isea", "fuller"),
                           polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {
  cpp_icosa_forward(projection_icosa(projection, polyhedron), as.numeric(lon),
                    as.numeric(lat))
}

#' Forward projection to specific face
#'
#' Projects to a known face (skips face detection).
#'
#' @param face Face index (0-19 on the icosahedron)
#' @param lon Longitude in degrees
#' @param lat Latitude in degrees
#' @inheritParams hexify_forward
#'
#' @return Named numeric vector: c(icosa_triangle_x, icosa_triangle_y)
#'
#' @family projection
#' @export
hexify_forward_to_face <- function(face, lon, lat, projection = c("isea", "fuller"),
                                   polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {
  cpp_project_to_icosa_triangle(projection_icosa(projection, polyhedron),
                                as.integer(face), as.numeric(lon), as.numeric(lat))
}

# =============================================================================
# INVERSE PROJECTION (face coordinates -> lon/lat)
# =============================================================================

#' Inverse face projection
#'
#' Converts face plane coordinates back to geographic coordinates.
#'
#' @param x X coordinate on face plane
#' @param y Y coordinate on face plane
#' @param face Face index (0-19 on the icosahedron)
#' @param tol Convergence tolerance (NULL for default)
#' @param max_iters Maximum iterations (NULL for default)
#' @inheritParams hexify_forward
#'
#' @return Named numeric vector: c(lon_deg, lat_deg)
#'
#' @family projection
#' @export
#' @examples
#' coords <- hexify_inverse(0.5, 0.3, face = 2)
hexify_inverse <- function(x, y, face, tol = NULL, max_iters = NULL,
                           projection = c("isea", "fuller"),
                           polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {
  stopifnot(length(x) == 1L, length(y) == 1L, length(face) == 1L)
  cpp_face_xy_to_ll(projection_icosa(projection, polyhedron), as.numeric(x),
                    as.numeric(y), as.integer(face), tol, max_iters)
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

