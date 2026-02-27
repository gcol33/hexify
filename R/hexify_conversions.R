# hexify_conversions.R
# Coordinate system conversion functions
#
# This file contains functions for converting between geographic coordinates
# (longitude/latitude) and hexagonal cell indices.

#' @title Coordinate Conversions
#' @description Functions for converting between coordinate systems
#' @name hexify-conversions
NULL

#' Convert longitude/latitude to hexagonal cell hierarchical index
#'
#' Converts geographic coordinates (longitude, latitude) to hexagonal cell
#' hierarchical index strings. These strings encode the face, resolution, and
#' cell location in a Z-order (Morton code) format.
#'
#' @param grid Grid specification from hexify_grid()
#' @param lon Longitude vector in degrees (numeric, -180 to 180)
#' @param lat Latitude vector in degrees (numeric, -90 to 90)
#'
#' @return Data frame with columns:
#'   \item{h_index}{Hierarchical index (character string)}
#'   \item{face}{Icosahedron face number (integer, 0-19)}
#'
#' @details
#' Most users should use \code{\link{hexify_lonlat_to_cell}} or
#' \code{\link{hexify_grid_to_cell}} which return DGGRID-compatible
#' integer cell IDs.
#'
#' This function returns hierarchical index strings useful for:
#' - Understanding the cell's position in the hierarchy
#' - Prefix-based spatial queries
#' - Parent/child cell relationships
#'
#' @keywords internal
hexify_lonlat_to_h_index <- function(grid, lon, lat) {
  
  # Validate grid object
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }
  
  # Validate inputs
  if (length(lon) != length(lat)) {
    stop("lon and lat must have the same length")
  }
  
  if (!is.numeric(lon) || !is.numeric(lat)) {
    stop("lon and lat must be numeric vectors")
  }
  
  # Check for out-of-range coordinates
  if (any(!is.na(lon) & (lon < -180 | lon > 180))) {
    warning("Some longitude values are outside valid range [-180, 180]")
  }
  
  if (any(!is.na(lat) & (lat < -90 | lat > 90))) {
    warning("Some latitude values are outside valid range [-90, 90]")
  }
  
  n <- length(lon)
  
  # Preallocate result vectors
  cell_indices <- character(n)
  faces <- integer(n)
  
  # Process each point
  for (i in seq_along(lon)) {
    # Handle NA values
    if (is.na(lon[i]) || is.na(lat[i])) {
      cell_indices[i] <- NA_character_
      faces[i] <- NA_integer_
      next
    }
    
    # Call unified C++ conversion function
    # This function:
    # 1. Projects lon/lat to icosahedron face
    # 2. Converts to face plane coordinates
    # 3. Finds hex cell containing the point
    # 4. Encodes cell as index string
    cell_indices[i] <- cpp_lonlat_to_index(
      lon[i], lat[i], 
      grid$resolution,
      grid$aperture,
      grid$index_type
    )
    
    # Extract face number from index (first 2 characters)
    # Index format: "FF..." where FF is 2-digit face number
    faces[i] <- as.integer(substr(cell_indices[i], 1, 2))
  }
  
  return(data.frame(
    h_index = cell_indices,
    face = faces
  ))
}

#' Convert hierarchical index strings to longitude/latitude centers
#'
#' Converts hierarchical cell index strings back to geographic coordinates,
#' returning the center point of each cell. This is the inverse operation
#' of hexify_lonlat_to_h_index().
#'
#' @param grid Grid specification from hexify_grid()
#' @param h_index Hierarchical index strings (character vector)
#'
#' @return Data frame with columns:
#'   \item{lon}{Longitude in degrees}
#'   \item{lat}{Latitude in degrees}
#'
#' @details
#' Most users should use \code{\link{hexify_cell_to_lonlat}} or
#' \code{\link{hexify_grid_cell_to_lonlat}} which work with DGGRID-compatible
#' integer cell IDs.
#'
#' @keywords internal
hexify_h_index_to_lonlat <- function(grid, h_index) {

  # Validate grid object
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }

  # Validate input
  if (!is.character(h_index)) {
    stop("h_index must be a character vector")
  }

  n <- length(h_index)

  # Preallocate result vectors
  lon <- numeric(n)
  lat <- numeric(n)
  
  # Process each cell index
  for (i in seq_along(h_index)) {
    # Handle NA values
    if (is.na(h_index[i])) {
      lon[i] <- NA_real_
      lat[i] <- NA_real_
      next
    }

    # Call C++ conversion function
    # This function:
    # 1. Decodes the index string to face, i, j, resolution
    # 2. Calculates hex center in face plane coordinates
    # 3. Projects back to lon/lat using inverse ISEA projection
    result <- cpp_index_to_lonlat(
      h_index[i],
      grid$aperture,
      grid$index_type
    )

    # Extract lon/lat from named vector
    lon[i] <- result["lon"]
    lat[i] <- result["lat"]
  }

  return(data.frame(
    lon = lon,
    lat = lat
  ))
}

# =============================================================================
# GRID-BASED WRAPPERS FOR CELL IDs (DGGRID-COMPATIBLE)
# =============================================================================

#' Convert longitude/latitude to cell ID using a grid object
#'
#' Grid-based wrapper for \code{\link{hexify_lonlat_to_cell}}. Converts
#' geographic coordinates to DGGRID-compatible cell IDs using
#' the resolution and aperture from a grid object.
#'
#' @param grid Grid specification from hexify_grid()
#' @param lon Numeric vector of longitudes in degrees
#' @param lat Numeric vector of latitudes in degrees
#'
#' @return Numeric vector of cell IDs (1-based)
#'
#' @family coordinate conversion
#' @seealso \code{\link{lonlat_to_cell}} for the recommended S4 interface,
#'   \code{\link{hexify_grid_cell_to_lonlat}} for the inverse operation
#'
#' @keywords internal
#' @export
#' @examples
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' cell_ids <- hexify_grid_to_cell(grid, lon = c(0, 10), lat = c(45, 50))
hexify_grid_to_cell <- function(grid, lon, lat) {
  # Validate grid object

  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }

  hexify_lonlat_to_cell(lon, lat, grid$resolution, grid$aperture)
}

#' Convert cell ID to longitude/latitude using a grid object
#'
#' Grid-based wrapper for \code{\link{hexify_cell_to_lonlat}}. Converts
#' DGGRID-compatible cell IDs back to cell center coordinates
#' using the resolution and aperture from a grid object.
#'
#' @param grid Grid specification from hexify_grid()
#' @param cell_id Numeric vector of cell IDs (1-based)
#'
#' @return Data frame with lon_deg and lat_deg columns
#'
#' @family coordinate conversion
#' @seealso \code{\link{cell_to_lonlat}} for the recommended S4 interface,
#'   \code{\link{hexify_grid_to_cell}} for the forward operation
#'
#' @keywords internal
#' @export
#' @examples
#' grid <- hexify_grid(area = 1000, aperture = 3)
#' cell_ids <- hexify_grid_to_cell(grid, lon = 5, lat = 45)
#' coords <- hexify_grid_cell_to_lonlat(grid, cell_ids)
hexify_grid_cell_to_lonlat <- function(grid, cell_id) {
  # Validate grid object
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }

  hexify_cell_to_lonlat(cell_id, grid$resolution, grid$aperture)
}

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

#' Decode a cell index to face, i, j, and resolution
#' 
#' Internal function to decode a cell index string into its constituent
#' components: face number, grid coordinates (i, j), and resolution level.
#' 
#' @param index Cell index string
#' @param aperture Grid aperture (3, 4, or 7)
#' @param index_type Index encoding type ("z3", "z7", or "zorder")
#' 
#' @return List with components:
#'   \item{face}{Face number (integer)}
#'   \item{i}{Grid coordinate i (integer)}
#'   \item{j}{Grid coordinate j (integer)}
#'   \item{resolution}{Resolution level (integer)}
#'   
#' @keywords internal
index_to_cell_internal <- function(index, aperture, index_type) {
  # Extract face from first 2 characters
  face <- as.integer(substr(index, 1, 2))

  # Extract the rest of the index
  index_body <- substr(index, 3, nchar(index))

  # Decode based on index type
  if (index_type == "z3") {
    # Z3 indexing for aperture 3
    result <- cpp_decode_z3(index_body, aperture)
  } else if (index_type == "z7") {
    # Z7 indexing for aperture 7 - expects the FULL index string
    result <- cpp_decode_z7(index, aperture)
    # cpp_decode_z7 returns quad, but we already extracted face
    return(list(
      face = result$quad,
      i = result$i,
      j = result$j,
      resolution = result$resolution
    ))
  } else {
    # Z-order (Morton) indexing for aperture 4
    result <- cpp_decode_zorder(index_body, aperture)
  }

  return(list(
    face = face,
    i = result$i,
    j = result$j,
    resolution = result$resolution
  ))
}

#' Round-trip accuracy test
#' 
#' Tests the accuracy of the coordinate conversion functions by converting
#' coordinates to cells and back, measuring the distance between original
#' and reconstructed coordinates.
#' 
#' @param grid Grid specification
#' @param lon Longitude to test
#' @param lat Latitude to test
#' @param units Distance units ("km" or "degrees")
#' 
#' @return List with:
#'   \item{original}{Original coordinates}
#'   \item{cell}{Cell index}
#'   \item{reconstructed}{Reconstructed coordinates}
#'   \item{error}{Distance between original and reconstructed}
#'
#' @family coordinate conversion
#' @export
hexify_roundtrip_test <- function(grid, lon, lat, units = "km") {
  # Convert to hierarchical index
  result <- hexify_lonlat_to_h_index(grid, lon, lat)

  # Convert back to coordinates
  coords <- hexify_h_index_to_lonlat(grid, result$h_index)

  # Calculate distance
  if (units == "km") {
    # Haversine distance
    R <- EARTH_RADIUS_KM
    dlat <- (coords$lat - lat) * pi / 180
    dlon <- (coords$lon - lon) * pi / 180
    lat_rad <- lat * pi / 180
    cen_lat_rad <- coords$lat * pi / 180
    a <- sin(dlat/2)^2 + cos(lat_rad) * cos(cen_lat_rad) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    distance <- R * c
  } else {
    # Euclidean distance in degrees
    distance <- sqrt((coords$lon - lon)^2 + (coords$lat - lat)^2)
  }

  return(list(
    original = c(lon = lon, lat = lat),
    h_index = result$h_index,
    reconstructed = c(lon = coords$lon, lat = coords$lat),
    error = distance,
    units = units
  ))
}

# =============================================================================
# QUAD COORDINATE ACCESS FUNCTIONS
# =============================================================================
#
# These functions expose the intermediate Quad IJ coordinate system used
# internally by DGGS. The coordinate pipeline is:
#
#   lon/lat → Icosa Triangle (face, tx, ty) → Quad IJ (quad, i, j) → Cell ID
#
# Most users should use the high-level hexify() function or cell ID-based
# functions. These Quad IJ functions are useful for:
# - Understanding the DGGS internals
# - Debugging coordinate transformations
# - Advanced applications requiring direct grid access
# =============================================================================

#' Convert longitude/latitude to Quad IJ coordinates
#'
#' Converts geographic coordinates to the intermediate Quad IJ representation
#' used internally by ISEA DGGS. Returns the quad number (0-11) and integer
#' cell indices (i, j) within that quad.
#'
#' @param lon Longitude in degrees (-180 to 180)
#' @param lat Latitude in degrees (-90 to 90)
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return List with components:
#'   \item{quad}{Quad number (0-11)}
#'   \item{i}{Integer cell index along first axis}
#'   \item{j}{Integer cell index along second axis}
#'   \item{icosa_triangle_face}{Source icosahedral face (0-19)}
#'   \item{icosa_triangle_x}{X coordinate on triangle face}
#'   \item{icosa_triangle_y}{Y coordinate on triangle face}
#'
#' @details
#' The 20 icosahedral triangle faces are grouped into 12 quads:
#' - Quad 0: North polar region
#' - Quads 1-5: Upper hemisphere rhombi
#' - Quads 6-10: Lower hemisphere rhombi
#' - Quad 11: South polar region
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # Get Quad IJ coordinates for Paris
#' result <- hexify_lonlat_to_quad_ij(lon = 2.35, lat = 48.86,
#'                                     resolution = 10, aperture = 3)
#' print(result)
hexify_lonlat_to_quad_ij <- function(lon, lat, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_lonlat_to_quad_ij(
    lon_deg = as.numeric(lon),
    lat_deg = as.numeric(lat),
    aperture = as.integer(aperture),
    resolution = as.integer(resolution)
  )
}

#' Convert Quad IJ coordinates to cell ID
#'
#' Converts Quad IJ coordinates to a global cell identifier.
#' This is the final step in the coordinate pipeline.
#'
#' @param quad Quad number (0-11), integer or vector
#' @param i Cell index along first axis, integer or vector
#' @param j Cell index along second axis, integer or vector
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Numeric vector of cell IDs
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # Convert Quad IJ to cell ID
#' cell_id <- hexify_quad_ij_to_cell(quad = 1, i = 100, j = 50,
#'                                    resolution = 10, aperture = 3)
#' print(cell_id)
hexify_quad_ij_to_cell <- function(quad, i, j, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_quad_ij_to_cell(
    quad = as.integer(quad),
    i = as.numeric(i),
    j = as.numeric(j),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

#' Convert Quad IJ to Quad XY (continuous coordinates)
#'
#' Converts discrete cell indices to continuous quad coordinates.
#' Useful for computing cell centers or understanding the cell geometry.
#'
#' @param quad Quad number (0-11)
#' @param i Cell index along first axis
#' @param j Cell index along second axis
#' @param resolution Grid resolution level
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return List with components:
#'   \item{quad_x}{Continuous X coordinate in quad space}
#'   \item{quad_y}{Continuous Y coordinate in quad space}
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # Get continuous quad coordinates for a cell
#' xy <- hexify_quad_ij_to_xy(quad = 1, i = 100, j = 50,
#'                            resolution = 10, aperture = 3)
#' print(xy)
hexify_quad_ij_to_xy <- function(quad, i, j, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }

  cpp_quad_ij_to_xy(
    quad = as.integer(quad),
    i = as.numeric(i),
    j = as.numeric(j),
    aperture = as.integer(aperture),
    resolution = as.integer(resolution)
  )
}

#' Convert Icosa Triangle to Quad XY coordinates
#'
#' Converts icosahedral triangle coordinates (from Snyder projection)
#' to quad XY coordinates. This is an intermediate step in the pipeline.
#'
#' @param icosa_triangle_face Triangle face number (0-19)
#' @param icosa_triangle_x X coordinate on triangle face
#' @param icosa_triangle_y Y coordinate on triangle face
#'
#' @return List with components:
#'   \item{quad}{Quad number (0-11)}
#'   \item{quad_x}{Continuous X coordinate in quad space}
#'   \item{quad_y}{Continuous Y coordinate in quad space}
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # First get triangle coordinates from lon/lat
#' fwd <- hexify_forward(lon = 2.35, lat = 48.86)
#'
#' # Then convert to quad XY
#' quad_xy <- hexify_icosa_tri_to_quad_xy(
#'   icosa_triangle_face = fwd["face"],
#'   icosa_triangle_x = fwd["icosa_triangle_x"],
#'   icosa_triangle_y = fwd["icosa_triangle_y"]
#' )
#' print(quad_xy)
hexify_icosa_tri_to_quad_xy <- function(icosa_triangle_face,
                                         icosa_triangle_x,
                                         icosa_triangle_y) {
  cpp_icosa_tri_to_quad_xy(
    icosa_triangle_face = as.integer(icosa_triangle_face),
    icosa_triangle_x = as.numeric(icosa_triangle_x),
    icosa_triangle_y = as.numeric(icosa_triangle_y)
  )
}

#' Convert Icosa Triangle to Quad IJ coordinates
#'
#' Converts icosahedral triangle coordinates directly to Quad IJ,
#' combining the transformation and quantization steps.
#'
#' @param icosa_triangle_face Triangle face number (0-19)
#' @param icosa_triangle_x X coordinate on triangle face
#' @param icosa_triangle_y Y coordinate on triangle face
#' @param resolution Grid resolution level
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return List with components:
#'   \item{quad}{Quad number (0-11)}
#'   \item{i}{Integer cell index along first axis}
#'   \item{j}{Integer cell index along second axis}
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # First get triangle coordinates from lon/lat
#' fwd <- hexify_forward(lon = 2.35, lat = 48.86)
#'
#' # Then convert to quad IJ
#' quad_ij <- hexify_icosa_tri_to_quad_ij(
#'   icosa_triangle_face = fwd["face"],
#'   icosa_triangle_x = fwd["icosa_triangle_x"],
#'   icosa_triangle_y = fwd["icosa_triangle_y"],
#'   resolution = 10,
#'   aperture = 3
#' )
#' print(quad_ij)
hexify_icosa_tri_to_quad_ij <- function(icosa_triangle_face,
                                         icosa_triangle_x,
                                         icosa_triangle_y,
                                         resolution,
                                         aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }

  cpp_icosa_tri_to_quad_ij(
    icosa_triangle_face = as.integer(icosa_triangle_face),
    icosa_triangle_x = as.numeric(icosa_triangle_x),
    icosa_triangle_y = as.numeric(icosa_triangle_y),
    aperture = as.integer(aperture),
    resolution = as.integer(resolution)
  )
}

#' Convert Quad XY to Icosa Triangle coordinates
#'
#' Inverse transformation from quad coordinates back to icosahedral
#' triangle coordinates. Useful for projecting cell centers back to lon/lat.
#'
#' @param quad Quad number (0-11)
#' @param quad_x Continuous X coordinate in quad space
#' @param quad_y Continuous Y coordinate in quad space
#'
#' @return List with components:
#'   \item{icosa_triangle_face}{Triangle face number (0-19)}
#'   \item{icosa_triangle_x}{X coordinate on triangle face}
#'   \item{icosa_triangle_y}{Y coordinate on triangle face}
#'
#' @family coordinate conversion
#' @keywords internal
#' @export
#' @examples
#' # Convert quad XY back to triangle coordinates
#' tri <- hexify_quad_xy_to_icosa_tri(quad = 1, quad_x = 0.5, quad_y = 0.3)
#' print(tri)
hexify_quad_xy_to_icosa_tri <- function(quad, quad_x, quad_y) {
  cpp_quad_xy_to_icosa_tri(
    quad = as.integer(quad),
    quad_x = as.numeric(quad_x),
    quad_y = as.numeric(quad_y)
  )
}

# =============================================================================
# CELL ID TO COORDINATE SYSTEM CONVERSIONS
# =============================================================================
#
# These functions convert from Cell IDs back to intermediate
# coordinate representations, compatible with 'dggridR' conversion functions.
# =============================================================================

#' Convert Cell ID to Quad IJ coordinates
#'
#' Converts DGGRID-compatible cell IDs to Quad IJ coordinates.
#' This is the inverse of hexify_quad_ij_to_cell().
#'
#' Compatible with 'dggridR' dgSEQNUM_to_Q2DI().
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Data frame with columns:
#'   \item{quad}{Quad number (0-11)}
#'   \item{i}{Integer cell index along first axis}
#'   \item{j}{Integer cell index along second axis}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_quad_ij_to_cell}} for the forward operation,
#'   \code{\link{hexify_cell_to_icosa_tri}} for conversion to triangle coords
#'
#' @keywords internal
#' @export
#' @examples
#' # Get Quad IJ coordinates for a cell
#' result <- hexify_cell_to_quad_ij(cell_id = 1000, resolution = 10, aperture = 3)
#' print(result)
#'
#' # Round-trip test
#' cell_id <- hexify_quad_ij_to_cell(result$quad, result$i, result$j,
#'                                    resolution = 10, aperture = 3)
#' # Should equal original cell_id
hexify_cell_to_quad_ij <- function(cell_id, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_cell_to_quad_ij(
    cell_id = as.numeric(cell_id),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

#' Convert Cell ID to Icosa Triangle coordinates
#'
#' Converts DGGRID-compatible cell IDs to icosahedral triangle
#' coordinates (face, x, y). These are the coordinates produced by the
#' 'Snyder' 'ISEA' forward projection.
#'
#' Compatible with 'dggridR' dgSEQNUM_to_PROJTRI().
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Data frame with columns:
#'   \item{icosa_triangle_face}{Triangle face number (0-19)}
#'   \item{icosa_triangle_x}{X coordinate on triangle face}
#'   \item{icosa_triangle_y}{Y coordinate on triangle face}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_cell_to_quad_ij}} for conversion to Quad IJ,
#'   \code{\link{hexify_cell_to_lonlat}} for conversion to lon/lat
#'
#' @keywords internal
#' @export
#' @examples
#' # Get triangle coordinates for a cell
#' result <- hexify_cell_to_icosa_tri(cell_id = 1000, resolution = 10, aperture = 3)
#' print(result)
#'
#' # Convert back to lon/lat using inverse projection
#' coords <- hexify_inverse(result$icosa_triangle_face,
#'                          result$icosa_triangle_x,
#'                          result$icosa_triangle_y)
hexify_cell_to_icosa_tri <- function(cell_id, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_cell_to_icosa_tri(
    cell_id = as.numeric(cell_id),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

#' Convert Quad IJ to Icosa Triangle coordinates
#'
#' Converts Quad IJ coordinates to icosahedral triangle coordinates.
#' This is useful for understanding where a cell is located on the
#' icosahedral projection.
#'
#' Equivalent to 'dggridR' dgQ2DI_to_PROJTRI().
#'
#' @param quad Quad number (0-11), integer or vector
#' @param i Cell index along first axis, integer or vector
#' @param j Cell index along second axis, integer or vector
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Data frame with columns:
#'   \item{icosa_triangle_face}{Triangle face number (0-19)}
#'   \item{icosa_triangle_x}{X coordinate on triangle face}
#'   \item{icosa_triangle_y}{Y coordinate on triangle face}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_icosa_tri_to_quad_ij}} for the inverse,
#'   \code{\link{hexify_cell_to_icosa_tri}} for conversion from cell ID
#'
#' @keywords internal
#' @export
#' @examples
#' # Get triangle coordinates for a Quad IJ position
#' result <- hexify_quad_ij_to_icosa_tri(quad = 1, i = 100, j = 50,
#'                                        resolution = 10, aperture = 3)
#' print(result)
hexify_quad_ij_to_icosa_tri <- function(quad, i, j, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_quad_ij_to_icosa_tri(
    quad = as.integer(quad),
    i = as.numeric(i),
    j = as.numeric(j),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

# =============================================================================
# QUAD XY (CONTINUOUS) COORDINATE CONVERSIONS
# =============================================================================
#
# These functions work with Quad XY coordinates - the continuous (non-quantized)
# representation in quad space. Equivalent to 'dggridR' Q2DD coordinate system.
# =============================================================================

#' Convert Cell ID to Quad XY coordinates
#'
#' Converts DGGRID-compatible cell IDs to Quad XY coordinates
#' (continuous quad space). This is the cell center in quad coordinates.
#'
#' Compatible with 'dggridR' dgSEQNUM_to_Q2DD().
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Data frame with columns:
#'   \item{quad}{Quad number (0-11)}
#'   \item{quad_x}{Continuous X coordinate in quad space}
#'   \item{quad_y}{Continuous Y coordinate in quad space}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_quad_xy_to_cell}} for the inverse operation,
#'   \code{\link{hexify_cell_to_quad_ij}} for integer grid coordinates
#'
#' @keywords internal
#' @export
#' @examples
#' # Get Quad XY coordinates for a cell
#' result <- hexify_cell_to_quad_xy(cell_id = 1000, resolution = 10, aperture = 3)
#' print(result)
#'
#' # Round-trip test
#' cell_id <- hexify_quad_xy_to_cell(result$quad, result$quad_x, result$quad_y,
#'                                    resolution = 10, aperture = 3)
#' # Should equal original cell_id
hexify_cell_to_quad_xy <- function(cell_id, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_cell_to_quad_xy(
    cell_id = as.numeric(cell_id),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

#' Convert Quad XY coordinates to Cell ID
#'
#' Converts Quad XY coordinates (continuous quad space) to DGGRID-compatible
#' cell IDs. The coordinates are quantized to the nearest cell.
#'
#' Compatible with 'dggridR' dgQ2DD_to_SEQNUM().
#'
#' @param quad Quad number (0-11), integer or vector
#' @param quad_x Continuous X coordinate in quad space
#' @param quad_y Continuous Y coordinate in quad space
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Numeric vector of cell IDs
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_cell_to_quad_xy}} for the inverse operation,
#'   \code{\link{hexify_quad_ij_to_cell}} for integer grid coordinates
#'
#' @keywords internal
#' @export
#' @examples
#' # Convert Quad XY to cell ID
#' cell_id <- hexify_quad_xy_to_cell(quad = 1, quad_x = 0.5, quad_y = 0.3,
#'                                    resolution = 10, aperture = 3)
#' print(cell_id)
hexify_quad_xy_to_cell <- function(quad, quad_x, quad_y, resolution,
                                   aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_quad_xy_to_cell(
    quad = as.integer(quad),
    quad_x = as.numeric(quad_x),
    quad_y = as.numeric(quad_y),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

# =============================================================================
# PLANE COORDINATE CONVERSIONS
# =============================================================================
#
# PLANE coordinates represent the unfolded icosahedron in 2D space.
# Each triangle face is rotated and translated to form a contiguous layout
# covering approximately 5.5 x 1.73 units.
#
# These functions are useful for:
# - Visualizing the entire DGGS on a flat surface
# - Understanding spatial relationships between triangles
# - Creating custom map projections
# =============================================================================

#' Convert Icosa Triangle coordinates to PLANE coordinates
#'
#' Transforms icosahedral triangle coordinates to the 2D PLANE representation
#' (unfolded icosahedron). Each triangle is rotated and translated to its
#' position in the unfolded layout.
#'
#' Equivalent to 'dggridR' dgPROJTRI_to_PLANE().
#'
#' @param icosa_triangle_face Triangle face number (0-19), integer or vector
#' @param icosa_triangle_x X coordinate on triangle face
#' @param icosa_triangle_y Y coordinate on triangle face
#'
#' @return Data frame with columns:
#'   \item{plane_x}{X coordinate in PLANE space (range ~0 to 5.5)}
#'   \item{plane_y}{Y coordinate in PLANE space (range ~0 to 1.73)}
#'
#' @details
#' The PLANE layout arranges all 20 icosahedral faces into a roughly
#' rectangular region. Faces 0-4 and 5-9 form the upper row, while
#' faces 10-14 and 15-19 form the lower row. Adjacent faces share
#' edges in this representation.
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_cell_to_plane}} for direct cell ID conversion,
#'   \code{\link{hexify_lonlat_to_plane}} for lon/lat to PLANE
#'
#' @keywords internal
#' @export
#' @examples
#' # Get PLANE coordinates from triangle coordinates
#' fwd <- hexify_forward(lon = 2.35, lat = 48.86)
#' plane <- hexify_icosa_tri_to_plane(
#'   icosa_triangle_face = fwd["face"],
#'   icosa_triangle_x = fwd["icosa_triangle_x"],
#'   icosa_triangle_y = fwd["icosa_triangle_y"]
#' )
#' print(plane)
hexify_icosa_tri_to_plane <- function(icosa_triangle_face,
                                       icosa_triangle_x,
                                       icosa_triangle_y) {
  cpp_icosa_tri_to_plane(
    icosa_triangle_face = as.integer(icosa_triangle_face),
    icosa_triangle_x = as.numeric(icosa_triangle_x),
    icosa_triangle_y = as.numeric(icosa_triangle_y)
  )
}

#' Convert Cell ID to PLANE coordinates
#'
#' Converts DGGRID-compatible cell IDs directly to PLANE coordinates.
#' Returns the cell center in the unfolded icosahedron layout.
#'
#' Compatible with 'dggridR' dgSEQNUM_to_PLANE().
#'
#' @param cell_id Numeric vector of cell IDs (1-based)
#' @param resolution Grid resolution level (0-30)
#' @param aperture Grid aperture: 3, 4, or 7
#'
#' @return Data frame with columns:
#'   \item{plane_x}{X coordinate in PLANE space (range ~0 to 5.5)}
#'   \item{plane_y}{Y coordinate in PLANE space (range ~0 to 1.73)}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_icosa_tri_to_plane}} for triangle conversion,
#'   \code{\link{hexify_lonlat_to_plane}} for lon/lat conversion
#'
#' @keywords internal
#' @export
#' @examples
#' # Get PLANE coordinates for cells
#' plane <- hexify_cell_to_plane(cell_id = c(100, 200, 300),
#'                                resolution = 5, aperture = 3)
#' plot(plane$plane_x, plane$plane_y)
hexify_cell_to_plane <- function(cell_id, resolution, aperture = 3L) {

  if (!aperture %in% c(3L, 4L, 7L)) {
    stop("aperture must be 3, 4, or 7")
  }
  if (resolution < 0 || resolution > 30) {
    stop("resolution must be between 0 and 30")
  }

  cpp_cell_to_plane(
    cell_id = as.numeric(cell_id),
    resolution = as.integer(resolution),
    aperture = as.integer(aperture)
  )
}

#' Convert longitude/latitude to PLANE coordinates
#'
#' Converts geographic coordinates directly to PLANE coordinates
#' (unfolded icosahedron). Combines forward 'Snyder' projection with
#' the PLANE transformation.
#'
#' Equivalent to 'dggridR' dgGEO_to_PLANE().
#'
#' @param lon Longitude in degrees (-180 to 180)
#' @param lat Latitude in degrees (-90 to 90)
#'
#' @return Data frame with columns:
#'   \item{plane_x}{X coordinate in PLANE space (range ~0 to 5.5)}
#'   \item{plane_y}{Y coordinate in PLANE space (range ~0 to 1.73)}
#'
#' @family coordinate conversion
#' @seealso \code{\link{hexify_cell_to_plane}} for cell ID conversion,
#'   \code{\link{hexify_icosa_tri_to_plane}} for triangle conversion
#'
#' @keywords internal
#' @export
#' @examples
#' # Plot world cities in PLANE coordinates
#' cities <- data.frame(
#'   lon = c(2.35, -74.00, 139.69, 151.21),
#'   lat = c(48.86, 40.71, 35.69, -33.87)
#' )
#' plane <- hexify_lonlat_to_plane(cities$lon, cities$lat)
#' plot(plane$plane_x, plane$plane_y)
hexify_lonlat_to_plane <- function(lon, lat) {
  cpp_lonlat_to_plane(
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )
}
