# hexify_conversions.R
# Coordinate system conversion functions
#
# This file contains functions for converting between geographic coordinates
# (longitude/latitude) and hexagonal cell indices.

#' @title Coordinate Conversions
#' @description Functions for converting between coordinate systems
#' @name hexify-conversions
NULL

#' Convert longitude/latitude to hexagonal cell indices
#' 
#' Converts geographic coordinates (longitude, latitude) to hexagonal cell
#' indices in the discrete global grid system. This is the primary function
#' for geocoding points to grid cells.
#' 
#' @param grid Grid specification from hexify_construct()
#' @param lon Longitude vector in degrees (numeric, -180 to 180)
#' @param lat Latitude vector in degrees (numeric, -90 to 90)
#' 
#' @return Data frame with columns:
#'   \item{seqnum}{Cell index (character string)}
#'   \item{face}{Icosahedron face number (integer, 0-19)}
#'   
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_construct(area = 1000, aperture = 3)
#' 
#' # Single point
#' cell <- hexify_geo_to_cell(grid, lon = 0, lat = 0)
#' print(cell$seqnum)
#' 
#' # Multiple points
#' lons <- c(0, 10, 20)
#' lats <- c(0, 45, -30)
#' cells <- hexify_geo_to_cell(grid, lons, lats)
#' print(cells)
#' }
hexify_geo_to_cell <- function(grid, lon, lat) {
  
  # Validate grid object
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_construct()")
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
    seqnum = cell_indices,
    face = faces,
    stringsAsFactors = FALSE
  ))
}

#' Convert hexagonal cell indices to longitude/latitude centers
#' 
#' Converts hexagonal cell indices back to geographic coordinates,
#' returning the center point of each cell. This is the inverse operation
#' of hexify_geo_to_cell().
#' 
#' @param grid Grid specification from hexify_construct()
#' @param seqnum Cell indices (character vector)
#' 
#' @return Data frame with columns:
#'   \item{lon_deg}{Longitude in degrees}
#'   \item{lat_deg}{Latitude in degrees}
#'   
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_construct(area = 1000, aperture = 3)
#' 
#' # Convert coordinates to cells
#' cells <- hexify_geo_to_cell(grid, lon = 5, lat = 45)
#' 
#' # Convert cells back to coordinates
#' coords <- hexify_cell_to_geo(grid, cells$seqnum)
#' print(coords)
#' 
#' # Should be close to original coordinates (5, 45)
#' # Exact match depends on cell size and resolution
#' }
hexify_cell_to_geo <- function(grid, seqnum) {
  
  # Validate grid object
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_construct()")
  }
  
  # Validate input
  if (!is.character(seqnum)) {
    stop("seqnum must be a character vector")
  }
  
  n <- length(seqnum)
  
  # Preallocate result vectors
  lon_deg <- numeric(n)
  lat_deg <- numeric(n)
  
  # Process each cell index
  for (i in seq_along(seqnum)) {
    # Handle NA values
    if (is.na(seqnum[i])) {
      lon_deg[i] <- NA_real_
      lat_deg[i] <- NA_real_
      next
    }
    
    # Call C++ conversion function
    # This function:
    # 1. Decodes the index string to face, i, j, resolution
    # 2. Calculates hex center in face plane coordinates
    # 3. Projects back to lon/lat using inverse ISEA projection
    result <- cpp_index_to_lonlat(
      seqnum[i],
      grid$aperture,
      grid$index_type
    )
    
    # Extract lon/lat from named vector
    lon_deg[i] <- result["lon"]
    lat_deg[i] <- result["lat"]
  }
  
  return(data.frame(
    lon_deg = lon_deg,
    lat_deg = lat_deg,
    stringsAsFactors = FALSE
  ))
}

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
#' @examples
#' \dontrun{
#' # This is an internal function, typically called by other functions
#' cell_info <- index_to_cell("0112345", aperture = 3, index_type = "z3")
#' print(cell_info)
#' }
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
#' @export
#' @examples
#' \dontrun{
#' grid <- hexify_construct(area = 1000, aperture = 3)
#' accuracy <- hexify_roundtrip_test(grid, lon = 0, lat = 45)
#' print(accuracy)
#' }
hexify_roundtrip_test <- function(grid, lon, lat, units = "km") {
  # Convert to cell
  cell <- hexify_geo_to_cell(grid, lon, lat)
  
  # Convert back to coordinates
  coords <- hexify_cell_to_geo(grid, cell$seqnum)
  
  # Calculate distance
  if (units == "km") {
    # Haversine distance
    R <- 6371  # Earth radius in km
    dlat <- (coords$lat_deg - lat) * pi / 180
    dlon <- (coords$lon_deg - lon) * pi / 180
    a <- sin(dlat/2)^2 + cos(lat * pi/180) * cos(coords$lat_deg * pi/180) * sin(dlon/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    distance <- R * c
  } else {
    # Euclidean distance in degrees
    distance <- sqrt((coords$lon_deg - lon)^2 + (coords$lat_deg - lat)^2)
  }
  
  return(list(
    original = c(lon = lon, lat = lat),
    cell = cell$seqnum,
    reconstructed = c(lon = coords$lon_deg, lat = coords$lat_deg),
    error = distance,
    units = units
  ))
}
