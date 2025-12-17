# classes.R
# S4 class definitions for HexGrid and HexData
#
# This file defines the core S4 classes that provide stateful return objects
# for the hexify package, enabling cleaner workflows without repeated parameters.

#' @import methods
#' @importFrom methods setClass setMethod setGeneric setValidity
#' @importFrom methods new slot slotNames validObject
NULL

# =============================================================================
# S4 CLASS: HexGrid
# =============================================================================
#
# HexGrid stores grid specification parameters (aperture, resolution, etc.)
# so downstream functions don't need them repeated.
# =============================================================================

#' HexGrid Class
#'
#' An S4 class representing a hexagonal grid specification. Stores all
#' parameters needed for grid operations.
#'
#' @slot aperture Character. Grid aperture: "3", "4", "7", or "4/3" for mixed.
#' @slot resolution Integer. Grid resolution level (0-30).
#' @slot area_km2 Numeric. Cell area in square kilometers.
#' @slot diagonal_km Numeric. Cell diagonal (long diagonal) in kilometers.
#' @slot crs Integer. Coordinate reference system (default 4326 = WGS84).
#'
#' @details
#' Create HexGrid objects using the \code{\link{hex_grid}} constructor function.
#' Do not use \code{new("HexGrid", ...)} directly.
#'
#' The aperture can be "3", "4", "7" for standard grids, or "4/3" for mixed
#' aperture grids that start with aperture 4 and switch to aperture 3.
#'
#' @seealso \code{\link{hex_grid}} for the constructor function,
#'   \code{\link{HexData-class}} for hexified data objects
#'
#' @exportClass HexGrid
setClass(
  "HexGrid",
  slots = c(
    aperture = "character",
    resolution = "integer",
    area_km2 = "numeric",
    diagonal_km = "numeric",
    crs = "integer"
  ),
  prototype = list(
    aperture = "3",
    resolution = 0L,
    area_km2 = NA_real_,
    diagonal_km = NA_real_,
    crs = 4326L
  )
)

# =============================================================================
# S4 CLASS: HexData
# =============================================================================
#
# HexData wraps user data with cell assignments from hexification.
# Original data is preserved; cell info stored separately.
# =============================================================================

#' HexData Class
#'
#' An S4 class representing hexified data. Contains the original user data
#' plus cell assignments from the hexification process.
#'
#' @slot data Data frame or sf object. The original user data (untouched).
#' @slot grid HexGrid object. The grid specification used.
#' @slot cell_id Numeric vector. Cell IDs for each row of data.
#' @slot cell_center Matrix. Two-column matrix (lon, lat) of cell centers.
#'
#' @details
#' HexData objects are created by \code{\link{hexify}}. The original data
#' is preserved in the \code{data} slot, while cell assignments are stored
#' separately in \code{cell_id} and \code{cell_center}.
#'
#' Use \code{as.data.frame()} to get a combined data frame with cell columns.
#'
#' @seealso \code{\link{hexify}} for creating HexData objects,
#'   \code{\link{HexGrid-class}} for grid specifications
#'
#' @exportClass HexData
setClass(
  "HexData",
  slots = c(
    data = "ANY",  # data.frame or sf
    grid = "HexGrid",
    cell_id = "numeric",
    cell_center = "matrix"
  ),
  prototype = list(
    data = data.frame(),
    grid = new("HexGrid"),
    cell_id = numeric(0),
    cell_center = matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("lon", "lat")))
  )
)

# =============================================================================
# VALIDITY METHODS
# =============================================================================

#' @noRd
setValidity("HexGrid", function(object) {

  errors <- character()

  # Validate aperture
  if (!object@aperture %in% c("3", "4", "7", "4/3")) {
    errors <- c(errors, "aperture must be '3', '4', '7', or '4/3'")
  }

  # Validate resolution
  if (object@resolution < 0L || object@resolution > 30L) {
    errors <- c(errors, "resolution must be between 0 and 30")
  }

  # Validate area_km2 (must be positive if provided)
  if (!is.na(object@area_km2) && object@area_km2 <= 0) {
    errors <- c(errors, "area_km2 must be positive")
  }

  # Validate diagonal_km (must be positive if provided)
  if (!is.na(object@diagonal_km) && object@diagonal_km <= 0) {
    errors <- c(errors, "diagonal_km must be positive")
  }

  # Validate crs (must be positive integer)
  if (object@crs <= 0L) {
    errors <- c(errors, "crs must be a positive integer EPSG code")
  }

  if (length(errors) == 0) TRUE else errors
})

#' @noRd
setValidity("HexData", function(object) {
  errors <- character()

  # Check data is valid type
  if (!inherits(object@data, "data.frame") && !inherits(object@data, "sf")) {
    errors <- c(errors, "data must be a data.frame or sf object")
  }

  # Check cell_id length matches data rows
  n_rows <- nrow(object@data)
  if (length(object@cell_id) != n_rows && length(object@cell_id) > 0) {
    errors <- c(errors, "cell_id length must match number of data rows")
  }

  # Check cell_center dimensions
  if (nrow(object@cell_center) != n_rows && nrow(object@cell_center) > 0) {
    errors <- c(errors, "cell_center rows must match number of data rows")
  }
  if (ncol(object@cell_center) != 2 && nrow(object@cell_center) > 0) {
    errors <- c(errors, "cell_center must have exactly 2 columns (lon, lat)")
  }

  if (length(errors) == 0) TRUE else errors
})

# =============================================================================
# GENERICS
# =============================================================================

#' Get Grid Specification
#'
#' Extract the grid specification from a HexData object.
#'
#' @param x A HexData object
#' @return A HexGrid object
#'
#' @export
#' @examples
#' \dontrun{
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#' grid_spec <- grid(result)
#' }
setGeneric("grid", function(x) standardGeneric("grid"))

#' Get Cell IDs
#'
#' Extract the unique cell IDs present in a HexData object.
#'
#' @param x A HexData object
#' @return A vector of cell IDs
#'
#' @export
setGeneric("cells", function(x) standardGeneric("cells"))

#' Get Number of Cells
#'
#' Get the number of unique cells in a HexData object.
#'
#' @param x A HexData object
#' @return Integer count of unique cells
#'
#' @export
setGeneric("n_cells", function(x) standardGeneric("n_cells"))

# =============================================================================
# ACCESSORS FOR HexGrid
# =============================================================================

#' @describeIn HexGrid-class Get aperture value
#' @param x HexGrid object
#' @export
setMethod("$", "HexGrid", function(x, name) {
  slot(x, name)
})

#' @describeIn HexGrid-class Get slot names
#' @export
setMethod("names", "HexGrid", function(x) {
  slotNames(x)
})

# =============================================================================
# ACCESSORS FOR HexData
# =============================================================================

#' @describeIn HexData-class Extract grid specification
#' @export
setMethod("grid", "HexData", function(x) {
  x@grid
})

#' @describeIn HexData-class Extract unique cell IDs
#' @export
setMethod("cells", "HexData", function(x) {
  unique(x@cell_id)
})

#' @describeIn HexData-class Count unique cells
#' @export
setMethod("n_cells", "HexData", function(x) {
  length(unique(x@cell_id))
})

#' @describeIn HexData-class Get number of rows
#' @export
setMethod("nrow", "HexData", function(x) {
  nrow(x@data)
})

#' @describeIn HexData-class Get number of columns (includes virtual cell columns)
#' @export
setMethod("ncol", "HexData", function(x) {
  ncol(x@data) + 5L  # +5 for cell_id, cell_cen_lon, cell_cen_lat, cell_area_km2, cell_diag_km
})

#' @describeIn HexData-class Get dimensions
#' @export
setMethod("dim", "HexData", function(x) {
  dim(x@data)
})

#' @describeIn HexData-class Get column names (includes virtual cell columns)
#' @export
setMethod("names", "HexData", function(x) {
  c(names(x@data), "cell_id", "cell_cen_lon", "cell_cen_lat", "cell_area_km2", "cell_diag_km")
})

#' @describeIn HexData-class Access columns via $ (includes virtual cell columns)
#' @export
setMethod("$", "HexData", function(x, name) {
  # Virtual cell columns
  if (name == "cell_id") {
    return(x@cell_id)
  }
  if (name == "cell_cen_lon") {
    return(x@cell_center[, "lon"])
  }
  if (name == "cell_cen_lat") {
    return(x@cell_center[, "lat"])
  }
  if (name == "cell_area_km2") {
    return(rep(x@grid@area_km2, nrow(x@data)))
  }
  if (name == "cell_diag_km") {
    return(rep(x@grid@diagonal_km, nrow(x@data)))
  }
  # Regular data columns
  x@data[[name]]
})

#' @describeIn HexData-class Set columns via $<-
#' @export
setMethod("$<-", "HexData", function(x, name, value) {
  x@data[[name]] <- value
  x
})

#' @describeIn HexData-class Subset rows/columns
#' @export
setMethod("[", c("HexData", "ANY", "ANY"), function(x, i, j, ..., drop = FALSE) {
  # Create new HexData with subsetted data
  new_data <- x@data[i, j, ..., drop = drop]

  # If result is still a data.frame/sf, return HexData
  if (inherits(new_data, "data.frame") || inherits(new_data, "sf")) {
    # Subset cell_id and cell_center if row indices provided
    if (!missing(i)) {
      new_cell_id <- x@cell_id[i]
      new_cell_center <- x@cell_center[i, , drop = FALSE]
    } else {
      new_cell_id <- x@cell_id
      new_cell_center <- x@cell_center
    }

    new("HexData",
        data = new_data,
        grid = x@grid,
        cell_id = new_cell_id,
        cell_center = new_cell_center)
  } else {
    # If subset extracted a vector, return it directly
    new_data
  }
})

#' @describeIn HexData-class Subset single column (includes virtual cell columns)
#' @export
setMethod("[[", c("HexData", "ANY"), function(x, i) {
  # Virtual cell columns by name
  if (is.character(i)) {
    if (i == "cell_id") return(x@cell_id)
    if (i == "cell_cen_lon") return(x@cell_center[, "lon"])
    if (i == "cell_cen_lat") return(x@cell_center[, "lat"])
    if (i == "cell_area_km2") return(rep(x@grid@area_km2, nrow(x@data)))
    if (i == "cell_diag_km") return(rep(x@grid@diagonal_km, nrow(x@data)))
  }
  x@data[[i]]
})

#' @describeIn HexData-class Set single column
#' @export
setMethod("[[<-", c("HexData", "ANY", "missing", "ANY"), function(x, i, j, value) {
  x@data[[i]] <- value
  x
})

# =============================================================================
# SHOW / PRINT METHODS
# =============================================================================

#' @describeIn HexGrid-class Print summary
#' @export
setMethod("show", "HexGrid", function(object) {
  cat("HexGrid Specification\n")
  cat("---------------------\n")
  cat(sprintf("Aperture:    %s\n", object@aperture))
  cat(sprintf("Resolution:  %d\n", object@resolution))

  if (!is.na(object@area_km2)) {
    cat(sprintf("Area:        %.2f km^2\n", object@area_km2))
  }
  if (!is.na(object@diagonal_km)) {
    cat(sprintf("Diagonal:    %.2f km\n", object@diagonal_km))
  }

  cat(sprintf("CRS:         EPSG:%d\n", object@crs))

  # Calculate total cells based on aperture
  if (object@aperture == "4/3") {
    # Mixed aperture: default to res/2 for level calculation
    level <- as.integer(object@resolution / 2)
    n_cells <- 10 * (4^level) * (3^(object@resolution - level)) + 2
  } else {
    ap <- as.integer(object@aperture)
    n_cells <- 10 * (ap^object@resolution) + 2
  }
  cat(sprintf("Total Cells: %.0f\n", n_cells))

  invisible(object)
})

#' @describeIn HexData-class Print summary
#' @export
setMethod("show", "HexData", function(object) {
  cat("HexData Object\n")
  cat("--------------\n")
  cat(sprintf("Rows:    %d\n", nrow(object@data)))
  cat(sprintf("Columns: %d\n", ncol(object@data)))
  cat(sprintf("Cells:   %d unique\n", n_cells(object)))

  if (inherits(object@data, "sf")) {
    cat("Type:    sf (spatial features)\n")
  } else {
    cat("Type:    data.frame\n")
  }

  cat("\nGrid:\n")
  cat(sprintf("  Aperture %s, Resolution %d",
              object@grid@aperture, object@grid@resolution))
  if (!is.na(object@grid@area_km2)) {
    cat(sprintf(" (~%.1f km^2)", object@grid@area_km2))
  }
  cat("\n")

  # Show column preview
  cat("\nColumns: ")
  col_names <- names(object@data)
  if (length(col_names) > 8) {
    cat(paste(col_names[1:8], collapse = ", "), ", ...\n")
  } else {
    cat(paste(col_names, collapse = ", "), "\n")
  }

  # Show first few rows with cell info

  if (nrow(object@data) > 0) {
    cat("\nData preview (with cell assignments):\n")
    # Combine data with cell info for preview
    preview_df <- data.frame(
      object@data[1:min(3, nrow(object@data)), 1:min(3, ncol(object@data)), drop = FALSE],
      cell_id = object@cell_id[1:min(3, length(object@cell_id))],
      check.names = FALSE
    )
    print(preview_df, row.names = FALSE)

    if (nrow(object@data) > 3) {
      cat(sprintf("... with %d more rows\n", nrow(object@data) - 3))
    }
  }

  invisible(object)
})

# =============================================================================
# COERCION METHODS
# =============================================================================

#' @describeIn HexData-class Convert to data.frame (includes cell columns)
#' @export
setMethod("as.data.frame", "HexData", function(x, row.names = NULL,
                                                optional = FALSE, ...) {
  df <- x@data
  if (inherits(df, "sf")) {
    df <- as.data.frame(sf::st_drop_geometry(df))
  }

  # Add cell columns
  df$cell_id <- x@cell_id
  df$cell_cen_lon <- x@cell_center[, "lon"]
  df$cell_cen_lat <- x@cell_center[, "lat"]
  df$cell_area_km2 <- x@grid@area_km2
  df$cell_diag_km <- x@grid@diagonal_km

  if (!is.null(row.names)) {
    rownames(df) <- row.names
  }
  df
})

#' @describeIn HexGrid-class Convert to list
#' @export
setMethod("as.list", "HexGrid", function(x, ...) {
  list(
    aperture = x@aperture,
    resolution = x@resolution,
    area_km2 = x@area_km2,
    diagonal_km = x@diagonal_km,
    crs = x@crs
  )
})

#' @describeIn HexData-class Convert to list
#' @export
setMethod("as.list", "HexData", function(x, ...) {
  list(
    data = x@data,
    grid = as.list(x@grid),
    cell_id = x@cell_id,
    cell_center = x@cell_center
  )
})

# =============================================================================
# HELPER FUNCTIONS FOR CLASS CONSTRUCTION
# =============================================================================

#' Check if object is HexGrid
#'
#' @param x Object to check
#' @return Logical
#' @export
is_hex_grid <- function(x) {
  inherits(x, "HexGrid")
}

#' Check if object is HexData
#'
#' @param x Object to check
#' @return Logical
#' @export
is_hex_data <- function(x) {
  inherits(x, "HexData")
}

#' Extract grid from various objects
#'
#' Internal function to extract a HexGrid from different input types.
#' Accepts HexGrid, HexData, or legacy hexify_grid objects.
#'
#' @param x Object containing grid info
#' @param allow_null If TRUE, return NULL when x is NULL
#' @return HexGrid object
#' @keywords internal
extract_grid <- function(x, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop("grid specification required")
  }

  if (is_hex_grid(x)) {
    return(x)
  }

  if (is_hex_data(x)) {
    return(x@grid)
  }

  # Handle legacy hexify_grid objects (S3 class)
  if (inherits(x, "hexify_grid")) {
    return(hexify_grid_to_HexGrid(x))
  }

  stop("Cannot extract grid from object of class ", class(x)[1])
}

#' Convert legacy hexify_grid to HexGrid
#'
#' @param x A hexify_grid object (S3)
#' @return A HexGrid object (S4)
#' @keywords internal
hexify_grid_to_HexGrid <- function(x) {
  area <- if (!is.null(x$area)) as.numeric(x$area) else NA_real_
  diagonal <- if (!is.na(area)) sqrt(area * 2 / sqrt(3)) else NA_real_

  new("HexGrid",
      aperture = as.character(x$aperture),
      resolution = as.integer(x$resolution),
      area_km2 = area,
      diagonal_km = diagonal,
      crs = 4326L)
}

#' Convert HexGrid to legacy hexify_grid
#'
#' For backwards compatibility with existing functions.
#'
#' @param x A HexGrid object (S4)
#' @return A hexify_grid object (S3)
#' @keywords internal
HexGrid_to_hexify_grid <- function(x) {
  # Determine legacy index_type based on aperture
  ap <- x@aperture
  legacy_index <- if (ap == "3") {
    "z3"
  } else if (ap == "7") {
    "z7"
  } else {
    "zorder"
  }

  # Convert aperture to numeric for legacy
  aperture_num <- if (ap == "4/3") 3L else as.integer(ap)

  grid <- list(
    area = x@area_km2,
    resolution = x@resolution,
    aperture = aperture_num,
    topology = "HEXAGON",
    projection = "ISEA",
    metric = TRUE,
    index_type = legacy_index,
    res = x@resolution,
    topology_family = "HEXAGON",
    metric_radius = if (!is.na(x@area_km2)) sqrt(x@area_km2 / pi) else NULL,
    pole_lon_deg = ISEA_VERT0_LON_DEG,
    pole_lat_deg = ISEA_VERT0_LAT_DEG,
    azimuth_deg = ISEA_AZIMUTH_DEG,
    aperture_type = if (ap == "4/3") "MIXED43" else "SEQUENCE",
    res_spec = x@resolution,
    precision = 7
  )

  class(grid) <- c("hexify_grid", "dggs", "list")
  grid
}
