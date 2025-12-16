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
#' parameters needed for grid operations so they don't need to be repeated
#' in downstream function calls.
#'
#' @slot aperture Character. Grid aperture: "3", "4", "7", or "4/3" for mixed.
#' @slot resolution Integer. Grid resolution level (0-30).
#' @slot area_km2 Numeric. Target cell area in square kilometers.
#' @slot grid_system Character. Grid system identifier (default "ISEA").
#' @slot topology Character. Grid topology (default "H" for hexagon).
#' @slot index_type Character. Index encoding type ("integer" or "character").
#' @slot crs_input Integer. Input coordinate reference system (default 4326).
#' @slot crs_work Integer. Working CRS for internal operations.
#' @slot meta List. Additional metadata for future extensions.
#'
#' @param name Slot name for $ access
#' @param object A HexGrid object (for show method)
#' @param ... Additional arguments (ignored)
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
    grid_system = "character",
    topology = "character",
    index_type = "character",
    crs_input = "integer",
    crs_work = "integer",
    meta = "list"
  ),
  prototype = list(
    aperture = "3",
    resolution = 0L,
    area_km2 = NA_real_,
    grid_system = "ISEA",
    topology = "H",
    index_type = "integer",
    crs_input = 4326L,
    crs_work = 4326L,
    meta = list()
  )
)

# =============================================================================
# S4 CLASS: HexData
# =============================================================================
#
# HexData wraps user data with a reference to the HexGrid specification.
# The underlying data remains tabular and compatible with dplyr workflows.
# =============================================================================

#' HexData Class
#'
#' An S4 class representing hexified data. Wraps the user's data with a
#' reference to the grid specification used, enabling downstream operations
#' without repeated parameter specification.
#'
#' @slot data Data frame or sf object. The underlying data with cell assignments.
#' @slot grid HexGrid object. The grid specification used for hexification.
#' @slot mapping List. Column name mappings (lon, lat, geometry columns used).
#' @slot kind Character. Data type: "points", "cells", or "unknown".
#' @slot meta List. Additional metadata (e.g., cached polygons).
#'
#' @param x A HexData object
#' @param name Column name for $ access
#' @param i Row indices or logical vector for subsetting
#' @param j Column indices, names, or logical vector for subsetting
#' @param ... Additional arguments passed to underlying methods
#' @param drop Whether to drop dimensions when subsetting (default FALSE)
#' @param object A HexData object (for show/as.data.frame methods)
#' @param row.names Row names for as.data.frame conversion
#' @param optional Logical; if TRUE, row.names may be omitted
#' @param value Value to assign
#'
#' @details
#' HexData objects are created by \code{\link{hexify}} and can be used with
#' standard R functions. Use \code{as.data.frame()} to extract the underlying
#' data as a plain data frame.
#'
#' @section Compatibility:
#' HexData objects preserve the structure of the underlying data:
#' \itemize{
#'   \item If input was data.frame, output data slot is data.frame
#'   \item If input was sf, output data slot is sf (geometry preserved)
#'   \item Subsetting operations work transparently via `[` method
#' }
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
    mapping = "list",
    kind = "character",
    meta = "list"
  ),
  prototype = list(
    data = data.frame(),
    grid = new("HexGrid"),
    mapping = list(),
    kind = "unknown",
    meta = list()
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

  # Validate topology
  if (!object@topology %in% c("H", "HEXAGON")) {
    errors <- c(errors, "topology must be 'H' or 'HEXAGON'")
  }

  # Validate grid_system
  if (!object@grid_system %in% c("ISEA")) {
    errors <- c(errors, "grid_system must be 'ISEA'")
  }

  # Validate index_type
  if (!object@index_type %in% c("integer", "character")) {
    errors <- c(errors, "index_type must be 'integer' or 'character'")
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

  # Check kind is valid
  if (!object@kind %in% c("points", "cells", "unknown")) {
    errors <- c(errors, "kind must be 'points', 'cells', or 'unknown'")
  }

  # Check required columns exist for hexified data
  if (nrow(object@data) > 0 && !"cell_id" %in% names(object@data)) {
    errors <- c(errors, "data must contain 'cell_id' column")
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
  unique(x@data$cell_id)
})

#' @describeIn HexData-class Count unique cells
#' @export
setMethod("n_cells", "HexData", function(x) {
  length(unique(x@data$cell_id))
})

#' @describeIn HexData-class Get number of rows
#' @export
setMethod("nrow", "HexData", function(x) {
  nrow(x@data)
})

#' @describeIn HexData-class Get number of columns
#' @export
setMethod("ncol", "HexData", function(x) {
  ncol(x@data)
})

#' @describeIn HexData-class Get dimensions
#' @export
setMethod("dim", "HexData", function(x) {
  dim(x@data)
})

#' @describeIn HexData-class Get column names
#' @export
setMethod("names", "HexData", function(x) {
  names(x@data)
})

#' @describeIn HexData-class Access columns via $
#' @export
setMethod("$", "HexData", function(x, name) {
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
    new("HexData",
        data = new_data,
        grid = x@grid,
        mapping = x@mapping,
        kind = x@kind,
        meta = x@meta)
  } else {
    # If subset extracted a vector, return it directly
    new_data
  }
})

#' @describeIn HexData-class Subset single column
#' @export
setMethod("[[", c("HexData", "ANY"), function(x, i) {
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

  cat(sprintf("Grid System: %s\n", object@grid_system))
  cat(sprintf("Topology:    %s\n", object@topology))
  cat(sprintf("Index Type:  %s\n", object@index_type))

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
  cat(sprintf("Kind:    %s\n", object@kind))

  if (inherits(object@data, "sf")) {
    cat("Type:    sf (spatial features)\
")
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

  # Show first few rows
  if (nrow(object@data) > 0) {
    cat("\nData preview:\n")
    # Select key columns for preview
    preview_cols <- intersect(
      c("cell_id", "cell_cen_lon", "cell_cen_lat", names(object@data)[1:3]),
      names(object@data)
    )
    preview_cols <- unique(preview_cols)[1:min(5, length(unique(preview_cols)))]

    preview <- head(object@data[, preview_cols, drop = FALSE], 3)
    print(preview, row.names = FALSE)

    if (nrow(object@data) > 3) {
      cat(sprintf("... with %d more rows\n", nrow(object@data) - 3))
    }
  }

  invisible(object)
})

# =============================================================================
# COERCION METHODS
# =============================================================================

#' @describeIn HexData-class Convert to data.frame
#' @export
setMethod("as.data.frame", "HexData", function(x, row.names = NULL,
                                                optional = FALSE, ...) {
  df <- x@data
  if (inherits(df, "sf")) {
    df <- as.data.frame(sf::st_drop_geometry(df))
  }
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
    grid_system = x@grid_system,
    topology = x@topology,
    index_type = x@index_type,
    crs_input = x@crs_input,
    crs_work = x@crs_work,
    meta = x@meta
  )
})

#' @describeIn HexData-class Convert to list
#' @export
setMethod("as.list", "HexData", function(x, ...) {
  list(
    data = x@data,
    grid = as.list(x@grid),
    mapping = x@mapping,
    kind = x@kind,
    meta = x@meta
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
  # Determine index_type
  idx_type <- if (!is.null(x$index_type)) {
    if (x$index_type %in% c("z3", "z7", "zorder")) "integer" else x$index_type
  } else {
    "integer"
  }

  new("HexGrid",
      aperture = as.character(x$aperture),
      resolution = as.integer(x$resolution),
      area_km2 = if (!is.null(x$area)) as.numeric(x$area) else NA_real_,
      grid_system = if (!is.null(x$projection)) x$projection else "ISEA",
      topology = "H",
      index_type = idx_type,
      crs_input = 4326L,
      crs_work = 4326L,
      meta = list(legacy = TRUE))
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
    projection = x@grid_system,
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
