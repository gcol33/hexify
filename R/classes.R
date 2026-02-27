# classes.R
# S4 class definitions for HexGridInfo and HexData
#
# This file defines the core S4 classes that provide stateful return objects
# for the hexify package, enabling cleaner workflows without repeated parameters.

#' @import methods
#' @importFrom methods setClass setMethod setGeneric setValidity
#' @importFrom methods new slot slotNames validObject
NULL

# =============================================================================
# S4 CLASS: HexGridInfo
# =============================================================================
#
# HexGridInfo stores grid specification parameters (aperture, resolution, etc.)
# so downstream functions don't need them repeated.
# =============================================================================

#' HexGridInfo Class
#'
#' An S4 class representing a hexagonal grid specification. Stores all
#' parameters needed for grid operations.
#'
#' @slot aperture Character. Grid aperture: "3", "4", "7", or "4/3" for mixed.
#' @slot resolution Integer. Grid resolution level (0-30 for ISEA, 0-15 for H3).
#' @slot area_km2 Numeric. Cell area in square kilometers.
#' @slot diagonal_km Numeric. Cell diagonal (long diagonal) in kilometers.
#' @slot crs Integer. Coordinate reference system (default 4326 = 'WGS84').
#' @slot grid_type Character. Grid system: "isea" (default) or "h3".
#'
#' @details
#' Create HexGridInfo objects using the \code{\link{hex_grid}} constructor function.
#' Do not use \code{new("HexGridInfo", ...)} directly.
#'
#' The aperture can be "3", "4", "7" for standard grids, or "4/3" for mixed
#' aperture grids that start with aperture 4 and switch to aperture 3.
#'
#' For H3 grids, the aperture is fixed at "7" and resolution ranges from 0 to 15.
#'
#' @seealso \code{\link{hex_grid}} for the constructor function,
#'   \code{\link{HexData-class}} for hexified data objects
#'
#' @exportClass HexGridInfo
setClass(
  "HexGridInfo",
  slots = c(
    aperture = "character",
    resolution = "integer",
    area_km2 = "numeric",
    diagonal_km = "numeric",
    crs = "integer",
    grid_type = "character"
  ),
  prototype = list(
    aperture = "3",
    resolution = 0L,
    area_km2 = NA_real_,
    diagonal_km = NA_real_,
    crs = 4326L,
    grid_type = "isea"
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
#' @slot grid HexGridInfo object. The grid specification used.
#' @slot cell_id Cell IDs for each row of data. Numeric for ISEA grids,
#'   character for H3 grids.
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
#'   \code{\link{HexGridInfo-class}} for grid specifications
#'
#' @exportClass HexData
setClass(
  "HexData",
  slots = c(
    data = "ANY",  # data.frame or sf
    grid = "HexGridInfo",
    cell_id = "ANY",  # numeric for ISEA, character for H3
    cell_center = "matrix"
  ),
  prototype = list(
    data = data.frame(),
    grid = new("HexGridInfo"),
    cell_id = numeric(0),
    cell_center = matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("lon", "lat")))
  )
)

# =============================================================================
# VALIDITY METHODS
# =============================================================================

#' @noRd
setValidity("HexGridInfo", function(object) {

  errors <- character()

  # Validate grid_type
  gt <- object@grid_type
  if (!gt %in% c("isea", "h3")) {
    errors <- c(errors, "grid_type must be 'isea' or 'h3'")
  }

  if (gt == "h3") {
    # H3 validation: aperture fixed at "7", resolution 0-15
    if (object@aperture != "7") {
      errors <- c(errors, "H3 grids must have aperture '7'")
    }
    if (object@resolution < 0L || object@resolution > 15L) {
      errors <- c(errors, "H3 resolution must be between 0 and 15")
    }
  } else {
    # ISEA validation
    if (!object@aperture %in% c("3", "4", "7", "4/3")) {
      errors <- c(errors, "aperture must be '3', '4', '7', or '4/3'")
    }
    if (object@resolution < 0L || object@resolution > 30L) {
      errors <- c(errors, "resolution must be between 0 and 30")
    }
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

  # Check cell_id type matches grid_type
  gt <- tryCatch(object@grid@grid_type, error = function(e) "isea")
  if (gt == "h3") {
    if (!is.character(object@cell_id) && length(object@cell_id) > 0) {
      errors <- c(errors, "H3 cell_id must be character")
    }
  } else {
    if (!is.numeric(object@cell_id) && length(object@cell_id) > 0) {
      errors <- c(errors, "ISEA cell_id must be numeric")
    }
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
#' @return A HexGridInfo object
#'
#' @export
#' @examples
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#' grid_spec <- grid_info(result)
setGeneric("grid_info", function(x) standardGeneric("grid_info"))

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
# ACCESSORS FOR HexGridInfo
# =============================================================================

#' HexGridInfo S4 Methods
#'
#' S4 methods for HexGridInfo objects. These provide standard R operations
#' like `$`, `names()`, `show()`, and `as.list()`.
#'
#' @name HexGridInfo-methods
#' @param x HexGridInfo object
#' @param name Slot name
#' @param object HexGridInfo object (for show)
#' @param ... Additional arguments
#' @return
#' - `$`: The value of the requested slot
#' - `names`: Character vector of slot names
#' - `show`: The object, invisibly (called for side effect of printing)
#' - `as.list`: A named list of slot values
#' @keywords internal
NULL

#' @rdname HexGridInfo-methods
#' @export
setMethod("$", "HexGridInfo", function(x, name) {
  slot(x, name)
})

#' @rdname HexGridInfo-methods
#' @keywords internal
#' @export
setMethod("names", "HexGridInfo", function(x) {
  slotNames(x)
})

# =============================================================================
# ACCESSORS FOR HexData
# =============================================================================

#' HexData S4 Methods
#'
#' S4 methods for HexData objects. These provide standard R operations
#' for accessing data, subsetting, and conversion.
#'
#' @name HexData-methods
#' @param x HexData object
#' @param name Column name
#' @param object HexData object (for show)
#' @param i,j Row/column indices
#' @param value Replacement value
#' @param drop Logical, whether to drop dimensions
#' @param row.names Optional row names
#' @param optional Logical (ignored)
#' @param ... Additional arguments
#' @return
#' - `grid_info`: HexGridInfo object containing grid specification
#' - `cells`: Numeric vector of unique cell IDs
#' - `n_cells`: Integer count of unique cells
#' - `nrow`, `ncol`, `dim`: Integer dimensions
#' - `names`: Character vector of column names (including virtual cell columns)
#' - `$`, `[[`: The requested column or cell data as a vector
#' - `$<-`, `[[<-`: The modified HexData object
#' - `[`: Subsetted HexData object or extracted data
#' - `show`: The object, invisibly (called for side effect of printing)
#' - `as.data.frame`: Data frame with original data plus cell columns
#' - `as.list`: Named list containing data, grid, cell_id, and cell_center
#' @keywords internal
NULL

#' @rdname HexData-methods
#' @export
setMethod("grid_info", "HexData", function(x) {
  x@grid
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("cells", "HexData", function(x) {
  unique(x@cell_id)
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("n_cells", "HexData", function(x) {
  length(unique(x@cell_id))
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("nrow", "HexData", function(x) {
  nrow(x@data)
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("ncol", "HexData", function(x) {
  ncol(x@data) + 5L  # +5 for cell_id, cell_cen_lon, cell_cen_lat, cell_area_km2, cell_diag_km
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("dim", "HexData", function(x) {
  dim(x@data)
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("names", "HexData", function(x) {
  c(names(x@data), "cell_id", "cell_cen_lon", "cell_cen_lat", "cell_area_km2", "cell_diag_km")
})

#' @rdname HexData-methods
#' @keywords internal
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
    if (is_h3_grid(x@grid)) {
      return(as.numeric(cell_area(grid = x)))
    }
    return(rep(x@grid@area_km2, nrow(x@data)))
  }
  if (name == "cell_diag_km") {
    return(rep(x@grid@diagonal_km, nrow(x@data)))
  }
  # Regular data columns
  x@data[[name]]
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("$<-", "HexData", function(x, name, value) {
  x@data[[name]] <- value
  x
})

#' @rdname HexData-methods
#' @keywords internal
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

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("[[", c("HexData", "ANY"), function(x, i) {
  # Virtual cell columns by name
  if (is.character(i)) {
    if (i == "cell_id") return(x@cell_id)
    if (i == "cell_cen_lon") return(x@cell_center[, "lon"])
    if (i == "cell_cen_lat") return(x@cell_center[, "lat"])
    if (i == "cell_area_km2") {
      if (is_h3_grid(x@grid)) return(as.numeric(cell_area(grid = x)))
      return(rep(x@grid@area_km2, nrow(x@data)))
    }
    if (i == "cell_diag_km") return(rep(x@grid@diagonal_km, nrow(x@data)))
  }
  x@data[[i]]
})

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("[[<-", c("HexData", "ANY", "missing", "ANY"), function(x, i, j, value) {
  x@data[[i]] <- value
  x
})

# =============================================================================
# SHOW / PRINT METHODS
# =============================================================================

#' @rdname HexGridInfo-methods
#' @keywords internal
#' @export
setMethod("show", "HexGridInfo", function(object) {
  gt <- tryCatch(object@grid_type, error = function(e) "isea")

  if (gt == "h3") {
    cat("HexGridInfo Specification [H3]\n")
    cat("-------------------------------\n")
    cat(sprintf("Grid Type:   H3 (Uber)\n"))
    cat(sprintf("Resolution:  %d\n", object@resolution))

    if (!is.na(object@area_km2)) {
      cat(sprintf("Avg Area:    %.4f km^2 (varies by location)\n", object@area_km2))
    }
    if (!is.na(object@diagonal_km)) {
      cat(sprintf("Avg Diagonal:%.2f km\n", object@diagonal_km))
    }

    cat(sprintf("CRS:         EPSG:%d\n", object@crs))

    h3_n_cells <- 2 + 120 * 7^object@resolution
    cat(sprintf("Total Cells: %.0f\n", h3_n_cells))
    cat("Note: H3 cells are NOT exactly equal-area\n")
  } else {
    cat("HexGridInfo Specification\n")
    cat("-------------------------\n")
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
      level <- as.integer(object@resolution / 2)
      n_cells <- 10 * (4^level) * (3^(object@resolution - level)) + 2
    } else {
      ap <- as.integer(object@aperture)
      n_cells <- 10 * (ap^object@resolution) + 2
    }
    cat(sprintf("Total Cells: %.0f\n", n_cells))
  }

  invisible(object)
})

#' @rdname HexData-methods
#' @keywords internal
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
  gt <- tryCatch(object@grid@grid_type, error = function(e) "isea")
  if (gt == "h3") {
    cat(sprintf("  H3 Resolution %d", object@grid@resolution))
    if (!is.na(object@grid@area_km2)) {
      cat(sprintf(" (~%.4f km^2 avg)", object@grid@area_km2))
    }
  } else {
    cat(sprintf("  Aperture %s, Resolution %d",
                object@grid@aperture, object@grid@resolution))
    if (!is.na(object@grid@area_km2)) {
      cat(sprintf(" (~%.1f km^2)", object@grid@area_km2))
    }
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

#' @rdname HexData-methods
#' @keywords internal
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
  if (is_h3_grid(x@grid)) {
    df$cell_area_km2 <- as.numeric(cell_area(grid = x))
  } else {
    df$cell_area_km2 <- x@grid@area_km2
  }
  df$cell_diag_km <- x@grid@diagonal_km

  if (!is.null(row.names)) {
    rownames(df) <- row.names
  }
  df
})

#' @rdname HexGridInfo-methods
#' @keywords internal
#' @export
setMethod("as.list", "HexGridInfo", function(x, ...) {
  list(
    aperture = x@aperture,
    resolution = x@resolution,
    area_km2 = x@area_km2,
    diagonal_km = x@diagonal_km,
    crs = x@crs,
    grid_type = x@grid_type
  )
})

#' @rdname HexData-methods
#' @keywords internal
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

#' Check if object is HexGridInfo
#'
#' @param x Object to check
#' @return Logical
#' @export
is_hex_grid <- function(x) {
  inherits(x, "HexGridInfo")
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
#' Internal function to extract a HexGridInfo from different input types.
#' Accepts HexGridInfo, HexData, or legacy hexify_grid objects.
#'
#' @param x Object containing grid info
#' @param allow_null If TRUE, return NULL when x is NULL
#' @return HexGridInfo object
#' @keywords internal
extract_grid <- function(x, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop("grid specification required")
  }

  if (is_hex_grid(x)) {
    # Handle deserialized old objects without grid_type slot
    if (!.hasSlot(x, "grid_type")) {
      x@grid_type <- "isea"
    }
    return(x)
  }

  if (is_hex_data(x)) {
    g <- x@grid
    if (!.hasSlot(g, "grid_type")) {
      g@grid_type <- "isea"
    }
    return(g)
  }

  # Handle legacy hexify_grid objects (S3 class)
  if (inherits(x, "hexify_grid")) {
    return(hexify_grid_to_HexGridInfo(x))
  }

  stop("Cannot extract grid from object of class ", class(x)[1])
}

#' Convert legacy hexify_grid to HexGridInfo
#'
#' @param x A hexify_grid object (S3)
#' @return A HexGridInfo object (S4)
#' @keywords internal
hexify_grid_to_HexGridInfo <- function(x) {
  area <- if (!is.null(x$area)) as.numeric(x$area) else NA_real_
  diagonal <- if (!is.na(area)) sqrt(area * 2 / sqrt(3)) else NA_real_

  new("HexGridInfo",
      aperture = as.character(x$aperture),
      resolution = as.integer(x$resolution),
      area_km2 = area,
      diagonal_km = diagonal,
      crs = 4326L)
}

#' Convert HexGridInfo to legacy hexify_grid
#'
#' For backwards compatibility with existing functions.
#'
#' @param x A HexGridInfo object (S4)
#' @return A hexify_grid object (S3)
#' @keywords internal
HexGridInfo_to_hexify_grid <- function(x) {
  # H3 grids cannot be converted to legacy format
 gt <- tryCatch(x@grid_type, error = function(e) "isea")
  if (gt == "h3") {
    stop("H3 grids cannot be converted to legacy hexify_grid format")
  }

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
