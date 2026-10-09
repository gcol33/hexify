# classes.R
# S4 class definitions for HexGridInfo and HexData
#
# This file defines the core S4 classes that provide stateful return objects
# for the hexify package, enabling cleaner workflows without repeated parameters.

#' @import methods
#' @importFrom methods setClass setMethod setGeneric setValidity
#' @importFrom methods new slot slotNames validObject
#' @include constants.R
NULL

# =============================================================================
# S4 CLASS: HexGridInfo
# =============================================================================
#
# HexGridInfo stores grid specification parameters (aperture, resolution, etc.)
# so downstream functions don't need them repeated.
# =============================================================================

# A grid CRS: an 'EPSG' code, or a 'PROJ' or 'WKT' string
setClassUnion("HexCRS", c("integer", "character"))

#' HexGridInfo Class
#'
#' An S4 class representing a hexagonal grid specification. Stores all
#' parameters needed for grid operations.
#'
#' @slot aperture Character. Grid aperture: "3", "4", "7", a mixed family such
#'   as "4/3" or "4/7", or one aperture per resolution level ("4,4,7,3").
#' @slot resolution Integer. Grid resolution level: for ISEA 0-30 for aperture 3, 0-29 for aperture 4, 0-21 for aperture 7
#'   on the icosahedron (cell IDs fit in 64 bits); 0-15 for H3.
#' @slot area_km2 Numeric. Cell area in square kilometers.
#' @slot diagonal_km Numeric. Centre spacing in kilometers: the short
#'   (flat-to-flat) diagonal of a regular hexagon of area \code{area_km2},
#'   \eqn{\sqrt{2A/\sqrt{3}}}.
#' @slot crs Integer or character. Coordinate reference system: an EPSG code,
#'   or a 'PROJ' or 'WKT' string. Defaults to 'WGS84' on Earth, and to a longlat
#'   CRS on the sphere of \code{radius_km} on any other body.
#' @slot grid_type Character. Grid system: "isea" (default) or "h3".
#' @slot radius_km Numeric. Radius of the body the grid covers, in kilometers.
#'   \code{NA} reads as Earth's mean radius.
#' @slot orientation Numeric. Where an ISEA grid's solid sits:
#'   \code{c(vert0_lon, vert0_lat, azimuth)} in degrees, vertex 0 and the
#'   azimuth of vertex 1 seen from it. Empty for H3 grids, whose orientation
#'   H3 fixes.
#' @slot projection Character. How an ISEA-family grid projects each face of
#'   its solid onto its plane triangle: "isea" (Snyder's equal-area
#'   projection), "ivea" (van Leeuwen and Strebe's vertex-oriented equal-area
#'   projection) or "fuller" (Fuller's projection). Empty for H3 grids.
#' @slot polyhedron Character. The solid an ISEA-family grid is built on:
#'   "icosahedron" or "octahedron". Empty for H3 grids.
#' @slot ellipsoid Numeric. The ellipsoid of revolution an ISEA-family grid
#'   reads geodetic latitude on, \code{c(a_km, f)}: semi-major axis in km and
#'   flattening. Empty for a grid that reads latitude on the sphere, and for
#'   H3 grids.
#'
#' @details
#' Create HexGridInfo objects using the \code{\link{hex_grid}} constructor function.
#' Do not use \code{new("HexGridInfo", ...)} directly.
#'
#' The aperture can be "3", "4", "7" for grids that refine by one aperture at
#' every level; a family name such as "4/3" or "4/7", which refines by the first
#' aperture for the first floor(resolution / 2) levels and by the second for the
#' rest; or one aperture per level, "4,4,7,3". The separator tells the two
#' apart: "4,7" and the family "4/7" are the same grid at resolution 2, but
#' their grids one level coarser are "4" and "7".
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
    crs = "HexCRS",
    grid_type = "character",
    radius_km = "numeric",
    orientation = "numeric",
    projection = "character",
    polyhedron = "character",
    ellipsoid = "numeric"
  ),
  prototype = list(
    aperture = "3",
    resolution = 0L,
    area_km2 = NA_real_,
    diagonal_km = NA_real_,
    crs = 4326L,
    grid_type = "isea",
    radius_km = NA_real_,
    orientation = ISEA_ORIENTATION,
    projection = "isea",
    polyhedron = "icosahedron",
    ellipsoid = numeric(0)
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
#' @slot cell_id Cell IDs for each row of data: \code{bit64::integer64} for
#'   ISEA grids, character for H3 grids.
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
    cell_id = "ANY",  # integer64 for ISEA, character for H3
    cell_center = "matrix"
  ),
  prototype = list(
    data = data.frame(),
    grid = new("HexGridInfo"),
    cell_id = bit64::integer64(0),
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
    if (length(object@orientation) != 0L) {
      errors <- c(errors, "H3 grids carry no orientation; H3 fixes its own")
    }
    if (length(object@projection) != 0L) {
      errors <- c(errors, "H3 grids carry no face projection; H3 fixes its own")
    }
    if (length(object@polyhedron) != 0L) {
      errors <- c(errors, "H3 grids carry no polyhedron; H3 is built on the icosahedron")
    }
    if (length(grid_ellipsoid(object)) != 0L) {
      errors <- c(errors, "H3 grids carry no ellipsoid; H3 reads latitude on the sphere")
    }
  } else {
    o <- object@orientation
    if (length(o) != 3L || !all(is.finite(o)) || o[2] < -90 || o[2] > 90) {
      errors <- c(errors, paste0("orientation must be c(vert0_lon, vert0_lat, ",
                                 "azimuth) in degrees, vert0_lat in [-90, 90]"))
    }
    if (length(object@projection) != 1L ||
        !object@projection %in% names(FACE_PROJECTIONS)) {
      errors <- c(errors, "projection must be \"isea\", \"fuller\" or \"ivea\"")
    }
    poly <- grid_polyhedron(object)
    if (length(poly) != 1L || !poly %in% GRID_POLYHEDRA) {
      errors <- c(errors, "polyhedron must be \"icosahedron\" or \"octahedron\"")
    } else if (length(object@projection) == 1L &&
               object@projection %in% ICOSAHEDRON_PROJECTIONS &&
               poly != "icosahedron") {
      errors <- c(errors, "Fuller's projection is defined on the icosahedron only")
    }
    e <- if (.hasSlot(object, "ellipsoid")) object@ellipsoid else numeric(0)
    if (length(e) != 0L &&
        (length(e) != 2L || !all(is.finite(e)) || e[1] <= 0 || e[2] < 0 || e[2] >= 1)) {
      errors <- c(errors, "ellipsoid must be empty or c(a_km, f) with a_km > 0 and 0 <= f < 1")
    }
    # ISEA validation
    ap_ok <- tryCatch({
      parse_aperture_seq(object@aperture, object@resolution)
      TRUE
    }, error = function(e) FALSE)
    if (!ap_ok) {
      errors <- c(errors, paste0("aperture must be 3, 4, 7, a family such as \"4/3\", ",
                                 "or one aperture per resolution level"))
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

  # Validate crs (an EPSG code, or a CRS string sf can read)
  if (is.character(object@crs)) {
    if (length(object@crs) != 1L || is.na(object@crs) || !nzchar(object@crs)) {
      errors <- c(errors, "crs must be a single non-empty CRS string")
    } else if (is.na(parse_crs(object@crs))) {
      errors <- c(errors, sprintf(
        "crs \"%s\" is not a coordinate reference system sf can read",
        object@crs
      ))
    }
  } else if (length(object@crs) != 1L || is.na(object@crs) || object@crs <= 0L) {
    errors <- c(errors, "crs must be a positive integer EPSG code")
  }

  # Validate radius_km (must be a positive finite scalar if provided)
  if (length(object@radius_km) != 1L) {
    errors <- c(errors, "radius_km must be a single number")
  } else if (!is.na(object@radius_km) &&
             (!is.finite(object@radius_km) || object@radius_km <= 0)) {
    errors <- c(errors, "radius_km must be positive")
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
    if (!bit64::is.integer64(object@cell_id)) {
      errors <- c(errors, "ISEA cell_id must be integer64")
    }
  }

  # Check cell_id length matches data rows. `length(cell_id) != n_rows` alone
  # already correctly allows the valid empty-prototype case (0 != 0 is
  # FALSE), so no extra "> 0" guard is needed -- and such a guard would hide
  # exactly the corrupt case this check exists to catch (data has rows but
  # cell_id is empty).
  n_rows <- nrow(object@data)
  if (length(object@cell_id) != n_rows) {
    errors <- c(errors, "cell_id length must match number of data rows")
  }

  # Check cell_center dimensions (same reasoning as cell_id above)
  if (nrow(object@cell_center) != n_rows) {
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
#' Counts cells: those a dataset occupies, or those a grid contains.
#'
#' @param x A HexData or HexGridInfo object
#' @return For a HexData, the number of distinct cells its rows fall in. For a
#'   HexGridInfo, the number of cells the grid divides the body into, which runs
#'   past integer range at fine resolutions and so comes back as a double.
#'
#' @export
#' @examples
#' grid <- hex_grid(area_km2 = 100000)
#' n_cells(grid)
#'
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' n_cells(hexify(df, lon = "lon", lat = "lat", grid = grid))
setGeneric("n_cells", function(x) standardGeneric("n_cells"))

#' Number of cells a grid contains
#'
#' Read by both `n_cells()` and the grid's own `summary()`, so the number a grid
#' reports and the number it prints are the same one.
#'
#' @param grid A HexGridInfo object
#' @return The cell count, as a double
#' @noRd
grid_n_cells <- function(grid) {
  gt <- tryCatch(grid@grid_type, error = function(e) "isea")

  if (gt == "h3") {
    h3_n_cells(grid@resolution)
  } else {
    aperture_n_cells(grid@aperture, grid@resolution, grid_polyhedron(grid))
  }
}

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
#' - `n_cells`: The number of cells the grid contains
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
setMethod("n_cells", "HexGridInfo", function(x) {
  grid_n_cells(x)
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
    return(unname(cell_area(grid = x)))
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
#' @details
#' Unlike \code{[.data.frame}, \code{drop} defaults to \code{FALSE}: selecting
#' a single column returns a \code{HexData} object (preserving \code{grid}/
#' \code{cell_id}/\code{cell_center}) rather than dropping to a bare vector.
#' Pass \code{drop = TRUE} explicitly to get data.frame-style dropping.
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
    if (i == "cell_area_km2") return(unname(cell_area(grid = x)))
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
  print(summary(object))
  invisible(object)
})

#' Summary of a grid or of gridded data
#'
#' Reports what printing the object reports, as a list, so that a caller can
#' read the cell count, area, diagonal or column names without reaching into
#' slots. Printing an object prints its summary, so the two agree.
#'
#' @param object A HexGridInfo or HexData object
#' @param x A summary, as \code{summary()} returns one
#' @param ... Ignored
#'
#' @return For a HexGridInfo, a list of class \code{hexify_grid_summary}
#'   carrying \code{grid_type}, \code{aperture}, \code{resolution},
#'   \code{area_km2}, \code{diagonal_km}, \code{crs}, \code{radius_km},
#'   \code{earth}, \code{orientation} (\code{c(vert0_lon, vert0_lat, azimuth)},
#'   empty for H3), \code{projection} (\code{"isea"}, \code{"ivea"} or \code{"fuller"},
#'   \code{NA} for H3), \code{polyhedron} (\code{"icosahedron"} or
#'   \code{"octahedron"}, \code{NA} for H3), \code{ellipsoid}
#'   (\code{c(a_km, f)}, empty for a grid on the sphere) and \code{n_cells}.
#'   For a HexData, a list of class
#'   \code{hexify_data_summary} carrying \code{rows}, \code{columns},
#'   \code{column_names}, \code{n_cells}, \code{type}, the \code{grid} summary
#'   and a \code{preview} of the first rows. The print methods return their
#'   input invisibly.
#'
#' @name hexify-summary
#' @examples
#' grid <- hex_grid(area_km2 = 100000)
#' summary(grid)$n_cells
#'
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' summary(hexify(df, lon = "lon", lat = "lat", grid = grid))$column_names
NULL

#' @rdname hexify-summary
#' @export
setMethod("summary", "HexGridInfo", function(object, ...) {
  gt <- tryCatch(object@grid_type, error = function(e) "isea")

  structure(
    list(
      grid_type = gt,
      aperture = if (gt == "h3") NA_character_ else as.character(object@aperture),
      resolution = object@resolution,
      area_km2 = object@area_km2,
      diagonal_km = object@diagonal_km,
      crs = object@crs,
      radius_km = grid_radius_km(object),
      earth = is_earth_grid(object),
      orientation = grid_orientation(object),
      projection = grid_projection(object),
      polyhedron = grid_polyhedron(object),
      ellipsoid = grid_ellipsoid(object),
      n_cells = grid_n_cells(object)
    ),
    class = "hexify_grid_summary"
  )
})

#' @rdname hexify-summary
#' @export
print.hexify_grid_summary <- function(x, ...) {
  if (x$grid_type == "h3") {
    cat("HexGridInfo Specification [H3]\n")
    cat("-------------------------------\n")
    cat(sprintf("Grid Type:   H3 (Uber)\n"))
    cat(sprintf("Resolution:  %d\n", x$resolution))

    if (!is.na(x$area_km2)) {
      cat(sprintf("Avg Area:    %.4f km^2 (varies by location)\n", x$area_km2))
    }
    if (!is.na(x$diagonal_km)) {
      cat(sprintf("Avg Diagonal:%.2f km\n", x$diagonal_km))
    }
  } else {
    cat("HexGridInfo Specification\n")
    cat("-------------------------\n")
    if (!is.null(x$polyhedron) && !identical(x$polyhedron, "icosahedron")) {
      cat(sprintf("Polyhedron:  %s\n", x$polyhedron))
    }
    if (identical(x$projection, "fuller")) {
      cat("Projection:  Fuller (cells not equal-area)\n")
    } else if (identical(x$projection, "ivea")) {
      cat("Projection:  IVEA (vertex-oriented equal-area)\n")
    }
    cat(sprintf("Aperture:    %s\n", x$aperture))
    cat(sprintf("Resolution:  %d\n", x$resolution))

    if (!is.na(x$area_km2)) {
      cat(sprintf(if (!is_equal_area_projection(x$projection)) "Mean Area:   %.2f km^2\n"
                  else "Area:        %.2f km^2\n", x$area_km2))
    }
    if (!is.na(x$diagonal_km)) {
      cat(sprintf("Diagonal:    %.2f km\n", x$diagonal_km))
    }
  }

  cat(sprintf("CRS:         %s\n", format_crs(x$crs)))

  if (length(x$ellipsoid) == 2L) {
    cat(sprintf("Ellipsoid:   %s, authalic radius %.4f km\n",
                ellipsoid_label(x$ellipsoid), x$radius_km))
  } else if (!x$earth) {
    cat(sprintf("Radius:      %.2f km\n", x$radius_km))
  }

  poly <- if (is.null(x$polyhedron) || is.na(x$polyhedron)) "icosahedron" else x$polyhedron
  if (length(x$orientation) == 3L && !is_standard_orientation(x$orientation, poly)) {
    cat(sprintf("Orientation: vertex 0 at %.6f, %.6f; azimuth %.6f\n",
                x$orientation[1], x$orientation[2], x$orientation[3]))
  }

  cat(sprintf("Total Cells: %.0f\n", x$n_cells))

  if (x$grid_type == "h3") {
    cat("Note: H3 cells are NOT exactly equal-area\n")
  }

  invisible(x)
}

#' @rdname HexData-methods
#' @keywords internal
#' @export
setMethod("show", "HexData", function(object) {
  print(summary(object))
  invisible(object)
})

#' @rdname hexify-summary
#' @export
setMethod("summary", "HexData", function(object, ...) {
  preview <- NULL
  if (nrow(object@data) > 0) {
    rows <- seq_len(min(3, nrow(object@data)))
    preview <- data.frame(
      object@data[rows, seq_len(min(3, ncol(object@data))), drop = FALSE],
      cell_id = object@cell_id[rows],
      check.names = FALSE
    )
  }

  structure(
    list(
      rows = nrow(object@data),
      columns = ncol(object@data),
      column_names = names(object@data),
      n_cells = n_cells(object),
      type = if (inherits(object@data, "sf")) "sf" else "data.frame",
      grid = summary(object@grid),
      preview = preview
    ),
    class = "hexify_data_summary"
  )
})

#' @rdname hexify-summary
#' @export
print.hexify_data_summary <- function(x, ...) {
  cat("HexData Object\n")
  cat("--------------\n")
  cat(sprintf("Rows:    %d\n", x$rows))
  cat(sprintf("Columns: %d\n", x$columns))
  cat(sprintf("Cells:   %d unique\n", x$n_cells))

  if (x$type == "sf") {
    cat("Type:    sf (spatial features)\n")
  } else {
    cat("Type:    data.frame\n")
  }

  cat("\nGrid:\n")
  if (x$grid$grid_type == "h3") {
    cat(sprintf("  H3 Resolution %d", x$grid$resolution))
    if (!is.na(x$grid$area_km2)) {
      cat(sprintf(" (~%.4f km^2 avg)", x$grid$area_km2))
    }
  } else {
    cat(sprintf("  Aperture %s, Resolution %d",
                x$grid$aperture, x$grid$resolution))
    if (!is.na(x$grid$area_km2)) {
      cat(sprintf(" (~%.1f km^2)", x$grid$area_km2))
    }
  }
  cat("\n")

  cat("\nColumns: ")
  if (length(x$column_names) > 8) {
    cat(paste(x$column_names[1:8], collapse = ", "), ", ...\n")
  } else {
    cat(paste(x$column_names, collapse = ", "), "\n")
  }

  if (!is.null(x$preview)) {
    cat("\nData preview (with cell assignments):\n")
    print(x$preview, row.names = FALSE)

    if (x$rows > 3) {
      cat(sprintf("... with %d more rows\n", x$rows - 3))
    }
  }

  invisible(x)
}

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
  df$cell_area_km2 <- unname(cell_area(grid = x))
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
    grid_type = x@grid_type,
    radius_km = grid_radius_km(x),
    orientation = grid_orientation(x),
    projection = grid_projection(x),
    polyhedron = grid_polyhedron(x),
    ellipsoid = grid_ellipsoid(x)
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

  if (is_hex_grid(x)) return(upgrade_grid(x))
  if (is_hex_data(x)) return(upgrade_grid(x@grid))

  # Handle legacy hexify_grid objects (S3 class)
  if (inherits(x, "hexify_grid")) {
    return(hexify_grid_to_HexGridInfo(x))
  }

  stop("Cannot extract grid from object of class ", class(x)[1])
}

#' Cell IDs and grid of a call that takes a grid or a HexData object
#'
#' A HexData object supplies its own cells when none are given; a grid needs
#' them, unless `all` lets it supply every cell it has.
#' @param cell_id Cell IDs, or NULL to read a HexData object's own
#' @param grid A HexGridInfo, HexData or legacy hexify_grid object
#' @param all Whether a grid given without cells stands for all of its cells
#' @return List with `cell_id` and `grid`, the HexGridInfo from extract_grid()
#' @noRd
resolve_cells_grid <- function(cell_id, grid, all = FALSE) {
  g <- extract_grid(grid)
  if (is.null(cell_id)) {
    if (is_hex_data(grid)) {
      cell_id <- grid@cell_id
    } else if (all) {
      cell_id <- grid_cells(g)
    } else {
      stop("cell_id required when grid is not HexData")
    }
  }
  if (!is_h3_grid(g)) cell_id <- as_cell_id(cell_id)
  list(cell_id = cell_id, grid = g)
}

#' The cells given, or every cell of the grid
#' @param g HexGridInfo object
#' @param cells Cell IDs, or NULL for all of them
#' @noRd
grid_cells <- function(g, cells = NULL) {
  if (!is.null(cells)) return(if (is_h3_grid(g)) cells else as_cell_id(cells))
  if (is_h3_grid(g)) {
    h3_all_cells(g@resolution)
  } else {
    as_cell_id(seq_len(grid_n_cells(g)))
  }
}

#' Stop unless an object is a legacy hexify_grid
#' @param grid Object to check
#' @noRd
check_hexify_grid <- function(grid) {
  if (!inherits(grid, "hexify_grid")) {
    stop("grid must be a hexify_grid object from hexify_grid()")
  }
  invisible(grid)
}

#' Fill the slots a grid saved by an older hexify lacks
#'
#' A grid deserialized from before a slot existed reads the slot's original
#' meaning: an ISEA grid on Earth in the standard orientation, on the ISEA
#' projection.
#' @param g HexGridInfo object
#' @return HexGridInfo object
#' @noRd
upgrade_grid <- function(g) {
  if (!.hasSlot(g, "grid_type")) {
    g@grid_type <- "isea"
  }
  if (!.hasSlot(g, "radius_km")) {
    g@radius_km <- EARTH_RADIUS_KM
  }
  if (!.hasSlot(g, "orientation")) {
    g@orientation <- if (g@grid_type == "h3") numeric(0) else ISEA_ORIENTATION
  }
  if (!.hasSlot(g, "projection")) {
    g@projection <- if (g@grid_type == "h3") character(0) else "isea"
  }
  if (!.hasSlot(g, "ellipsoid")) {
    g@ellipsoid <- numeric(0)
  }
  g
}

#' Convert legacy hexify_grid to HexGridInfo
#'
#' @param x A hexify_grid object (S3)
#' @return A HexGridInfo object (S4)
#' @keywords internal
hexify_grid_to_HexGridInfo <- function(x) {
  area <- if (!is.null(x$area)) as.numeric(x$area) else NA_real_
  diagonal <- hex_spacing_km(area)

  new("HexGridInfo",
      aperture = as.character(x$aperture),
      resolution = as.integer(x$resolution),
      area_km2 = area,
      diagonal_km = diagonal,
      crs = resolve_crs(x$crs, grid_radius_km(x), grid_ellipsoid(x)),
      radius_km = grid_radius_km(x),
      orientation = grid_orientation(x),
      projection = grid_projection(x),
      polyhedron = grid_polyhedron(x),
      ellipsoid = grid_ellipsoid(x))
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

  ap <- x@aperture
  legacy_index <- index_type_for_aperture(ap)

  # Convert aperture to numeric for legacy
  aperture_num <- aperture_to_int(ap)
  orientation <- grid_orientation(x)

  grid <- list(
    area = x@area_km2,
    resolution = x@resolution,
    aperture = aperture_num,
    topology = "HEXAGON",
    projection = toupper(grid_projection(x)),
    polyhedron = grid_polyhedron(x),
    ellipsoid = grid_ellipsoid(x),
    metric = TRUE,
    radius_km = grid_radius_km(x),
    index_type = legacy_index,
    res = x@resolution,
    topology_family = "HEXAGON",
    metric_radius = if (!is.na(x@area_km2)) sqrt(x@area_km2 / pi) else NULL,
    pole_lon_deg = orientation[["vert0_lon"]],
    pole_lat_deg = orientation[["vert0_lat"]],
    azimuth_deg = orientation[["azimuth"]],
    # MIXED43 is DGGRID's own name for the 4/3 arrangement; other sequences
    # have no DGGRID aperture type, so they carry the generic label.
    aperture_type = if (ap == "4/3") "MIXED43" else if (is_mixed_aperture(ap)) "MIXED" else "SEQUENCE",
    res_spec = x@resolution,
    precision = 7
  )

  class(grid) <- c("hexify_grid", "dggs", "list")
  grid
}
