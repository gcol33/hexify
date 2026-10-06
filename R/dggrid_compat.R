# R/dggrid_compat.R
# dggridR compatibility layer
#
# Functions for converting between hexify and dggridR grid objects.
# These enable interoperability with dggridR without copying its API.

# =============================================================================
# GRID OBJECT CONVERSION
# =============================================================================

#' Convert hexify grid to 'dggridR'-compatible grid object
#'
#' Creates a 'dggridR'-compatible grid specification from a hexify_grid object.
#' The resulting object can be used with 'dggridR' functions that accept a dggs
#' object.
#'
#' @param grid A hexify_grid object from hexify_grid()
#'
#' @details
#' A dggs carries no body radius, so a grid built on another body cannot be
#' expressed as one: 'dggridR' reads every dggs on Earth's radius, and its areas,
#' spacings and CLS come back on Earth. Converting such a grid warns.
#'
#' @return A list with 'dggridR'-compatible fields:
#'   \item{pole_lon_deg}{Longitude of icosahedron vertex 0 (standard 11.25)}
#'   \item{pole_lat_deg}{Latitude of icosahedron vertex 0 (standard
#'     58.282525588538995)}
#'   \item{azimuth_deg}{Azimuth of vertex 1 seen from vertex 0 (standard 0)}
#'   \item{aperture}{Grid aperture (3, 4, or 7)}
#'   \item{res}{Resolution level}
#'   \item{topology}{Grid topology ("HEXAGON")}
#'   \item{projection}{Face projection ('ISEA' or 'FULLER')}
#'   \item{precision}{Output decimal precision (default 7)}
#'
#' @family 'dggridR' compatibility
#' @export
as_dggrid <- function(grid) {

  # Accept the modern HexGridInfo (S4) object by converting to the legacy
  # list representation this function builds from, matching dgverify()'s
  # fix (#7) for the same gap.
  if (is_hex_grid(grid)) {
    if (grid@grid_type == "h3") {
      stop("as_dggrid() has no dggridR-compatible representation for H3 grids")
    }
    grid <- HexGridInfo_to_hexify_grid(grid)
  }

  check_hexify_grid(grid)

  if (!is_earth_grid(grid)) {
    warning("A dggs carries no body radius: this grid is sized on a radius of ",
            format(grid_radius_km(grid)), " km, and 'dggridR' will read it on ",
            "Earth's ", format(EARTH_RADIUS_KM), " km.", call. = FALSE)
  }

  orientation <- grid_orientation(grid)
  dggs <- list(
    pole_lon_deg = orientation[["vert0_lon"]],
    pole_lat_deg = orientation[["vert0_lat"]],
    azimuth_deg = orientation[["azimuth"]],
    aperture = grid$aperture,
    res = grid$resolution,
    topology = "HEXAGON",
    projection = toupper(grid_projection(grid)),
    precision = 7L
  )

  class(dggs) <- c("dggs", "list")
  dggs
}

#' Convert 'dggridR' grid object to hexify_grid
#'
#' Creates a hexify_grid object from a 'dggridR' dggs object. This allows
#' using hexify functions with grids created by 'dggridR' dgconstruct().
#'
#' @param dggs A 'dggridR' grid object from dgconstruct()
#' @param radius_km Radius the grid is sized on, in km or as a body name. A dggs
#'   carries no radius, so the default is Earth's; pass this when the dggs
#'   describes a grid on another body.
#'
#' @return A hexify_grid object
#'
#' @details
#' The 'ISEA' and 'FULLER' projections with HEXAGON topology are supported.
#' Other configurations will generate warnings.
#'
#' A dggs has no field for the body it is sized on, so the resolution is all
#' that carries over and the result is an Earth grid unless \code{radius_km}
#' says otherwise. Cell area follows from the radius, so it is read from the
#' radius given rather than from the dggs.
#'
#' The orientation carries over: \code{pole_lon_deg}, \code{pole_lat_deg} and
#' \code{azimuth_deg} place the icosahedron as DGGRID's \code{dggs_vert0_lon},
#' \code{dggs_vert0_lat} and \code{dggs_vert0_azimuth} do, and a field left out
#' takes its standard ISEA value.
#'
#' The function validates that the 'dggridR' grid uses compatible settings:
#' - Projection must be 'ISEA' or 'FULLER'
#' - Topology must be "HEXAGON" (DIAMOND, TRIANGLE not supported)
#' - Aperture must be 3, 4, or 7
#'
#' @family 'dggridR' compatibility
#' @export
from_dggrid <- function(dggs, radius_km = EARTH_RADIUS_KM) {
  # Validate dggridR object

  if (!is.list(dggs)) {
    stop("dggs must be a list (dggridR grid object)")
  }

  required <- c("res", "aperture", "topology", "projection")
  missing <- setdiff(required, names(dggs))
  if (length(missing) > 0) {
    stop(sprintf("dggs missing required fields: %s", paste(missing, collapse = ", ")))
  }

  # Check compatibility

  projection <- if (is.null(dggs$projection)) "ISEA" else dggs$projection
  if (!projection %in% DGGS_PROJECTIONS) {
    warning("Only the ISEA and FULLER projections are supported. ",
            "Results may differ from dggridR.")
    projection <- "ISEA"
  }

  if (!is.null(dggs$topology) && dggs$topology != "HEXAGON") {
    warning("Only HEXAGON topology is supported. Results may differ from dggridR.")
  }

  if (!dggs$aperture %in% c(3L, 4L, 7L)) {
    stop(sprintf("Aperture %d not supported. Must be 3, 4, or 7.", dggs$aperture))
  }

  orientation <- dggs_orientation(dggs)

  # Create hexify_grid. `area` is a throwaway placeholder -- both `resolution`
  # and `area` are overwritten below from `dggs$res` directly, but
  # hexify_grid() requires a valid positive area (see #43) to construct the
  # object at all.
  hexify_grid(
    area = 1000,
    topology = "HEXAGON",
    metric = TRUE,
    resround = "nearest",
    aperture = as.integer(dggs$aperture),
    projection = projection,
    radius_km = radius_km
  ) -> grid

  # Override resolution to match dggridR exactly
  validate_resolution(dggs$res)
  grid$resolution <- as.integer(dggs$res)
  grid$res <- as.integer(dggs$res)
  grid$pole_lon_deg <- orientation[["vert0_lon"]]
  grid$pole_lat_deg <- orientation[["vert0_lat"]]
  grid$azimuth_deg <- orientation[["azimuth"]]

  # Calculate actual area for this resolution
  n_cells <- max_cell_id(grid$resolution, grid$aperture, grid_polyhedron(grid))
  grid$area <- body_surface_km2(grid_radius_km(grid)) / n_cells

  grid
}

#' Validate 'dggridR' grid compatibility with hexify
#'
#' Checks whether a 'dggridR' grid object is compatible with hexify functions.
#' Returns TRUE if compatible, or throws an error describing incompatibilities.
#'
#' @param dggs A 'dggridR' grid object
#' @param strict If TRUE (default), throw errors for incompatibilities.
#'   If FALSE, return FALSE instead of throwing errors.
#'
#' @return TRUE if compatible, FALSE if not compatible (when strict=FALSE)
#'
#' @family 'dggridR' compatibility
#' @export
dggrid_is_compatible <- function(dggs, strict = TRUE) {
  issues <- character()

  if (!is.list(dggs)) {
    issues <- c(issues, "dggs must be a list")
  } else {
    if (is.null(dggs$projection) || !dggs$projection %in% DGGS_PROJECTIONS) {
      issues <- c(issues, "Projection must be ISEA or FULLER")
    }

    if (is.null(dggs$topology) || dggs$topology != "HEXAGON") {
      issues <- c(issues, "Only HEXAGON topology supported (not DIAMOND/TRIANGLE)")
    }

    if (is.null(dggs$aperture) || !dggs$aperture %in% c(3L, 4L, 7L)) {
      issues <- c(issues, "Aperture must be 3, 4, or 7")
    }

    orientation_ok <- tryCatch({
      dggs_orientation(dggs)
      TRUE
    }, error = function(e) FALSE)
    if (!orientation_ok) {
      issues <- c(issues, paste0("pole_lon_deg, pole_lat_deg and azimuth_deg must ",
                                 "be finite, pole_lat_deg in [-90, 90]"))
    }
  }

  if (length(issues) > 0) {
    if (strict) {
      stop(sprintf("dggridR grid not compatible with hexify:\n- %s",
                   paste(issues, collapse = "\n- ")))
    }
    return(FALSE)
  }

  TRUE
}

#' Orientation a dggs carries
#'
#' The DGGRID fields dggs_vert0_lon, dggs_vert0_lat and dggs_vert0_azimuth,
#' which 'dggridR' names pole_lon_deg, pole_lat_deg and azimuth_deg. A field
#' left out takes its standard ISEA value. DGGRID and 'dggridR' write the
#' standard vertex 0 latitude to eight decimals (58.28252559), so an
#' orientation within 1e-6 degrees of the standard one is read as it.
#' @noRd
dggs_orientation <- function(dggs) {
  field <- function(name, standard) if (is.null(dggs[[name]])) standard else dggs[[name]]
  o <- check_orientation(c(
    vert0_lon = field("pole_lon_deg", ISEA_VERT0_LON_DEG),
    vert0_lat = field("pole_lat_deg", ISEA_VERT0_LAT_DEG),
    azimuth = field("azimuth_deg", ISEA_AZIMUTH_DEG)
  ))
  if (all(abs(o - ISEA_ORIENTATION) <= 1e-6)) ISEA_ORIENTATION else o
}
