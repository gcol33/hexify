# constructors.R
# Constructor functions for HexGridInfo and HexData S4 classes
#
# These functions provide user-friendly interfaces for creating
# grid specifications and hexified data objects.

# =============================================================================
# HexGridInfo CONSTRUCTOR
# =============================================================================

#' Create a Hexagonal Grid Specification
#'
#' Creates a HexGridInfo object that stores all parameters needed for hexagonal
#' grid operations. Use this to define the grid once and pass it to all
#' downstream functions.
#'
#' @param area_km2 Target cell area in square kilometers. Mutually exclusive
#'   with \code{resolution}.
#' @param resolution Grid resolution level (0-30 for ISEA, 0-15 for H3).
#'   Mutually exclusive with \code{area_km2}. For H3, typical use cases by
#'   resolution:
#'   \itemize{
#'     \item 0-3: continental/country scale
#'     \item 4-7: regional/city scale
#'     \item 8-10: neighborhood/block scale (FCC uses 8-9)
#'     \item 11-15: building/sub-meter scale
#'   }
#' @param aperture Grid aperture: 3 (default), 4, 7, a mixed family such as
#'   "4/3", "4/7" or "7/4", or one aperture per resolution level as a vector,
#'   e.g. \code{c(4, 4, 7, 3)}. A family name refines by the first aperture for
#'   the first \code{floor(resolution / 2)} levels and by the second for the
#'   rest, which is how DGGRID arranges ISEA43H. A per-level vector needs
#'   \code{resolution} rather than \code{area_km2}, and is stored as a
#'   comma-separated string (\code{"4,4,7,3"}), which \code{aperture} also
#'   accepts. It names no level finer than \code{resolution}, so
#'   \code{get_children()} needs a family name. Ignored for H3 grids (fixed
#'   at 7).
#' @param type Grid type: "isea" (default) or "h3".
#' @param resround Resolution rounding when using \code{area_km2}:
#'   "nearest" (default), "up", or "down".
#' @param crs Coordinate reference system: an EPSG code, or a 'PROJ' or 'WKT'
#'   string. Defaults to 'WGS84' on Earth, and to a longlat CRS on the sphere of
#'   \code{radius_km} on any other body, which has no EPSG code to name it.
#' @param radius_km Radius of the body the grid covers, in kilometers, or the
#'   name of a body: "mercury", "venus", "earth" (default), "moon", "mars",
#'   "ceres", "jupiter", "io", "europa", "ganymede", "callisto", "saturn",
#'   "enceladus", "titan", "uranus", "neptune", "pluto".
#' @param orientation Where the icosahedron of an ISEA grid sits on the
#'   sphere: "standard" (default), "random", "region", "face", or
#'   \code{c(vert0_lon, vert0_lat, azimuth)} in degrees. See the Orientation
#'   section. H3 fixes its own orientation, so H3 grids take only "standard".
#' @param projection How an ISEA-family grid projects each icosahedron face
#'   onto its plane triangle: "isea" (default), Snyder's equal-area
#'   projection, or "fuller", Fuller's projection. See the Projection section.
#'   H3 fixes its own projection, so H3 grids take only "isea".
#' @param region For \code{orientation = "region"} or \code{"face"}, the area
#'   to centre the grid
#'   on: \code{c(lon, lat)} in degrees, or an sf, sfc or bbox object, whose
#'   spherical centroid is used.
#' @param polyhedron The solid an ISEA-family grid is built on:
#'   "icosahedron" (default) or "octahedron". See the Polyhedron section. H3
#'   grids are built on the icosahedron and take only the default.
#'
#' @return A HexGridInfo object containing the grid specification.
#'
#' @details
#' Exactly one of \code{area_km2} or \code{resolution} must be provided.
#'
#' When \code{area_km2} is provided, the resolution is calculated automatically
#' using the cell count formula: N = 10 * aperture^res + 2 (ISEA on the
#' icosahedron, 4 * aperture^res + 2 on the octahedron) or by matching the
#' closest H3 resolution.
#'
#' H3 grids use the Uber H3 hierarchical hexagonal system. Unlike ISEA grids,
#' H3 cells are NOT exactly equal-area: hexagon area varies by about a factor
#' of 2 within a resolution, and by about 2.4 once the twelve pentagons are
#' counted. The variation follows position on the icosahedron rather than
#' latitude, with the smallest cells near face centres.
#'
#' @section Other Bodies:
#'
#' A grid is a partition of the sphere, and \code{radius_km} sets the sphere it
#' is measured on. Cell geometry -- which cell a coordinate lands in, where cell
#' centres and corners sit, the hierarchy, the neighbours -- is angular and
#' identical on every body; the radius sets the kilometer figures: cell area,
#' diagonal, spacing, and the resolution that \code{area_km2} picks. Earth's
#' area comes from the 'WGS84' ellipsoid, every other radius gives the sphere
#' area 4*pi*r^2.
#'
#' \preformatted{
#' mars <- hex_grid(area_km2 = 1000, radius_km = "mars")
#' hex_grid(resolution = 8, radius_km = 3389.5)   # the same grid
#' }
#'
#' Both backends take a radius. 'H3' reports a cell's area as its solid angle
#' times Earth's radius squared, so another radius scales those areas by the
#' square of the radius ratio, exactly. One caveat carries: an 'H3' cell ID
#' names a position in 'H3''s topology, which 'Uber''s 'H3' reads on Earth, so
#' the IDs of a grid on another body are that topology on that body and are not
#' interchangeable with Earth 'H3' data.
#'
#' @section Orientation:
#'
#' An ISEA grid is built on a solid, an icosahedron unless \code{polyhedron}
#' names another, and its orientation places that solid on the sphere: vertex
#' 0 at longitude \code{vert0_lon} and latitude \code{vert0_lat}, and vertex 1
#' at azimuth \code{azimuth} seen from vertex 0, all in degrees, as DGGRID's
#' \code{dggs_vert0_lon}, \code{dggs_vert0_lat} and \code{dggs_vert0_azimuth}
#' do.
#'
#' \itemize{
#'   \item "standard" is the ISEA orientation on the icosahedron: vertex 0 at
#'     11.25E, 58.28N, azimuth 0, which places the twelve pentagons over the
#'     oceans. On the octahedron it puts vertices at both poles and at
#'     longitudes 0, 90E, 180 and 90W on the equator.
#'   \item "random" draws vertex 0 uniformly on the sphere and the azimuth
#'     uniformly in [0, 360), from R's random number generator, so
#'     \code{set.seed()} repeats it. Grids in several random orientations show
#'     how much a result depends on where the cell boundaries fall.
#'   \item "region" places the grid as DGGRID's \code{REGION_CENTER} does: the
#'     centre of \code{region} lands on the midpoint of an edge of the solid,
#'     the middle of the two faces sharing it, far from the vertices where the
#'     vertex cells and the largest distortion sit.
#'   \item "face" places the centre of \code{region} on the centre of a face.
#'     On the icosahedron that is about 20.9 degrees from the nearest face edge
#'     and 37.4 degrees from the nearest vertex, so a region up to about 40
#'     degrees across lies on one face, away from the face edges where the
#'     projection bends.
#' }
#'
#' Rotating the icosahedron rotates the grid with it, so cell IDs, the cell
#' hierarchy and neighbours are the same under every orientation; only where
#' each cell sits on the sphere changes. Every function taking the grid reads
#' its orientation.
#'
#' \preformatted{
#' alps <- hex_grid(area_km2 = 100, orientation = "region", region = c(10, 46.5))
#' hex_grid(resolution = 8, orientation = c(0, 90, 0))   # vertex 0 at the pole
#' }
#'
#' @section Projection:
#'
#' An ISEA-family grid lays its cells out on the plane triangles of the
#' icosahedron's faces, and the projection carries each spherical face onto its
#' triangle. With \code{"isea"}, Snyder's icosahedral equal-area projection,
#' every cell of a resolution has the same area. With \code{"fuller"}, Fuller's
#' projection keeps lengths along the face edges and is not equal-area, so
#' \code{area_km2} is the mean cell area and \code{cell_area()} reports each
#' cell's own. These are DGGRID's \code{dggs_proj ISEA} and \code{FULLER},
#' giving FULLER3H, FULLER4H, FULLER7H and FULLER43H. Cell IDs, the hierarchy
#' and neighbours are the same under both projections; where each cell's centre
#' and corners sit on the sphere differs.
#'
#' \preformatted{
#' hex_grid(resolution = 5, aperture = 4, projection = "fuller")   # FULLER4H
#' }
#'
#' @section Polyhedron:
#'
#' Snyder's equal-area projection is defined on every regular solid with
#' triangular faces, and a hexagonal grid lies on the icosahedron and on the
#' octahedron alike: each pair of faces sharing an edge forms a diamond of
#' cells, and the cell at each vertex has one side per face meeting there.
#' The icosahedron has twelve pentagons and the octahedron six squares; every
#' other cell is a hexagon, and every cell of a resolution has the same area
#' (the vertex cells two thirds or five sixths of it, by their sides). The
#' octahedron's cells are more distorted than the icosahedron's, so it suits
#' layouts that need its four-fold symmetry rather than analysis. Fuller's
#' projection is defined on the icosahedron only. The tetrahedron carries the
#' face projection (see \code{\link{hexify_forward}}) but no grid: its faces
#' do not pair into diamonds.
#'
#' \preformatted{
#' hex_grid(resolution = 5, aperture = 4, polyhedron = "octahedron")
#' }
#'
#' @seealso \code{\link{hexify}} for assigning points to cells,
#'   \code{\link{HexGridInfo-class}} for class documentation
#'
#' @section One Grid, Many Datasets:
#'
#' A HexGridInfo acts as a shared spatial reference system - like a CRS, but
#' discrete and equal-area. Define the grid once, then attach multiple
#' datasets without repeating parameters:
#'
#' \preformatted{
#' # Step 1: Define the grid once
#' grid <- hex_grid(area_km2 = 1000)
#'
#' # Step 2: Attach multiple datasets to the same grid
#' birds <- hexify(bird_obs, lon = "longitude", lat = "latitude", grid = grid)
#' mammals <- hexify(mammal_obs, lon = "lon", lat = "lat", grid = grid)
#' climate <- hexify(weather_stations, lon = "x", lat = "y", grid = grid)
#'
#' # No aperture, resolution, or area needed after step 1 - the grid
#' # travels with the data.
#'
#' # Step 3: Work at the cell level
#' # Once hexified, lon/lat no longer matter - cell_id is the shared key
#' bird_counts <- aggregate(species ~ cell_id, data = as.data.frame(birds), length)
#' mammal_richness <- aggregate(species ~ cell_id, data = as.data.frame(mammals),
#'                              function(x) length(unique(x)))
#'
#' # Join datasets by cell_id - guaranteed to align because same grid
#' combined <- merge(bird_counts, mammal_richness, by = "cell_id")
#'
#' # Step 4: Visual confirmation
#' # All datasets produce identical grid overlays
#' plot(birds)   # See the grid
#' plot(mammals) # Same grid, different data
#' }
#'
#' @export
#' @examples
#' # Create grid by target area
#' grid <- hex_grid(area_km2 = 1000)
#' print(grid)
#'
#' # Create grid by resolution
#' grid <- hex_grid(resolution = 8, aperture = 3)
#'
#' # Create grid with different aperture
#' grid4 <- hex_grid(area_km2 = 500, aperture = 4)
#'
#' # Create mixed aperture grid
#' grid43 <- hex_grid(area_km2 = 1000, aperture = "4/3")
#'
#' # Mix in aperture 7, either as a family or level by level
#' grid47 <- hex_grid(area_km2 = 1000, aperture = "4/7")
#' grid_seq <- hex_grid(resolution = 4, aperture = c(4, 4, 7, 3))
#'
#' # Grid on another body, by name or by radius
#' mars <- hex_grid(area_km2 = 1000, radius_km = "mars")
#' titan <- hex_grid(resolution = 6, radius_km = 2574.76)
#'
#' # Use grid in hexify
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
hex_grid <- function(area_km2 = NULL,
                     resolution = NULL,
                     aperture = 3,
                     type = c("isea", "h3"),
                     resround = "nearest",
                     crs = NULL,
                     radius_km = EARTH_RADIUS_KM,
                     orientation = "standard",
                     projection = c("isea", "fuller"),
                     region = NULL,
                     polyhedron = c("icosahedron", "octahedron", "tetrahedron")) {

  type <- match.arg(type)
  projection <- match.arg(projection)
  polyhedron <- match.arg(polyhedron)

  radius_km <- resolve_radius_km(radius_km)
  crs <- resolve_crs(crs, radius_km)

  # =========================================================================
  # H3 grid path
  # =========================================================================
  if (type == "h3") {
    if (!missing(aperture) && aperture != 3) {
      warning("aperture is ignored for H3 grids (H3 uses fixed aperture 7)")
    }
    if (!identical(orientation, "standard") || !is.null(region)) {
      stop("H3 fixes its own orientation; orientation applies to ISEA grids",
           call. = FALSE)
    }
    if (projection != "isea") {
      stop("H3 fixes its own projection; projection applies to ISEA grids",
           call. = FALSE)
    }
    if (polyhedron != "icosahedron") {
      stop("H3 is built on the icosahedron; polyhedron applies to ISEA grids",
           call. = FALSE)
    }

    # Validate: exactly one of area_km2 or resolution
    if (is.null(area_km2) && is.null(resolution)) {
      stop("Exactly one of 'area_km2' or 'resolution' must be provided")
    }
    if (!is.null(area_km2) && !is.null(resolution)) {
      stop("Provide either 'area_km2' or 'resolution', not both")
    }

    if (!is.null(area_km2)) {
      if (!is.numeric(area_km2) || area_km2 <= 0) {
        stop("area_km2 must be a positive number")
      }
      # Find closest H3 resolution by area
      resolution <- closest_h3_resolution(area_km2, radius_km)
      warning(sprintf(
        "H3 cells are not exactly equal-area. Closest resolution %d has average area ~%.3f km^2 (requested %.3f km^2)",
        resolution, h3_avg_area_km2(resolution, radius_km), area_km2
      ))
    } else {
      if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution)) {
        stop("resolution must be a single non-NA number")
      }
      resolution <- as.integer(resolution)
      if (resolution < H3_MIN_RESOLUTION || resolution > H3_MAX_RESOLUTION) {
        stop(sprintf("H3 resolution must be between %d and %d",
                     H3_MIN_RESOLUTION, H3_MAX_RESOLUTION))
      }
      rlang::inform(
        "H3 cells are not exactly equal-area; hexagon area varies ~2x within a resolution.",
        .frequency = "once",
        .frequency_id = "hexify_h3_not_equal_area"
      )
    }

    if (radius_km != EARTH_RADIUS_KM) {
      rlang::inform(
        paste0(
          "H3 cell IDs name a position in H3's topology, which Uber's H3 reads ",
          "on Earth. A grid on another body reuses that topology and its own ",
          "radius for areas; the IDs are not interchangeable with Earth H3 data."
        ),
        .frequency = "once",
        .frequency_id = "hexify_h3_other_body"
      )
    }

    actual_area <- h3_avg_area_km2(resolution, radius_km)
    actual_diagonal <- hex_spacing_km(actual_area)

    grid <- new("HexGridInfo",
                aperture = "7",
                resolution = as.integer(resolution),
                area_km2 = as.numeric(actual_area),
                diagonal_km = as.numeric(actual_diagonal),
                crs = crs,
                grid_type = "h3",
                radius_km = radius_km,
                orientation = numeric(0),
                projection = character(0),
                polyhedron = character(0))
    return(grid)
  }

  # =========================================================================
  # ISEA grid path
  # =========================================================================

  if (!polyhedron %in% GRID_POLYHEDRA) {
    stop("the ", polyhedron, " carries no hexagonal grid: a vertex shared by ",
         "three faces would start a diamond of cells with all three, so every ",
         "such vertex must be one of the two single-cell quads, and the ",
         "tetrahedron has four. Its face projection is available through ",
         "hexify_forward() and hexify_inverse().", call. = FALSE)
  }
  if (projection == "fuller" && polyhedron != "icosahedron") {
    stop("Fuller's projection is defined on the icosahedron only", call. = FALSE)
  }

  # -------------------------------------------------------------------------
  # Parse aperture: a single aperture, a family such as "4/3", or one aperture
  # per resolution level (see R/aperture_sequence.R)
  # -------------------------------------------------------------------------
  if (length(aperture) > 1L && is.null(resolution)) {
    stop("An aperture given per level needs 'resolution', not 'area_km2'")
  }
  aperture_str <- format_aperture(aperture, resolution)

  if (!is_mixed_aperture(aperture_str) &&
      !aperture_str %in% as.character(VALID_APERTURES)) {
    stop("Aperture must be 3, 4, 7, a family such as \"4/3\", or one aperture per level")
  }

  # -------------------------------------------------------------------------
  # Validate area_km2 / resolution (exactly one required)
  # -------------------------------------------------------------------------
  if (is.null(area_km2) && is.null(resolution)) {
    stop("Exactly one of 'area_km2' or 'resolution' must be provided")
  }
  if (!is.null(area_km2) && !is.null(resolution)) {
    stop("Provide either 'area_km2' or 'resolution', not both")
  }

  # -------------------------------------------------------------------------
  # Calculate resolution from area if needed
  # -------------------------------------------------------------------------
  if (!is.null(area_km2)) {
    if (!is.numeric(area_km2) || area_km2 <= 0) {
      stop("area_km2 must be a positive number")
    }
    resolution <- resolve_resolution_from_area(area_km2, aperture_str,
                                               radius_km, resround, polyhedron)
  } else {
    if (!is.numeric(resolution) || length(resolution) != 1 || is.na(resolution)) {
      stop("resolution must be a single non-NA number")
    }
    resolution <- as.integer(resolution)
    if (resolution < MIN_RESOLUTION || resolution > MAX_RESOLUTION) {
      stop(sprintf("Resolution must be between %d and %d",
                   MIN_RESOLUTION, MAX_RESOLUTION))
    }
  }

  # -------------------------------------------------------------------------
  # Calculate actual area and diagonal for this resolution
  # -------------------------------------------------------------------------
  actual_area <- mean_cell_area_km2(aperture_str, resolution, radius_km,
                                    polyhedron)
  actual_diagonal <- hex_spacing_km(actual_area)

  # -------------------------------------------------------------------------
  # Create and validate HexGridInfo object
  # -------------------------------------------------------------------------
  grid <- new("HexGridInfo",
              aperture = aperture_str,
              resolution = as.integer(resolution),
              area_km2 = as.numeric(actual_area),
              diagonal_km = as.numeric(actual_diagonal),
              crs = crs,
              grid_type = "isea",
              radius_km = radius_km,
              orientation = resolve_orientation(orientation, region, polyhedron),
              projection = projection,
              polyhedron = polyhedron)

  # Validation happens automatically via setValidity
  grid
}

# =============================================================================
# HexData CONSTRUCTOR (Internal)
# =============================================================================

#' Create a HexData Object (Internal)
#'
#' Internal constructor for HexData objects. Users should use \code{hexify()}
#' instead.
#'
#' @param data Data frame or sf object (original user data, untouched)
#' @param grid HexGridInfo object
#' @param cell_id Numeric vector of cell IDs for each row
#' @param cell_center Matrix with columns lon, lat for cell centers
#'
#' @return A HexData object
#' @keywords internal
new_hex_data <- function(data,
                         grid,
                         cell_id,
                         cell_center) {

  # Validate inputs
  if (!inherits(data, "data.frame") && !inherits(data, "sf")) {
    stop("data must be a data.frame or sf object")
  }

  if (!is_hex_grid(grid)) {
    stop("grid must be a HexGridInfo object")
  }

  # Ensure cell_center is a matrix with correct column names
  if (!is.matrix(cell_center)) {
    cell_center <- as.matrix(cell_center)
  }
  if (is.null(colnames(cell_center))) {
    colnames(cell_center) <- c("lon", "lat")
  }

  # For ISEA grids, coerce to numeric; for H3, keep as character
  if (is_h3_grid(grid)) {
    cell_id <- as.character(cell_id)
  } else {
    cell_id <- as.numeric(cell_id)
  }

  new("HexData",
      data = data,
      grid = grid,
      cell_id = cell_id,
      cell_center = cell_center)
}

# =============================================================================
# TIBBLE COERCION (if tibble is available)
# =============================================================================

#' Convert HexData to tibble
#'
#' @param x A HexData object
#' @param ... Additional arguments (ignored)
#' @return A tibble
#'
#' @exportS3Method tibble::as_tibble
as_tibble.HexData <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for as_tibble(). ",
         "Install with: install.packages('tibble')")
  }

  df <- x@data
  if (inherits(df, "sf")) {
    df <- as.data.frame(sf::st_drop_geometry(df))
  }

  tibble::as_tibble(df)
}

# =============================================================================
# SF COERCION
# =============================================================================

#' Convert gridded data or a grid to sf
#'
#' Methods for \code{sf::st_as_sf()}. Gridded data comes back as one feature per
#' row, carrying either the centre of the cell each row fell in or that cell's
#' boundary. A grid specification comes back as its global cell set.
#'
#' @param x A HexData or HexGridInfo object
#' @param ... For a HexGridInfo, passed on to \code{\link{grid_global}}
#' @param geometry Type of geometry: "point" (default) or "polygon"
#' @param shape,depth Cell shape for polygon geometry, as for
#'   \code{\link{cell_to_sf}}
#'
#' @return An sf object
#'
#' @details
#' For point geometry, cell centers (cell_cen_lon, cell_cen_lat) are used.
#' For polygon geometry, cell boundaries are computed using the grid specification.
#'
#' @seealso \code{\link{cell_to_sf}} for the cells of named cell IDs,
#'   \code{\link{grid_global}} for every cell of a grid
#'
#' @export
#' @examples
#' df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Get sf points
#' sf_pts <- st_as_sf(result)
#'
#' # Get sf polygons
#' sf_poly <- st_as_sf(result, geometry = "polygon")
#'
#' # Every cell of a coarse grid
#' cells <- st_as_sf(hex_grid(resolution = 2))
st_as_sf.HexData <- function(x, ..., geometry = c("point", "polygon"),
                             shape = c("hexagon", "gosper", "descendants"),
                             depth = 3L) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  geometry <- match.arg(geometry)
  data <- x@data
  grid <- x@grid

  # If already sf, check if we need to change geometry type
  if (inherits(data, "sf")) {
    if (geometry == "point") {
      return(data)  # Already has geometry
    }
    # For polygon, we need to replace geometry
    data <- sf::st_drop_geometry(data)
  }

  if (geometry == "point") {
    # Point geometry from cell centers (stored in cell_center slot)
    df_with_coords <- cbind(
      data,
      cell_id = x@cell_id,
      cell_cen_lon = x@cell_center[, "lon"],
      cell_cen_lat = x@cell_center[, "lat"]
    )
    sf::st_as_sf(
      df_with_coords,
      coords = c("cell_cen_lon", "cell_cen_lat"),
      crs = grid_crs(grid)
    )

  } else {
    # Polygon geometry from cell boundaries
    unique_ids <- unique(x@cell_id)

    polys_sf <- cell_to_sf(unique_ids, grid, shape = match.arg(shape),
                           depth = depth)

    # Add cell_id to data for merge
    data_with_id <- cbind(data, cell_id = x@cell_id)

    # Merge with original data
    result <- merge(polys_sf, data_with_id, by = "cell_id", all.y = TRUE)
    sf::st_as_sf(result)
  }
}

#' @rdname st_as_sf.HexData
#' @export
st_as_sf.HexGridInfo <- function(x, ...) {
  grid_global(x, ...)
}

#' Coordinate reference system of a grid
#'
#' Reads the CRS a grid's coordinates are in, as \code{sf::st_crs()} does for
#' any spatial object. An Earth grid returns 'WGS84'; a grid built on another
#' body returns a longlat CRS on the sphere of its radius.
#'
#' @param x A HexGridInfo or HexData object
#' @param ... Passed on to \code{sf::st_crs()}
#' @return An object of class \code{crs}, as \code{sf::st_crs()} returns:
#'   a list carrying the reference system in 'PROJ' and 'WKT' form. It says how
#'   to read the coordinates that the grid's cells, centres and sf exports come
#'   back in.
#'
#' @importFrom sf st_crs
#' @export
#' @examples
#' st_crs(hex_grid(resolution = 5))
#' st_crs(hex_grid(resolution = 5, radius_km = "mars"))
st_crs.HexGridInfo <- function(x, ...) {
  grid_crs(x)
}

#' @rdname st_crs.HexGridInfo
#' @export
st_crs.HexData <- function(x, ...) {
  grid_crs(x@grid)
}
