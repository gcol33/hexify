# DGGRS definitions in the JSON schema of OGC API - DGGS

#' Definition of a grid as a Discrete Global Grid Reference System
#'
#' Describes a grid in the DGGRS definition schema of OGC API - Discrete
#' Global Grid Systems: the hierarchy of grids (solid, refinement ratio,
#' refinement strategy, zone shapes, orientation), the zone identifiers and
#' the order of sub-zones.
#'
#' @param grid A HexGridInfo or HexData object.
#'
#' @return A named list with the schema's fields \code{title},
#'   \code{description}, \code{dggh}, \code{zirs} and \code{subZoneOrder}.
#'   \code{jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE)} writes it as
#'   JSON.
#'
#' @details
#' The definition describes the grid as hexify builds it. The zone
#' identifiers are the strings \code{\link{cell_to_index}} writes, which carry
#' their resolution and are unique across all resolutions, and the cell IDs,
#' which number the cells of one resolution from 1 and are unique only
#' together with it.
#'
#' The ISEA3H, ISEA7H, IVEA3H and IVEA7H definitions registered with OGC
#' (\url{https://www.opengis.net/def/dggrs/OGC/1.0/ISEA7H}) convert WGS84
#' geodetic latitude to authalic latitude before projecting, and place the
#' icosahedron's first vertex at 11.20 degrees E and authalic latitude
#' \eqn{\arctan(\phi)}: \code{hex_grid(aperture = 3, ellipsoid = "WGS84",
#' orientation = "ogc")}, its aperture-7 twin and the same with
#' \code{projection = "ivea"}. For such a grid the definition carries OGC's
#' zone identifiers (\code{cell_to_index(form = "textZIRS")}) and links to the
#' OGC definition, which it matches in its hierarchy and zone identifiers; it
#' carries no OGC URI, since hexify lists sub-zones in ascending order of
#' cell ID rather than in OGC's scanlines. On a grid with an ellipsoid the
#' orientation's latitude is geodetic, as OGC states it.
#'
#' @references
#' Open Geospatial Consortium. OGC API - Discrete Global Grid Systems - Part
#' 1: Core. OGC 21-038r1. \url{https://docs.ogc.org/is/21-038r1/21-038r1.html}
#'
#' @seealso \code{\link{hex_grid}}, \code{\link{cell_to_index}}
#'
#' @export
#' @examples
#' def <- dggrs_definition(hex_grid(resolution = 5, aperture = 7))
#' def$dggh$definition$refinementRatio
#' def$zirs$textZIRS$description
dggrs_definition <- function(grid) {
  g <- extract_grid(grid)
  if (is_h3_grid(g)) return(h3_dggrs_definition(g))

  steps <- if (is_mixed_aperture(g@aperture)) {
    isea_levels(g@aperture, g@resolution)$ap_seq[-1]
  } else {
    as.integer(g@aperture)
  }
  apertures <- unique(steps)
  projection <- grid_projection(g)
  polyhedron <- grid_polyhedron(g)
  o <- grid_orientation(g)

  strategy <- c("centredChildCell",
                c(`3` = "nodeCentredChildCell", `4` = "edgeCentredChildCell",
                  `7` = "nodeSharingChildCell")[as.character(apertures)])
  vertex_cell <- if (polyhedron == "icosahedron") "pentagon" else "square"
  e <- grid_ellipsoid(g)
  ogc <- ogc_dggrs_name(g)
  earth_model <- if (length(e) == 2L) {
    sprintf(paste("on the %s ellipsoid, geodetic latitude converted to authalic",
                  "latitude on the sphere of the ellipsoid's area, radius %s km."),
            ellipsoid_label(e), format(grid_radius_km(g), digits = 10))
  } else {
    sprintf("on a sphere of radius %s km taking latitude as spherical latitude.",
            format(grid_radius_km(g), digits = 10))
  }
  parameters <- list(sphere = list(radius_km = grid_radius_km(g)))
  if (length(e) == 2L) {
    parameters <- c(list(ellipsoid = if (identical(unname(e), unname(ELLIPSOIDS$wgs84))) {
      "[EPSG:7030]"
    } else {
      list(semiMajorAxis_km = e[["a_km"]], flattening = e[["f"]])
    }), parameters)
  }
  parameters$orientation <- list(
    latitude = geodetic_lat(unname(o[["vert0_lat"]]), icosa_arg(g)),
    longitude = unname(o[["vert0_lon"]]),
    azimuth = unname(o[["azimuth"]]),
    description = if (length(e) == 2L) {
      sprintf(paste("Geodetic latitude and longitude of the solid's first vertex,",
                    "at authalic latitude %s, and azimuth of its second vertex seen",
                    "from the first."), format(unname(o[["vert0_lat"]]), digits = 15))
    } else {
      "Spherical latitude and longitude of the solid's first vertex, and azimuth of its second vertex seen from the first."
    })

  zirs <- if (!is.null(ogc)) {
    list(textZIRS = list(description = paste(
      "OGC's", ogc, "textZIRS, as DGGAL writes it: a letter for the level,",
      "the root rhombus (0-9, A and B for the poles), the sub-rhombus counted",
      "row by row in hexadecimal, and a letter for the zone at that",
      "sub-rhombus, e.g. E6-317-A; cell_to_index(form = \"textZIRS\")."),
      type = "levelRootFaceHexRowMajorSubZone"),
      uint64ZIRS = list(description = paste(
        "OGC's", ogc, "uint64ZIRS: the fields of the textZIRS packed into a",
        "64-bit integer; cell_to_index(form = \"uint64ZIRS\").")))
  } else {
    list(textZIRS = list(description = isea_index_description(g),
                         type = "hierarchicalConcatenation"),
         uint64ZIRS = list(description = paste(
           "The cell ID: 1 to N in the order DGGRID numbers the cells of one",
           "resolution (SEQNUM), unique together with the resolution.")))
  }

  out <- list(
    title = dggrs_title(g),
    description = paste(c(sprintf(paste(
      "hexify %s grid at resolution %d: %s cells on the %s, refined by",
      "aperture %s,"),
      dggrs_title(g), g@resolution, toupper(projection), polyhedron,
      paste(steps, collapse = ", ")), earth_model,
      if (!is.null(ogc)) paste0("Its hierarchy and zone identifiers are those of OGC's ", ogc, ".")),
      collapse = " "),
    dggh = list(
      definition = list(
        spatialDimensions = 2L,
        temporalDimensions = 0L,
        crs = dggrs_crs(g),
        basePolyhedron = polyhedron,
        refinementRatio = as.integer(steps),
        refinementStrategy = unname(strategy),
        constraints = list(cellEqualSized = is_equal_area_projection(projection)),
        zoneTypes = c("hexagon", vertex_cell)
      ),
      parameters = parameters
    ),
    zirs = zirs,
    subZoneOrder = list(
      description = "Sub-zones at any depth below a zone are listed in ascending order of cell ID.",
      type = "cellIdAscending"
    )
  )
  if (!is.null(ogc)) {
    out <- append(out, list(links = list(list(
      rel = "related", href = paste0("https://www.opengis.net/def/dggrs/OGC/1.0/", ogc)))),
      after = 2L)
  }
  out
}

#' The name of the OGC DGGRS definition a grid is, or NULL
#'
#' OGC's ISEA3H, ISEA7H, IVEA3H and IVEA7H: aperture 3 or 7 on the
#' icosahedron, Snyder's or the vertex-oriented projection, WGS84 read through
#' authalic latitude, and the OGC orientation.
#' @noRd
ogc_dggrs_name <- function(g) {
  if (is_h3_grid(g) || grid_polyhedron(g) != "icosahedron" ||
      !g@aperture %in% c("3", "7") || !grid_projection(g) %in% c("isea", "ivea") ||
      !identical(unname(grid_ellipsoid(g)), unname(ELLIPSOIDS$wgs84)) ||
      !isTRUE(all.equal(unname(grid_orientation(g)), unname(OGC_ORIENTATION),
                        tolerance = 0))) {
    return(NULL)
  }
  paste0(toupper(grid_projection(g)), g@aperture, "H")
}

#' Name of an ISEA-family grid in DGGRID's style
#' @noRd
dggrs_title <- function(g) {
  ap <- as.character(g@aperture)
  body <- if (is_per_level_aperture(ap)) sprintf("[%s]", ap) else gsub("/", "", ap, fixed = TRUE)
  title <- paste0(toupper(grid_projection(g)), body, "H")
  if (grid_polyhedron(g) != "icosahedron") title <- paste0(title, "-", grid_polyhedron(g))
  title
}

#' A grid's CRS as an identifier string
#' @noRd
dggrs_crs <- function(g) {
  crs <- sf::st_crs(grid_crs(g))
  if (!is.na(crs$epsg)) sprintf("EPSG:%d", crs$epsg) else crs$proj4string
}

#' The index string cell_to_index() writes for a grid, in words
#' @noRd
isea_index_description <- function(g) {
  if (is_mixed_aperture(g@aperture)) {
    return(paste(
      "Two digits naming the resolution-0 cell the zone descends from, then two",
      "digits per resolution giving which child of its parent the zone is,",
      "children counted from 00 in ascending order of cell ID. A zone's parent",
      "is the coarser zone holding its centre; dropping the last two digits",
      "gives the parent's identifier."))
  }
  n_quad <- sprintf("00-%02d", polyhedron_diamonds(grid_polyhedron(g)) + 1L)
  switch(index_type_for_aperture(g@aperture),
    z7 = if (grid_polyhedron(g) == "icosahedron") paste(
      "Z7 (IGEO7): two digits naming the resolution-0 cell (00-11) the zone",
      "descends from, then one digit 0-6 per resolution naming the child of",
      "the parent: 0 the centred child, 1-6 the six around it; the string",
      "DGGRID writes. Dropping the last digit gives the parent's identifier.")
    else paste(
      "Z7: a leading field quad + 6 * s, s the digit of the neighbour of the",
      "quad's vertex the zone descends from, then one digit 0-6 per resolution",
      "naming the child of the parent: 0 the centred child, 1-6 the six around",
      "it. Dropping the last digit gives the parent's identifier."),
    z3 = paste0(
      "Z3: two digits naming the quad (", n_quad, "), then one digit pair per ",
      "two resolutions encoding the aperture-3 path; dropping the last digits ",
      "gives the ancestor's identifier."),
    zorder = paste0(
      "Z-order: two digits naming the quad (", n_quad, "), then the base-2 ",
      "interleaved digits of the zone's (i, j) lattice coordinates, one per ",
      "resolution; dropping the last digit gives the parent's identifier."))
}

#' DGGRS definition of an H3 grid
#' @noRd
h3_dggrs_definition <- function(g) {
  list(
    title = "H3",
    description = sprintf(
      "H3 grid at resolution %d: aperture-7 hexagons and twelve pentagons on a gnomonic projection of each icosahedron face, on a sphere of radius %s km.",
      g@resolution, format(grid_radius_km(g), digits = 10)),
    dggh = list(
      definition = list(
        spatialDimensions = 2L,
        temporalDimensions = 0L,
        crs = dggrs_crs(g),
        basePolyhedron = "icosahedron",
        refinementRatio = 7L,
        refinementStrategy = c("centredChildCell", "nodeSharingChildCell"),
        constraints = list(cellEqualSized = FALSE),
        zoneTypes = c("hexagon", "pentagon")
      ),
      parameters = list(sphere = list(radius_km = grid_radius_km(g)))
    ),
    zirs = list(
      textZIRS = list(description = "The H3 index as 15 lowercase hexadecimal digits.",
                      type = "hierarchicalConcatenation"),
      uint64ZIRS = list(description = "The H3 index as a 64-bit integer, unique across all resolutions.")
    ),
    subZoneOrder = list(
      description = "Sub-zones at any depth below a zone are listed in ascending order of H3 index.",
      type = "cellIdAscending"
    )
  )
}
