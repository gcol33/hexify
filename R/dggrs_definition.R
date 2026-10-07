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
#' The ISEA3H and ISEA7H definitions registered with OGC
#' (\url{https://www.opengis.net/def/dggrs/OGC/1.0/ISEA7H}) are other grids:
#' they convert WGS84 geodetic latitude to authalic latitude before
#' projecting, and place the icosahedron's first vertex at 11.20 degrees E
#' and authalic latitude \eqn{\arctan(\phi)}, where hexify, like DGGRID, reads
#' latitude on the sphere and places it at 11.25 degrees E. A hexify grid is
#' therefore described without an OGC URI.
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

  list(
    title = dggrs_title(g),
    description = sprintf(paste(
      "hexify %s grid at resolution %d: %s cells on the %s, refined by",
      "aperture %s, on a sphere of radius %s km taking latitude as spherical",
      "latitude."),
      dggrs_title(g), g@resolution, toupper(projection), polyhedron,
      paste(steps, collapse = ", "),
      format(grid_radius_km(g), digits = 10)),
    dggh = list(
      definition = list(
        spatialDimensions = 2L,
        temporalDimensions = 0L,
        crs = dggrs_crs(g),
        basePolyhedron = polyhedron,
        refinementRatio = as.integer(steps),
        refinementStrategy = unname(strategy),
        constraints = list(cellEqualSized = projection != "fuller"),
        zoneTypes = c("hexagon", vertex_cell)
      ),
      parameters = list(
        sphere = list(radius_km = grid_radius_km(g)),
        orientation = list(
          latitude = unname(o[["vert0_lat"]]),
          longitude = unname(o[["vert0_lon"]]),
          azimuth = unname(o[["azimuth"]]),
          description = "Spherical latitude and longitude of the solid's first vertex, and azimuth of its second vertex seen from the first."
        )
      )
    ),
    zirs = list(
      textZIRS = list(description = isea_index_description(g),
                      type = "hierarchicalConcatenation"),
      uint64ZIRS = list(description = paste(
        "The cell ID: 1 to N in the order DGGRID numbers the cells of one",
        "resolution (SEQNUM), unique together with the resolution."))
    ),
    subZoneOrder = list(
      description = "Sub-zones at any depth below a zone are listed in ascending order of cell ID.",
      type = "cellIdAscending"
    )
  )
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
