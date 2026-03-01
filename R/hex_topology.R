# R/hex_topology.R
# Pentagon detection and topological queries

#' Detect Pentagon Cells
#'
#' Identifies which cells are pentagons. Any hexagonal tiling of the sphere
#' must contain exactly 12 pentagons (at the icosahedron vertices).
#' Pentagon cells have 5 neighbors instead of 6.
#'
#' @param cell_id Cell IDs to check. Numeric for ISEA, character for H3.
#' @param grid A HexGridInfo or HexData object specifying the grid.
#'
#' @return A logical vector. `TRUE` for pentagon cells, `FALSE` for hexagons.
#'
#' @details
#' **H3 backend:** Uses the vendored H3 `isPentagon` function.
#'
#' **ISEA backend:** The 12 pentagons are located at icosahedron vertices.
#' At any resolution, the pentagon cells are those whose IJ coordinates
#' are (0, 0) within the 12 vertex quads (quads 0 and 11 for the poles,
#' and specific cells in quads 1-10).
#'
#' @seealso [get_neighbors()] for neighbor finding (pentagons have 5 neighbors)
#'
#' @export
#' @examples
#' \donttest{
#' # H3 pentagon detection
#' g <- hex_grid(resolution = 1, type = "h3")
#' cells <- grid_global(g)
#' pent <- is_pentagon(cells$cell_id, g)
#' sum(pent)  # Should be 12
#' }
is_pentagon <- function(cell_id, grid) {
  g <- extract_grid(grid)

  if (is_h3_grid(g)) {
    return(as.logical(cpp_h3_isPentagon(as.character(cell_id))))
  }

  # ISEA backend: pentagon cells are at icosahedron vertices
  # The 12 pentagons correspond to cells at (i=0, j=0) in:
  #   - Quad 0 (north pole)
  #   - Quads 1-10: the cell at origin of each quad
  #   - Quad 11 (south pole)
  # At resolution 0, cell_id 1 through 12 are the 12 quads (all pentagons)
  # At higher resolutions, pentagon = cell_id corresponding to (quad, 0, 0)
  ap <- as.integer(g@aperture)
  res <- g@resolution

  if (res == 0L) {
    # At resolution 0, all 12 cells are pentagons
    return(rep(TRUE, length(cell_id)))
  }

  # Calculate which cell_ids correspond to (quad, 0, 0)
  nCells <- 0
  offsetPerQuad <- 0
  # Re-derive pentagon cell IDs
  pentagon_ids <- vapply(0:11, function(q) {
    ll <- cell_to_lonlat(q + 1L, hex_grid(resolution = 0L, aperture = g@aperture))
    cid <- lonlat_to_cell(ll$lon_deg, ll$lat_deg, g)
    cid
  }, numeric(1))

  cell_id %in% pentagon_ids
}
