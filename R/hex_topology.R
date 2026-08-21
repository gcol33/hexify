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
#' **ISEA backend:** The 12 pentagons are located at icosahedron vertices,
#' which are always the (i, j) = (0, 0) cell of their quad. Pentagon status
#' is checked by decoding each input cell ID's own (quad, i, j) coordinates
#' via [cell_to_lonlat()]'s underlying grid math and testing whether i and j
#' are both zero, rather than by re-deriving each vertex's cell ID (the
#' forward direction has aperture-specific quirks -- e.g. aperture 7's
#' substrate/surrogate coordinate distinction -- that make a single fixed
#' formula for "the (0,0) cell ID of quad Q" unreliable across apertures).
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

  if (g@resolution == 0L) {
    # At resolution 0, all 12 cells (the icosahedron vertices) are pentagons
    return(rep(TRUE, length(cell_id)))
  }

  # Pentagon cells are exactly the (i, j) = (0, 0) cell of each quad.
  qij <- if (is_mixed_aperture(g@aperture)) {
    mixed_cell_qij(cell_id, g@resolution, g@aperture)
  } else {
    cpp_cell_to_quad_ij(cell_id, g@resolution, aperture_to_int(g@aperture))
  }
  qij$i == 0 & qij$j == 0
}
