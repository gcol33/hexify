# R/hex_compact.R
# Multi-resolution compaction and uncompaction

#' Compact Hex Cells
#'
#' Merges child cells into their parent when all children are present.
#' This is a lossless compression — no spatial information is lost.
#' The compact representation uses fewer cells to cover the same area.
#'
#' @param cell_ids Cell IDs to compact. For H3 grids, a character vector.
#'   For ISEA grids, a character vector of hierarchical index strings.
#' @param grid A HexGridInfo object specifying the grid.
#'
#' @return A character vector of compacted cell IDs. Cells that could be
#'   merged into parents appear as parent IDs at coarser resolution.
#'
#' @details
#' **H3 backend:** Uses the vendored H3 `compactCells` function.
#'
#' **ISEA backend (aperture 7, Z7 index):** Groups cells by parent index
#' (dropping the last digit). If all 7 children are present, replaces them
#' with the parent. Iterates until no further compaction is possible.
#'
#' @seealso [hex_uncompact()] for the inverse operation,
#'   [get_parent()], [get_children()] for hierarchical operations
#'
#' @export
#' @examples
#' \donttest{
#' # H3 compaction
#' g <- hex_grid(resolution = 3, type = "h3")
#' parent <- "832830fffffffff"
#' children <- get_children(parent, g)[[1]]
#' compact <- hex_compact(children, g)
#' compact  # Should return the parent
#' }
hex_compact <- function(cell_ids, grid) {
  g <- extract_grid(grid)

  if (is_h3_grid(g)) {
    return(cpp_h3_compactCells(as.character(cell_ids)))
  }

  # ISEA Z7 compaction (pure R)
  ap <- as.integer(g@aperture)
  if (ap != 7) {
    stop("hex_compact() for ISEA currently only supports aperture 7 (Z7 index)")
  }

  ids <- as.character(cell_ids)
  changed <- TRUE

  while (changed) {
    changed <- FALSE
    # Group by parent (drop last digit)
    parents <- substr(ids, 1, nchar(ids) - 1L)
    # Only cells with length > 2 can be compacted (resolution > 0)
    compactable <- nchar(ids) > 2L

    if (!any(compactable)) break

    parent_table <- table(parents[compactable])
    # Parents with all 7 children present
    full_parents <- names(parent_table[parent_table == 7L])

    if (length(full_parents) > 0) {
      changed <- TRUE
      # Remove children of full parents, add parent
      is_child_of_full <- compactable & parents %in% full_parents
      remaining <- ids[!is_child_of_full]
      ids <- c(remaining, full_parents)
    }
  }

  ids
}

#' Uncompact Hex Cells
#'
#' Expands compacted cells to a target resolution. All cells in the output
#' share the same resolution.
#'
#' @param cell_ids Character vector of (possibly mixed-resolution) cell IDs.
#' @param grid A HexGridInfo object specifying the grid.
#' @param target_resolution Integer. The resolution to expand all cells to.
#'
#' @return A character vector of cell IDs, all at `target_resolution`.
#'
#' @details
#' **H3 backend:** Uses the vendored H3 `uncompactCells` function.
#'
#' **ISEA backend (aperture 7, Z7 index):** Appends digits 0-6 to expand
#' each cell to its 7 children, repeating until the target resolution is
#' reached.
#'
#' @seealso [hex_compact()] for the inverse operation
#'
#' @export
#' @examples
#' \donttest{
#' g <- hex_grid(resolution = 3, type = "h3")
#' parent <- "832830fffffffff"
#' hex_uncompact(parent, g, target_resolution = 4L)
#' }
hex_uncompact <- function(cell_ids, grid, target_resolution) {
  g <- extract_grid(grid)
  target_resolution <- as.integer(target_resolution)

  if (is_h3_grid(g)) {
    return(cpp_h3_uncompactCells(as.character(cell_ids), target_resolution))
  }

  # ISEA Z7 uncompaction (pure R)
  ap <- as.integer(g@aperture)
  if (ap != 7) {
    stop("hex_uncompact() for ISEA currently only supports aperture 7 (Z7 index)")
  }

  ids <- as.character(cell_ids)

  repeat {
    # Check which cells need expansion
    current_res <- nchar(ids) - 2L  # Z7: first 2 chars are quad
    needs_expansion <- current_res < target_resolution

    if (!any(needs_expansion)) break

    # Expand cells that need it
    expanded <- unlist(lapply(ids[needs_expansion], function(id) {
      paste0(id, 0:6)
    }))

    ids <- c(ids[!needs_expansion], expanded)
  }

  ids
}
