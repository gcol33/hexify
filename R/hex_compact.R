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

    # Parents with all 7 distinct children (0-6) present. Group by distinct
    # last digit rather than raw row count, so duplicate cell IDs can't be
    # mistaken for a full sibling set.
    compactable_ids <- ids[compactable]
    compactable_parents <- parents[compactable]
    last_digit <- substr(compactable_ids, nchar(compactable_ids), nchar(compactable_ids))
    digit_sets <- split(last_digit, compactable_parents)
    full_parents <- names(digit_sets)[vapply(digit_sets, function(d) length(unique(d)) == 7L, logical(1))]

    if (length(full_parents) > 0) {
      changed <- TRUE
      # Remove children of full parents, add parent. If the input already
      # contained the parent id itself alongside all 7 children, that parent
      # entry's own `parents` value (its grandparent) never matches
      # `full_parents`, so it would survive untouched while the
      # newly-synthesized parent is also appended -- drop it explicitly to
      # avoid emitting the same parent id twice.
      is_child_of_full <- compactable & parents %in% full_parents
      is_full_parent_itself <- ids %in% full_parents
      remaining <- ids[!is_child_of_full & !is_full_parent_itself]
      ids <- c(remaining, full_parents)
    }
  }

  # Mixed-resolution input can also leave a coarser ancestor and one of its
  # own (not-fully-compactable) descendants both present in the result --
  # the ancestor's area already covers the descendant's, so the descendant
  # is redundant even though it's not a literal duplicate ID.
  if (length(ids) > 1L) {
    is_redundant <- vapply(seq_along(ids), function(k) {
      any(nchar(ids) < nchar(ids[k]) & substr(ids[k], 1, nchar(ids)) == ids)
    }, logical(1))
    ids <- ids[!is_redundant]
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

  initial_res <- nchar(ids) - 2L  # Z7: first 2 chars are quad
  if (any(initial_res > target_resolution)) {
    stop(sprintf(
      "target_resolution (%d) is coarser than some input cells (max resolution %d); hex_uncompact() cannot expand to a coarser resolution",
      target_resolution, max(initial_res)
    ))
  }

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
