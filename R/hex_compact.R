# R/hex_compact.R
# Multi-resolution compaction and uncompaction

#' A grid's own resolutions, for reading indexed cells at each of them
#' @noRd
isea_grid_at <- function(g, resolution) {
  hex_grid(resolution = resolution, aperture = g@aperture,
           radius_km = grid_radius_km(g))
}

#' The cells an ISEA index string names
#' @noRd
isea_index_cells <- function(indices, g) {
  isea_index_to_cells(indices, aperture_to_int(g@aperture),
                      index_type_for_aperture(g@aperture))
}

#' Stop unless every index names a cell of the grid
#'
#' An index string is arithmetic, so one can be written that names nothing --
#' by appending a digit to a cell at an icosahedron vertex, for instance. Saying
#' so here names the string, rather than letting a cell-ID range check further
#' down report a number the caller never wrote.
#'
#' @param indices Character vector of index strings
#' @param g The grid the cells belong to
#' @param what The calling function, for the message
#' @return The cells the indices name
#' @noRd
check_isea_indices <- function(indices, g, what) {
  resolutions <- nchar(indices) - 2L
  if (any(resolutions < 0L)) {
    stop(what, "(): ", "index strings carry a two-character quad and one digit ",
         "per resolution; these are shorter than that: ",
         paste(utils::head(indices[resolutions < 0L], 5), collapse = ", "),
         call. = FALSE)
  }

  cells <- isea_index_cells(indices, g)
  limit <- vapply(resolutions, function(r) aperture_n_cells(g@aperture, r),
                  numeric(1))
  named_nothing <- is.na(cells) | cells < 1 | cells > limit

  if (any(named_nothing)) {
    stop(what, "(): these index strings name no cell of the grid: ",
         paste(utils::head(indices[named_nothing], 5), collapse = ", "),
         if (sum(named_nothing) > 5) {
           sprintf(" (and %d more)", sum(named_nothing) - 5)
         } else {
           ""
         },
         call. = FALSE)
  }

  cells
}

#' The index strings of the children of indexed cells, one level down
#'
#' Appending a digit to the index would name the children, except at the twelve
#' icosahedron vertices: a vertex sits at the corner of several quads and
#' carries an index spelling in each, so a digit appended to one spelling can
#' name a cell held by another quad, or nothing at all. The cells are read
#' instead, and written back as the index each of them carries.
#'
#' @param indices Character vector of index strings, all at `resolution`
#' @param resolution Resolution the indices are at
#' @param g The grid the cells belong to
#' @return Character vector of child index strings
#' @noRd
isea_child_indices <- function(indices, resolution, g) {
  cells <- isea_index_cells(indices, g)
  children <- get_children(cells, isea_grid_at(g, resolution))
  cell_to_index(unique(unlist(children, use.names = FALSE)),
                isea_grid_at(g, resolution + 1L))
}

#' Compact Hex Cells
#'
#' Merges child cells into their parent when all children are present, so that
#' the same set of cells is named by fewer of them.
#'
#' @param cell_ids Cell IDs to compact. For H3 grids, a character vector.
#'   For ISEA grids, a character vector of hierarchical index strings, as
#'   \code{\link{cell_to_index}} returns for apertures 3, 4 and 7.
#' @param grid A HexGridInfo object specifying the grid.
#'
#' @return A character vector of compacted cell IDs. Cells that could be
#'   merged into parents appear as parent IDs at coarser resolution.
#'
#' @details
#' **H3 backend:** Uses the vendored H3 `compactCells` function.
#'
#' **ISEA backend:** Cells are grouped by [get_parent()], and a group holding as
#' many distinct cells as its parent has children replaces them, from the finest
#' resolution up until no further compaction is possible. That count is the
#' aperture away from the twelve icosahedron vertices and fewer at one, and it
#' is read rather than assumed: a vertex cell sits at the corner of several
#' quads and carries an index spelling in each, so appending a digit to any one
#' spelling names neither all of its children nor only its children.
#'
#' What is preserved is the set of cells: uncompacting the result at the
#' original resolution returns the input. The covered area is not identical,
#' because a hexagonal hierarchy is not congruent at any aperture -- a parent
#' does not tile exactly into its children -- so an area computed on the
#' compacted set differs from one computed on the original.
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
#'
#' # ISEA, on the default aperture
#' g <- hex_grid(resolution = 2, aperture = 3)
#' children <- cell_to_index(get_children(40L, g)[[1]],
#'                           hex_grid(resolution = 3, aperture = 3))
#' hex_compact(children, g)
hex_compact <- function(cell_ids, grid) {
  g <- extract_grid(grid)

  if (is_h3_grid(g)) {
    return(cpp_h3_compactCells(as.character(cell_ids)))
  }

  # ISEA compaction (pure R), on the per-level digit every ISEA index carries
  ap <- aperture_to_int(g@aperture)
  if (is_mixed_aperture(g@aperture)) {
    stop("hex_compact() has no index arithmetic for a mixed aperture sequence")
  }

  ids <- unique(as.character(cell_ids))
  check_isea_indices(ids, g, "hex_compact")

  # Levels are taken from the finest up: each level's full sibling sets become
  # parents at the next one, which the following pass reads.
  levels <- rev(seq_len(max(nchar(ids) - 2L)))

  for (level in levels) {
    at_level <- nchar(ids) - 2L == level
    if (!any(at_level)) next

    parent_grid <- isea_grid_at(g, level - 1L)
    cells <- isea_index_cells(ids[at_level], g)
    groups <- split(cells, get_parent(cells, isea_grid_at(g, level)))
    candidates <- as.numeric(names(groups))

    # Every cell in a group has this parent, so a group holding as many
    # distinct cells as the parent has children holds all of them. That is the
    # aperture, fewer at the twelve icosahedron vertices, which are read rather
    # than assumed.
    expected <- rep(ap, length(candidates))
    vertex <- is_pentagon(candidates, parent_grid)
    if (any(vertex)) {
      expected[vertex] <- lengths(get_children(candidates[vertex], parent_grid))
    }

    sizes <- vapply(groups, function(cells) length(unique(cells)), integer(1))
    full <- sizes == expected
    if (!any(full)) next

    merged <- cell_to_index(candidates[full], parent_grid)
    covered <- unlist(groups[full], use.names = FALSE)
    ids <- unique(c(setdiff(ids[!at_level], merged),
                    ids[at_level][!(cells %in% covered)],
                    merged))
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
#' **ISEA backend:** Cells are expanded through [get_children()], a level at a
#' time, until the target resolution is reached. Appending each digit of the
#' aperture names the children of a cell whose whole ancestry lies inside one
#' quad, but not those of a cell at an icosahedron vertex, whose children carry
#' index spellings in the several quads meeting there.
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
#'
#' # ISEA, on the default aperture
#' g <- hex_grid(resolution = 2, aperture = 3)
#' hex_uncompact(cell_to_index(40L, g), g, target_resolution = 3L)
hex_uncompact <- function(cell_ids, grid, target_resolution) {
  g <- extract_grid(grid)
  target_resolution <- as.integer(target_resolution)

  if (is_h3_grid(g)) {
    return(cpp_h3_uncompactCells(as.character(cell_ids), target_resolution))
  }

  # ISEA uncompaction (pure R)
  if (is_mixed_aperture(g@aperture)) {
    stop("hex_uncompact() has no index arithmetic for a mixed aperture sequence")
  }

  ids <- as.character(cell_ids)
  check_isea_indices(ids, g, "hex_uncompact")

  initial_res <- nchar(ids) - 2L  # the first 2 chars are the quad
  if (any(initial_res > target_resolution)) {
    stop(sprintf(
      "target_resolution (%d) is coarser than some input cells (max resolution %d); hex_uncompact() cannot expand to a coarser resolution",
      target_resolution, max(initial_res)
    ))
  }

  repeat {
    current_res <- nchar(ids) - 2L
    needs_expansion <- current_res < target_resolution

    if (!any(needs_expansion)) break

    # A level at a time, so each cell's children are read at its own resolution
    coarsest <- min(current_res[needs_expansion])
    expanding <- needs_expansion & current_res == coarsest

    ids <- c(ids[!expanding],
             isea_child_indices(ids[expanding], coarsest, g))
  }

  unique(ids)
}
