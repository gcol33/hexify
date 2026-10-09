# R/hex_compact.R
# Multi-resolution compaction and uncompaction

#' The resolution an ISEA index string carries
#'
#' A Z-order, Z3 or Z7 index spells its quad in two characters and then one
#' digit per resolution; a Hex9 label spells one digit per resolution from
#' resolution 0, then a dot and its key tail.
#' @param indices Character vector of index strings
#' @param g The grid the cells belong to
#' @noRd
index_resolution <- function(indices, g) {
  if (is_hex9_grid(g)) return(nchar(sub(".", "", indices, fixed = TRUE)) - 2L)
  nchar(indices) - 2L
}

#' The cells an ISEA index string names
#' @noRd
isea_index_cells <- function(indices, g) {
  if (is_hex9_grid(g)) return(cpp_hex9_parse_label(as.character(indices))$cell_id)
  isea_index_to_cells(indices, aperture_to_int(g@aperture),
                      index_type_for_aperture(g@aperture), icosa_arg(g))
}

#' Stop unless every index names a cell of the grid
#'
#' An index string is arithmetic, so one can be written that names nothing --
#' by appending a digit to a cell at a vertex of the solid, for instance. Saying
#' so here names the string, rather than letting a cell-ID range check further
#' down report a number the caller never wrote.
#'
#' @param indices Character vector of index strings
#' @param g The grid the cells belong to
#' @param what The calling function, for the message
#' @return The cells the indices name
#' @noRd
check_isea_indices <- function(indices, g, what) {
  resolutions <- index_resolution(indices, g)
  if (any(resolutions < 0L)) {
    stop(what, "(): ", "index strings carry a two-character quad and one digit ",
         "per resolution; these are shorter than that: ",
         paste(utils::head(indices[resolutions < 0L], 5), collapse = ", "),
         call. = FALSE)
  }

  cells <- isea_index_cells(indices, g)
  limit <- vapply(resolutions,
                  function(r) aperture_n_cells(g@aperture, r, grid_polyhedron(g)),
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
  children <- get_children(cells, grid_at_resolution(g, resolution))
  cell_to_index(unique(cell_id_unlist(children)),
                grid_at_resolution(g, resolution + 1L))
}

#' Compact Hex Cells
#'
#' Merges child cells into their parent when all children are present, so that
#' the same set of cells is named by fewer of them.
#'
#' @param cell_ids Cell IDs to compact. For H3 grids, a character vector.
#'   For ISEA grids, a character vector of hierarchical index strings, as
#'   \code{\link{cell_to_index}} returns for apertures 3, 4, 7 and 9.
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
  levels <- rev(seq_len(max(index_resolution(ids, g))))

  for (level in levels) {
    at_level <- index_resolution(ids, g) == level
    if (!any(at_level)) next

    parent_grid <- grid_at_resolution(g, level - 1L)
    cells <- isea_index_cells(ids[at_level], g)
    groups <- split(cells, get_parent(cells, grid_at_resolution(g, level)))
    candidates <- as_cell_id(names(groups))

    # Every cell in a group has this parent, so a group holding as many
    # distinct cells as the parent has children holds all of them. That is the
    # aperture, fewer at the vertex cells of the solid, which are read rather
    # than assumed.
    expected <- rep(ap, length(candidates))
    vertex <- isea_cell_sides(candidates, parent_grid) != 6L
    if (any(vertex)) {
      expected[vertex] <- lengths(get_children(candidates[vertex], parent_grid))
    }

    sizes <- vapply(groups, function(cells) length(unique(cells)), integer(1))
    full <- sizes == expected
    if (!any(full)) next

    merged <- cell_to_index(candidates[full], parent_grid)
    covered <- cell_id_unlist(groups[full])
    ids <- unique(c(setdiff(ids[!at_level], merged),
                    ids[at_level][!(cells %in% covered)],
                    merged))
  }

  # Mixed-resolution input can also leave a coarser ancestor and one of its
  # own (not-fully-compactable) descendants both present in the result --
  # the ancestor's area already covers the descendant's, so the descendant
  # is redundant even though it's not a literal duplicate ID.
  if (length(ids) > 1L) ids <- ids[!has_listed_ancestor(ids, g)]

  ids
}

#' Whether each index string has an ancestor among the others
#'
#' The ancestors are the iterated parents of get_parent(), read up to the
#' coarsest resolution present.
#' @param ids Character vector of index strings
#' @param g The grid the cells belong to
#' @return Logical vector
#' @noRd
has_listed_ancestor <- function(ids, g) {
  res <- index_resolution(ids, g)
  cells <- isea_index_cells(ids, g)
  listed <- paste(res, as.character(cells))
  out <- logical(length(ids))
  if (max(res) == min(res)) return(out)
  up <- cells
  for (r in seq.int(max(res), min(res) + 1L)) {
    moving <- res >= r
    up[moving] <- get_parent(up[moving], grid_at_resolution(g, r))
    hit <- moving & paste(r - 1L, as.character(up)) %in% listed
    out <- out | hit
  }
  out
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

  initial_res <- index_resolution(ids, g)
  if (any(initial_res > target_resolution)) {
    stop(sprintf(
      "target_resolution (%d) is coarser than some input cells (max resolution %d); hex_uncompact() cannot expand to a coarser resolution",
      target_resolution, max(initial_res)
    ))
  }

  repeat {
    current_res <- index_resolution(ids, g)
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
