# Compact Hex Cells

Merges child cells into their parent when all children are present, so
that the same set of cells is named by fewer of them.

## Usage

``` r
hex_compact(cell_ids, grid)
```

## Arguments

- cell_ids:

  Cell IDs to compact. For H3 grids, a character vector. For ISEA grids,
  a character vector of hierarchical index strings, as
  [`cell_to_index`](https://gillescolling.com/hexify/reference/cell_to_index.md)
  returns for apertures 3, 4 and 7.

- grid:

  A HexGridInfo object specifying the grid.

## Value

A character vector of compacted cell IDs. Cells that could be merged
into parents appear as parent IDs at coarser resolution.

## Details

**H3 backend:** Uses the vendored H3 `compactCells` function.

**ISEA backend:** Cells are grouped by
[`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md),
and a group holding as many distinct cells as its parent has children
replaces them, from the finest resolution up until no further compaction
is possible. That count is the aperture away from the twelve icosahedron
vertices and fewer at one, and it is read rather than assumed: a vertex
cell sits at the corner of several quads and carries an index spelling
in each, so appending a digit to any one spelling names neither all of
its children nor only its children.

What is preserved is the set of cells: uncompacting the result at the
original resolution returns the input. The covered area is not
identical, because a hexagonal hierarchy is not congruent at any
aperture – a parent does not tile exactly into its children – so an area
computed on the compacted set differs from one computed on the original.

## See also

[`hex_uncompact()`](https://gillescolling.com/hexify/reference/hex_uncompact.md)
for the inverse operation,
[`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md),
[`get_children()`](https://gillescolling.com/hexify/reference/get_children.md)
for hierarchical operations

## Examples

``` r
# \donttest{
# H3 compaction
g <- hex_grid(resolution = 3, type = "h3")
parent <- "832830fffffffff"
children <- get_children(parent, g)[[1]]
compact <- hex_compact(children, g)
compact  # Should return the parent
# }

# ISEA, on the default aperture
g <- hex_grid(resolution = 2, aperture = 3)
children <- cell_to_index(get_children(40L, g)[[1]],
                          hex_grid(resolution = 3, aperture = 3))
hex_compact(children, g)
```
