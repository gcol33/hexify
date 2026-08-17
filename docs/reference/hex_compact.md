# Compact Hex Cells

Merges child cells into their parent when all children are present. This
is a lossless compression — no spatial information is lost. The compact
representation uses fewer cells to cover the same area.

## Usage

``` r
hex_compact(cell_ids, grid)
```

## Arguments

- cell_ids:

  Cell IDs to compact. For H3 grids, a character vector. For ISEA grids,
  a character vector of hierarchical index strings.

- grid:

  A HexGridInfo object specifying the grid.

## Value

A character vector of compacted cell IDs. Cells that could be merged
into parents appear as parent IDs at coarser resolution.

## Details

**H3 backend:** Uses the vendored H3 `compactCells` function.

**ISEA backend (aperture 7, Z7 index):** Groups cells by parent index
(dropping the last digit). If all 7 children are present, replaces them
with the parent. Iterates until no further compaction is possible.

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
```
