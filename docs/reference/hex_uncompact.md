# Uncompact Hex Cells

Expands compacted cells to a target resolution. All cells in the output
share the same resolution.

## Usage

``` r
hex_uncompact(cell_ids, grid, target_resolution)
```

## Arguments

- cell_ids:

  Character vector of (possibly mixed-resolution) cell IDs.

- grid:

  A HexGridInfo object specifying the grid.

- target_resolution:

  Integer. The resolution to expand all cells to.

## Value

A character vector of cell IDs, all at `target_resolution`.

## Details

**H3 backend:** Uses the vendored H3 `uncompactCells` function.

**ISEA backend:** Cells are expanded through
[`get_children()`](https://gillescolling.com/hexify/reference/get_children.md),
a level at a time, until the target resolution is reached. Appending
each digit of the aperture names the children of a cell whose whole
ancestry lies inside one quad, but not those of a cell at an icosahedron
vertex, whose children carry index spellings in the several quads
meeting there.

## See also

[`hex_compact()`](https://gillescolling.com/hexify/reference/hex_compact.md)
for the inverse operation

## Examples

``` r
# \donttest{
g <- hex_grid(resolution = 3, type = "h3")
parent <- "832830fffffffff"
hex_uncompact(parent, g, target_resolution = 4L)
# }

# ISEA, on the default aperture
g <- hex_grid(resolution = 2, aperture = 3)
hex_uncompact(cell_to_index(40L, g), g, target_resolution = 3L)
```
