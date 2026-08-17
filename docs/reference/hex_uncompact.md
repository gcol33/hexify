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

**ISEA backend (aperture 7, Z7 index):** Appends digits 0-6 to expand
each cell to its 7 children, repeating until the target resolution is
reached.

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
```
