# Sibling cells

Returns the other children of each cell's parent.

## Usage

``` r
get_siblings(cell_id, grid, include_self = FALSE)
```

## Arguments

- cell_id:

  Cell IDs: integer64 for ISEA grids (whole numbers below 2^53 and digit
  strings are accepted), character for H3 grids

- grid:

  A HexGridInfo or HexData object

- include_self:

  If `TRUE`, keep the cell itself among its siblings.

## Value

List with one vector of cell IDs per input cell.

## Details

The parent is the one
[`get_parent`](https://gillescolling.com/hexify/reference/get_parent.md)
returns, the coarser cell holding the cell's centre, and the siblings
are its children as
[`get_children`](https://gillescolling.com/hexify/reference/get_children.md)
returns them: the `siblingOf` relation of OGC Topic 21 with `inheritID`
true.

## See also

[`get_parent`](https://gillescolling.com/hexify/reference/get_parent.md),
[`get_children`](https://gillescolling.com/hexify/reference/get_children.md)

## Examples

``` r
g <- hex_grid(resolution = 4, aperture = 7)
cell <- lonlat_to_cell(16.37, 48.21, g)
get_siblings(cell, g)
```
