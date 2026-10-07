# Convert cell ID to hierarchical index string

Advanced function for working with hierarchical index strings. Most
users don't need this - use cell IDs directly.

## Usage

``` r
cell_to_index(cell_id, grid)
```

## Arguments

- cell_id:

  Cell IDs: integer64 for ISEA grids (whole numbers below 2^53 and digit
  strings are accepted), character for H3 grids

- grid:

  A HexGridInfo or HexData object

## Value

Character vector of hierarchical index strings

## Details

A cell ID numbers the cells of one resolution from 1, so it names a cell
only together with its grid. The index string carries its resolution in
its length and names one cell among all resolutions of the grid's
family, and so serves as the zonal identifier of OGC Topic 21; an H3
index carries its resolution in its bits.
