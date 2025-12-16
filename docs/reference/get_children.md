# Get children cells

Returns the child cells at a finer resolution.

## Usage

``` r
get_children(cell_id, grid, levels = 1L)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs

- grid:

  A HexGrid or HexData object

- levels:

  Number of levels down (default 1)

## Value

List of numeric vectors containing child cell IDs
