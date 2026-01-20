# Get parent cell

Returns the parent cell at a coarser resolution.

## Usage

``` r
get_parent(cell_id, grid, levels = 1L)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs

- grid:

  A HexGridInfo or HexData object

- levels:

  Number of levels up (default 1)

## Value

Numeric vector of parent cell IDs

## Examples

``` r
grid <- hex_grid(resolution = 10)
child_cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
parent_cells <- get_parent(child_cells, grid)
```
