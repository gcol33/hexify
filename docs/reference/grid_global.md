# Generate a global hexagon grid

Creates hexagon polygons covering the entire Earth.

## Usage

``` r
grid_global(grid)
```

## Arguments

- grid:

  A HexGridInfo object specifying the grid parameters

## Value

sf object with hexagon polygons

## Details

This function generates a complete global grid by sampling points
densely across the globe. For large grids (many small cells), consider
using
[`grid_rect()`](https://gcol33.github.io/hexify/reference/grid_rect.md)
to generate regional subsets.

## See also

[`grid_rect`](https://gcol33.github.io/hexify/reference/grid_rect.md)
for regional grids

## Examples

``` r
# Coarse global grid
grid <- hex_grid(area_km2 = 100000)
global <- grid_global(grid)
plot(global)
```
