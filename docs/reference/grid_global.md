# Generate a global hexagon grid

Creates hexagon polygons covering the entire Earth.

## Usage

``` r
grid_global(grid, wrap_dateline = TRUE)
```

## Arguments

- grid:

  A HexGridInfo object specifying the grid parameters

- wrap_dateline:

  Logical. If TRUE (default), calls
  [`sf::st_wrap_dateline()`](https://r-spatial.github.io/sf/reference/st_transform.html)
  to split antimeridian-crossing polygons, which flat maps and planar
  (GEOS) operations need. Set to FALSE for orthographic/globe
  projections where wrapping creates gaps, and for spherical (s2)
  operations: each cell then stays one polygon whose corners lie in
  -180..180, so a cell crossing the antimeridian spans a flat lon/lat
  map.

## Value

sf object with hexagon polygons

## Details

Every cell of the grid, one row each, in cell ID order for ISEA grids.
For large grids (many small cells), consider using
[`grid_rect()`](https://gillescolling.com/hexify/reference/grid_rect.md)
to generate regional subsets.

## See also

[`grid_rect`](https://gillescolling.com/hexify/reference/grid_rect.md)
for regional grids

## Examples

``` r
# Coarse global grid
grid <- hex_grid(area_km2 = 100000)
global <- grid_global(grid)
plot(global)
```
