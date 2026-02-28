# Generate a rectangular grid of hexagons

Creates hexagon polygons covering a rectangular geographic region. For
H3 grids, all cells that overlap the bounding box are included (not just
cells whose center falls inside), ensuring full spatial coverage.

## Usage

``` r
grid_rect(bbox, grid)
```

## Arguments

- bbox:

  Bounding box as c(xmin, ymin, xmax, ymax), or an sf/sfc object

- grid:

  A HexGridInfo object specifying the grid parameters

## Value

sf object with hexagon polygons

## See also

[`grid_global`](https://gcol33.github.io/hexify/reference/grid_global.md)
for global grids

## Examples

``` r
grid <- hex_grid(area_km2 = 5000)
europe <- grid_rect(c(-10, 35, 30, 60), grid)
plot(europe)
```
