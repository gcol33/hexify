# Generate a rectangular grid of hexagons

Creates hexagon polygons covering a rectangular geographic region.

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
if (FALSE) { # \dontrun{
grid <- hex_grid(area_km2 = 5000)
europe <- grid_rect(c(-10, 35, 30, 60), grid)
plot(europe)
} # }
```
