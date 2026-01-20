# Convert cell IDs to sf polygons

Creates sf polygon geometries for hexagonal grid cells.

## Usage

``` r
cell_to_sf(cell_id = NULL, grid)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs. If NULL and x is HexData, uses cells from
  x.

- grid:

  A HexGridInfo or HexData object. If HexData and cell_id is NULL,
  polygons are generated for all cells in the data.

## Value

sf object with cell_id and geometry columns

## Details

When called with a HexData object and no cell_id argument, this function
generates polygons for all unique cells in the data, which is useful for
plotting.

## See also

[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
grid specifications,
[`as_sf`](https://gcol33.github.io/hexify/reference/as_sf.md) for
converting HexData to sf

## Examples

``` r
# From grid specification
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(c(0, 10, 20), c(45, 50, 55), grid)
polys <- cell_to_sf(cells, grid)

# From HexData (all cells)
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
polys <- cell_to_sf(grid = result)
```
