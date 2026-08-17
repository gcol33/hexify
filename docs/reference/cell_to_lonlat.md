# Convert cell ID to longitude/latitude

Converts DGGS cell IDs back to geographic coordinates (cell centers).

## Usage

``` r
cell_to_lonlat(cell_id, grid)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs

- grid:

  A HexGridInfo or HexData object

## Value

Data frame with lon_deg and lat_deg columns

## See also

[`lonlat_to_cell`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md)
for the forward operation

## Examples

``` r
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(c(0, 10), c(45, 50), grid)
coords <- cell_to_lonlat(cells, grid)
```
