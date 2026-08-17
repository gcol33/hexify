# Convert longitude/latitude to cell ID

Converts geographic coordinates to DGGS cell IDs using a grid
specification.

## Usage

``` r
lonlat_to_cell(lon, lat, grid)
```

## Arguments

- lon:

  Numeric vector of longitudes in degrees

- lat:

  Numeric vector of latitudes in degrees

- grid:

  A HexGridInfo or HexData object, or legacy hexify_grid

## Value

Numeric vector of cell IDs

## Details

This function accepts either a HexGridInfo object from
[`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md)
or a HexData object from
[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md). If a
HexData object is provided, its grid specification is extracted
automatically.

## See also

[`cell_to_lonlat`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
for the inverse operation,
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
creating grid specifications

## Examples

``` r
grid <- hex_grid(area_km2 = 1000)
cells <- lonlat_to_cell(lon = c(0, 10), lat = c(45, 50), grid = grid)

# Or use HexData object
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
cells <- lonlat_to_cell(lon = 5, lat = 48, grid = result)
```
