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

Cell IDs: a
[`bit64::integer64`](https://bit64.r-lib.org/reference/bit64-package.html)
vector for an ISEA grid, a character vector of H3 indices for an H3 grid

## Details

ISEA cell IDs number a grid's cells from 1 and pass 2^53, the largest
whole number a double holds exactly, at fine resolutions, so they are
returned as 64-bit integers. Functions taking cell IDs also accept whole
numbers below 2^53 and character strings of digits.

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
