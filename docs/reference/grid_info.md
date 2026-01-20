# Get Grid Specification

Extract the grid specification from a HexData object.

## Usage

``` r
grid_info(x)
```

## Arguments

- x:

  A HexData object

## Value

A HexGridInfo object

## Examples

``` r
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
grid_spec <- grid_info(result)
```
