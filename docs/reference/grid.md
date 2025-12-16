# Get Grid Specification

Extract the grid specification from a HexData object.

## Usage

``` r
grid(x)
```

## Arguments

- x:

  A HexData object

## Value

A HexGrid object

## Examples

``` r
if (FALSE) { # \dontrun{
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
grid_spec <- grid(result)
} # }
```
