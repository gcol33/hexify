# Get Cell IDs

Extract the unique cell IDs present in a HexData object.

## Usage

``` r
cells(x)
```

## Arguments

- x:

  A HexData object

## Value

A vector of cell IDs

## Examples

``` r
df <- data.frame(lon = c(0, 0.01, 20), lat = c(45, 45, 55))
cells(hexify(df, lon = "lon", lat = "lat", area_km2 = 1000))
```
