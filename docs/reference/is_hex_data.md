# Check If an Object Is a HexData

Check If an Object Is a HexData

## Usage

``` r
is_hex_data(x)
```

## Arguments

- x:

  Object to check

## Value

Logical

## Examples

``` r
df <- data.frame(lon = c(0, 10), lat = c(45, 50))
is_hex_data(hexify(df, lon = "lon", lat = "lat", area_km2 = 1000))
is_hex_data(df)
```
