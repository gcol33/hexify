# Extract plain data frame from hexify result

Convenience function to get a plain data.frame from hexify() result.
Alias for
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## Usage

``` r
hexify_df(x)
```

## Arguments

- x:

  Result from
  [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)

## Value

A data.frame

## Examples

``` r
if (FALSE) { # \dontrun{
result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
df_plain <- hexify_df(result)
} # }
```
