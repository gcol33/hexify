# Extract grid from various objects

Internal function to extract a HexGridInfo from different input types.
Accepts HexGridInfo, HexData, or legacy hexify_grid objects.

## Usage

``` r
extract_grid(x, allow_null = FALSE)
```

## Arguments

- x:

  Object containing grid info

- allow_null:

  If TRUE, return NULL when x is NULL

## Value

HexGridInfo object
