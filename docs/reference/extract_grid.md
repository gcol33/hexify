# Extract grid from various objects

Internal function to extract a HexGrid from different input types.
Accepts HexGrid, HexData, or legacy hexify_grid objects.

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

HexGrid object
