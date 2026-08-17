# Verify grid object

Validates that a grid object has all required fields and valid values.
This function is called internally by most hexify functions to ensure
grid integrity.

## Usage

``` r
dgverify(dggs)
```

## Arguments

- dggs:

  Grid object to verify (from hexify_grid() or hex_grid())

## Value

TRUE (invisibly) if valid, otherwise throws an error

## Examples

``` r
grid <- hexify_grid(area = 1000, aperture = 3)
dgverify(grid)  # Should pass silently

# Modern HexGridInfo objects are accepted too
dgverify(hex_grid(area_km2 = 1000))

# Invalid grid will throw error
bad_grid <- list(aperture = 5)
try(dgverify(bad_grid))  # Will error
```
