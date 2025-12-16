# Inverse Snyder projection

Converts face plane coordinates back to geographic coordinates.

## Usage

``` r
hexify_inverse(x, y, face, tol = NULL, max_iters = NULL)
```

## Arguments

- x:

  X coordinate on face plane

- y:

  Y coordinate on face plane

- face:

  Face index (0-19)

- tol:

  Convergence tolerance (NULL for default)

- max_iters:

  Maximum iterations (NULL for default)

## Value

Named numeric vector: c(lon_deg, lat_deg)

## Examples

``` r
coords <- hexify_inverse(0.5, 0.3, face = 2)
```
