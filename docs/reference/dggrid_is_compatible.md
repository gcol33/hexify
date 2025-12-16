# Validate dggridR grid compatibility with hexify

Checks whether a dggridR grid object is compatible with hexify
functions. Returns TRUE if compatible, or throws an error describing
incompatibilities.

## Usage

``` r
dggrid_is_compatible(dggs, strict = TRUE)
```

## Arguments

- dggs:

  A dggridR grid object

- strict:

  If TRUE (default), throw errors for incompatibilities. If FALSE,
  return FALSE instead of throwing errors.

## Value

TRUE if compatible, FALSE if not compatible (when strict=FALSE)

## See also

Other dggridR compatibility:
[`as_dggrid()`](https://gcol33.github.io/hexify/reference/as_dggrid.md),
[`dggrid_43h_sequence()`](https://gcol33.github.io/hexify/reference/dggrid_43h_sequence.md),
[`from_dggrid()`](https://gcol33.github.io/hexify/reference/from_dggrid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dggridR)

# Compatible grid
dggs <- dgconstruct(res = 5, aperture = 3, projection = "ISEA")
dggrid_is_compatible(dggs)  # TRUE

# Incompatible grid (FULLER projection)
dggs2 <- dgconstruct(res = 5, projection = "FULLER")
dggrid_is_compatible(dggs2, strict = FALSE)  # FALSE
} # }
```
