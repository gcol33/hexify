# Validate 'dggridR' Grid Compatibility with hexify

Checks whether a 'dggridR' grid object is compatible with hexify
functions. Returns TRUE if compatible, or throws an error describing
incompatibilities.

## Usage

``` r
dggrid_is_compatible(dggs, strict = TRUE)
```

## Arguments

- dggs:

  A 'dggridR' grid object

- strict:

  If TRUE (default), throw errors for incompatibilities. If FALSE,
  return FALSE instead of throwing errors.

## Value

TRUE if compatible, FALSE if not compatible (when strict=FALSE)

## See also

Other 'dggridR' compatibility:
[`as_dggrid()`](https://gillescolling.com/hexify/reference/as_dggrid.md),
[`dggrid_43h_sequence()`](https://gillescolling.com/hexify/reference/dggrid_43h_sequence.md),
[`from_dggrid()`](https://gillescolling.com/hexify/reference/from_dggrid.md)

## Examples

``` r
dggrid_is_compatible(as_dggrid(hex_grid(resolution = 4, aperture = 7)))
dggrid_is_compatible(list(projection = "ISEA", topology = "DIAMOND",
                          aperture = 4, res = 4), strict = FALSE)
```
