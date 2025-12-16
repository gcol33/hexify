# Get canonical form of Z7 index

For Z7 indices that form cycles during decode/encode, returns the
lexicographically smallest index in the cycle. Provides stable unique
identifiers for aperture 7 grids.

## Usage

``` r
hexify_z7_canonical(index, max_iterations = 128L)
```

## Arguments

- index:

  Z7 index string

- max_iterations:

  Maximum iterations for cycle detection (default 128)

## Value

Canonical form (lexicographically smallest in cycle)

## Examples

``` r
# These all return the same canonical form
hexify_z7_canonical("110001")
hexify_z7_canonical("110002")
```
