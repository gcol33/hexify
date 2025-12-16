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

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_index.md),
[`hexify_compare_indices()`](https://gcol33.github.io/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gcol33.github.io/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gcol33.github.io/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gcol33.github.io/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gcol33.github.io/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_index.md)

## Examples

``` r
# These all return the same canonical form
hexify_z7_canonical("110001")
hexify_z7_canonical("110002")
```
