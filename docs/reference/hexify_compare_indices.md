# Compare two indices

Lexicographic comparison of two index strings.

## Usage

``` r
hexify_compare_indices(idx1, idx2)
```

## Arguments

- idx1:

  First index strings

- idx2:

  Second index strings

## Value

Integer vector: -1 where idx1 \< idx2, 0 where equal, 1 where idx1 \>
idx2. The two are read in step, and either may be length one.

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gillescolling.com/hexify/reference/hexify_cell_to_index.md),
[`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gillescolling.com/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)
