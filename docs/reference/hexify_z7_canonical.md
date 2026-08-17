# Get canonical form of Z7 index

Decodes and re-encodes a Z7 index until it reaches a stable form.
Current Z7 indices are bijective, so every valid index is already
canonical and this function normally returns its input unchanged. It
remains available for validating or normalizing indices created by older
hexify versions.

## Usage

``` r
hexify_z7_canonical(index, max_iterations = 128L)
```

## Arguments

- index:

  A length-one Z7 index string.

- max_iterations:

  Maximum number of decode/encode iterations. This is a safety bound for
  legacy indices; the default is 128.

## Value

A length-one character string containing the stable index.

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gillescolling.com/hexify/reference/hexify_cell_to_index.md),
[`hexify_compare_indices()`](https://gillescolling.com/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gillescolling.com/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_index.md)

## Examples

``` r
# Valid Z7 indices are stable
hexify_z7_canonical("110001")

cell <- hexify_index_to_cell("110001", aperture = 7)
identical(
  hexify_cell_to_index(cell$face, cell$i, cell$j,
    cell$resolution, aperture = 7
  ),
  "110001"
)
```
