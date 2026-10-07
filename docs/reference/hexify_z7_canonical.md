# Get canonical form of Z7 index

Decodes and re-encodes a Z7 index until it reaches a stable form. Every
valid index is already canonical and comes back unchanged. A string
whose first nonzero digit is the direction its pentagon base cell lacks
names no cell; it decodes, as in DGGRID, to the lattice point its digits
reach, and this function returns the index of the cell there.

## Usage

``` r
hexify_z7_canonical(index, max_iterations = 128L)
```

## Arguments

- index:

  Character vector of Z7 index strings.

- max_iterations:

  Maximum number of decode/encode iterations; the default is 128.

## Value

A character vector of stable indices.

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
