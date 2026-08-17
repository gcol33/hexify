# Convert index string to cell coordinates

Decodes a hierarchical index string back to its cell coordinates and
resolution. For Z7, valid indices round-trip exactly through
[`hexify_cell_to_index()`](https://gillescolling.com/hexify/reference/hexify_cell_to_index.md).

## Usage

``` r
hexify_index_to_cell(
  index,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- index:

  Index string

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding used. Default "auto" infers from aperture.

## Value

A list with `face`, `i`, `j`, and `resolution`.

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gillescolling.com/hexify/reference/hexify_cell_to_index.md),
[`hexify_compare_indices()`](https://gillescolling.com/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)

## Examples

``` r
cell <- hexify_index_to_cell("0012012", aperture = 3)

z7_cell <- hexify_index_to_cell("110001", aperture = 7)
hexify_cell_to_index(z7_cell$face, z7_cell$i, z7_cell$j,
  z7_cell$resolution, aperture = 7
)
```
