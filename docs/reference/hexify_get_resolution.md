# Get index resolution

Returns the resolution level encoded in an index string.

## Usage

``` r
hexify_get_resolution(
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

  Index encoding. Default "auto".

## Value

Integer resolution level

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_index.md),
[`hexify_compare_indices()`](https://gcol33.github.io/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gcol33.github.io/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gcol33.github.io/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gcol33.github.io/hexify/reference/hexify_get_parent.md),
[`hexify_index_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gcol33.github.io/hexify/reference/hexify_z7_canonical.md)
