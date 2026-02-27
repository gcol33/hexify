# Get default index type for aperture

Returns the recommended index type:

- Aperture 3: "z3"

- Aperture 4: "zorder"

- Aperture 7: "z7"

## Usage

``` r
hexify_default_index_type(aperture)
```

## Arguments

- aperture:

  Aperture (3, 4, or 7)

## Value

String: "z3", "z7", or "zorder"

## See also

Other hierarchical index:
[`hexify_cell_to_index()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_index.md),
[`hexify_compare_indices()`](https://gcol33.github.io/hexify/reference/hexify_compare_indices.md),
[`hexify_get_children()`](https://gcol33.github.io/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gcol33.github.io/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gcol33.github.io/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gcol33.github.io/hexify/reference/hexify_z7_canonical.md)
