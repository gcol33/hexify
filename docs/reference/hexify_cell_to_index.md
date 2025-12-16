# Convert cell coordinates to index string

Converts DGGRID cell coordinates (face, i, j) to a hierarchical index
string. The index type is automatically selected based on aperture
unless specified.

## Usage

``` r
hexify_cell_to_index(
  face,
  i,
  j,
  resolution,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- face:

  Face/quad number (0-19)

- i:

  I coordinate

- j:

  J coordinate

- resolution:

  Resolution level

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding: "auto" (default), "z3", "z7", or "zorder"

## Value

Index string (e.g., "051223")

## Details

Default index types by aperture:

- Aperture 3: Z3 (optimized digit selection)

- Aperture 4: Z-order (Morton curve)

- Aperture 7: Z7 (hierarchical with Class III handling)

## See also

Other hierarchical index:
[`hexify_compare_indices()`](https://gcol33.github.io/hexify/reference/hexify_compare_indices.md),
[`hexify_default_index_type()`](https://gcol33.github.io/hexify/reference/hexify_default_index_type.md),
[`hexify_get_children()`](https://gcol33.github.io/hexify/reference/hexify_get_children.md),
[`hexify_get_parent()`](https://gcol33.github.io/hexify/reference/hexify_get_parent.md),
[`hexify_get_resolution()`](https://gcol33.github.io/hexify/reference/hexify_get_resolution.md),
[`hexify_index_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_index_to_cell.md),
[`hexify_index_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_index_to_lonlat.md),
[`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md),
[`hexify_lonlat_to_index()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_index.md),
[`hexify_z7_canonical()`](https://gcol33.github.io/hexify/reference/hexify_z7_canonical.md)

## Examples

``` r
idx <- hexify_cell_to_index(5, 10, 15, resolution = 3, aperture = 3)
```
