# Convert longitude/latitude to index string

Projects geographic coordinates to grid cells and returns their
hierarchical index strings. Inputs are vectorized over `lon` and `lat`.

## Usage

``` r
hexify_lonlat_to_index(
  lon,
  lat,
  resolution,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- lon:

  Longitude in degrees

- lat:

  Latitude in degrees

- resolution:

  Resolution level

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding: "auto" (default), "z3", "z7", or "zorder"

## Value

A character vector of index strings.

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
[`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)

## Examples

``` r
idx <- hexify_lonlat_to_index(16.37, 48.21, resolution = 5, aperture = 3)
idx7 <- hexify_lonlat_to_index(16.37, 48.21,
  resolution = 4, aperture = 7
)
```
