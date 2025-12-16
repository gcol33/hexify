# Convert Quad IJ coordinates to cell ID

Converts Quad IJ coordinates to a global cell identifier. This is the
final step in the coordinate pipeline.

## Usage

``` r
hexify_quad_ij_to_cell(quad, i, j, resolution, aperture = 3L)
```

## Arguments

- quad:

  Quad number (0-11), integer or vector

- i:

  Cell index along first axis, integer or vector

- j:

  Cell index along second axis, integer or vector

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Numeric vector of cell IDs

## See also

Other coordinate conversion:
[`hexify_cell_id_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_cell_id_to_quad_ij.md),
[`hexify_cell_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_icosa_tri.md),
[`hexify_cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md),
[`hexify_cell_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_plane.md),
[`hexify_cell_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_ij.md),
[`hexify_cell_to_quad_xy()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_xy.md),
[`hexify_grid_cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md),
[`hexify_grid_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md),
[`hexify_icosa_tri_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_plane.md),
[`hexify_icosa_tri_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_quad_ij.md),
[`hexify_icosa_tri_to_quad_xy()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_quad_xy.md),
[`hexify_lonlat_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md),
[`hexify_lonlat_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_plane.md),
[`hexify_lonlat_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_quad_ij.md),
[`hexify_quad_ij_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_icosa_tri.md),
[`hexify_quad_ij_to_xy()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_xy.md),
[`hexify_quad_xy_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md),
[`hexify_quad_xy_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_icosa_tri.md),
[`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Convert Quad IJ to cell ID
cell_id <- hexify_quad_ij_to_cell(quad = 1, i = 100, j = 50,
                                   resolution = 10, aperture = 3)
print(cell_id)
} # }
```
