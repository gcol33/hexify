# Convert Quad XY to Icosa Triangle coordinates

Inverse transformation from quad coordinates back to icosahedral
triangle coordinates. Useful for projecting cell centers back to
lon/lat.

## Usage

``` r
hexify_quad_xy_to_icosa_tri(quad, quad_x, quad_y)
```

## Arguments

- quad:

  Quad number (0-11)

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## Value

List with components:

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

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
[`hexify_quad_ij_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md),
[`hexify_quad_ij_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_icosa_tri.md),
[`hexify_quad_ij_to_xy()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_xy.md),
[`hexify_quad_xy_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md),
[`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
# Convert quad XY back to triangle coordinates
tri <- hexify_quad_xy_to_icosa_tri(quad = 1, quad_x = 0.5, quad_y = 0.3)
print(tri)
```
