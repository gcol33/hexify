# Convert Icosa Triangle to Quad XY coordinates

Converts icosahedral triangle coordinates (from Snyder projection) to
quad XY coordinates. This is an intermediate step in the pipeline.

## Usage

``` r
hexify_icosa_tri_to_quad_xy(
  icosa_triangle_face,
  icosa_triangle_x,
  icosa_triangle_y
)
```

## Arguments

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Value

List with components:

- quad:

  Quad number (0-11)

- quad_x:

  Continuous X coordinate in quad space

- quad_y:

  Continuous Y coordinate in quad space

## See also

Other coordinate conversion:
[`hexify_cell_id_to_quad_ij()`](https://gillescolling.com/hexify/reference/hexify_cell_id_to_quad_ij.md),
[`hexify_cell_to_icosa_tri()`](https://gillescolling.com/hexify/reference/hexify_cell_to_icosa_tri.md),
[`hexify_cell_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_cell_to_lonlat.md),
[`hexify_cell_to_plane()`](https://gillescolling.com/hexify/reference/hexify_cell_to_plane.md),
[`hexify_cell_to_quad_ij()`](https://gillescolling.com/hexify/reference/hexify_cell_to_quad_ij.md),
[`hexify_cell_to_quad_xy()`](https://gillescolling.com/hexify/reference/hexify_cell_to_quad_xy.md),
[`hexify_grid_cell_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_grid_cell_to_lonlat.md),
[`hexify_grid_to_cell()`](https://gillescolling.com/hexify/reference/hexify_grid_to_cell.md),
[`hexify_icosa_tri_to_plane()`](https://gillescolling.com/hexify/reference/hexify_icosa_tri_to_plane.md),
[`hexify_icosa_tri_to_quad_ij()`](https://gillescolling.com/hexify/reference/hexify_icosa_tri_to_quad_ij.md),
[`hexify_lonlat_to_cell()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_cell.md),
[`hexify_lonlat_to_plane()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_plane.md),
[`hexify_lonlat_to_quad_ij()`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_quad_ij.md),
[`hexify_quad_ij_to_cell()`](https://gillescolling.com/hexify/reference/hexify_quad_ij_to_cell.md),
[`hexify_quad_ij_to_icosa_tri()`](https://gillescolling.com/hexify/reference/hexify_quad_ij_to_icosa_tri.md),
[`hexify_quad_ij_to_xy()`](https://gillescolling.com/hexify/reference/hexify_quad_ij_to_xy.md),
[`hexify_quad_xy_to_cell()`](https://gillescolling.com/hexify/reference/hexify_quad_xy_to_cell.md),
[`hexify_quad_xy_to_icosa_tri()`](https://gillescolling.com/hexify/reference/hexify_quad_xy_to_icosa_tri.md),
[`hexify_roundtrip_test()`](https://gillescolling.com/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
# First get triangle coordinates from lon/lat
fwd <- hexify_forward(lon = 2.35, lat = 48.86)

# Then convert to quad XY
quad_xy <- hexify_icosa_tri_to_quad_xy(
  icosa_triangle_face = fwd["face"],
  icosa_triangle_x = fwd["icosa_triangle_x"],
  icosa_triangle_y = fwd["icosa_triangle_y"]
)
print(quad_xy)
```
