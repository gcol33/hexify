# Convert longitude/latitude to Quad IJ coordinates

Converts geographic coordinates to the intermediate Quad IJ
representation used internally by ISEA DGGS. Returns the quad number
(0-11) and integer cell indices (i, j) within that quad.

## Usage

``` r
hexify_lonlat_to_quad_ij(lon, lat, resolution, aperture = 3L)
```

## Arguments

- lon:

  Longitude in degrees (-180 to 180)

- lat:

  Latitude in degrees (-90 to 90)

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

List with components:

- quad:

  Quad number (0-11)

- i:

  Integer cell index along first axis

- j:

  Integer cell index along second axis

- icosa_triangle_face:

  Source icosahedral face (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Details

The 20 icosahedral triangle faces are grouped into 12 quads:

- Quad 0: North polar region

- Quads 1-5: Upper hemisphere rhombi

- Quads 6-10: Lower hemisphere rhombi

- Quad 11: South polar region

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
[`hexify_quad_ij_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md),
[`hexify_quad_ij_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_icosa_tri.md),
[`hexify_quad_ij_to_xy()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_xy.md),
[`hexify_quad_xy_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md),
[`hexify_quad_xy_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_icosa_tri.md),
[`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
# Get Quad IJ coordinates for Paris
result <- hexify_lonlat_to_quad_ij(lon = 2.35, lat = 48.86,
                                    resolution = 10, aperture = 3)
print(result)
```
