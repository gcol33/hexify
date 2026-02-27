# Convert Cell ID to Icosa Triangle coordinates

Converts DGGRID-compatible cell IDs to icosahedral triangle coordinates
(face, x, y). These are the coordinates produced by the 'Snyder' 'ISEA'
forward projection.

## Usage

``` r
hexify_cell_to_icosa_tri(cell_id, resolution, aperture = 3L)
```

## Arguments

- cell_id:

  Numeric vector of cell IDs (1-based)

- resolution:

  Grid resolution level (0-30)

- aperture:

  Grid aperture: 3, 4, or 7

## Value

Data frame with columns:

- icosa_triangle_face:

  Triangle face number (0-19)

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Details

Compatible with 'dggridR' dgSEQNUM_to_PROJTRI().

## See also

[`hexify_cell_to_quad_ij`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_ij.md)
for conversion to Quad IJ,
[`hexify_cell_to_lonlat`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md)
for conversion to lon/lat

Other coordinate conversion:
[`hexify_cell_id_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_cell_id_to_quad_ij.md),
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
[`hexify_quad_xy_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_icosa_tri.md),
[`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
# Get triangle coordinates for a cell
result <- hexify_cell_to_icosa_tri(cell_id = 1000, resolution = 10, aperture = 3)
print(result)

# Convert back to lon/lat using inverse projection
coords <- hexify_inverse(result$icosa_triangle_face,
                         result$icosa_triangle_x,
                         result$icosa_triangle_y)
```
