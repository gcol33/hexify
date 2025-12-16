# Convert longitude/latitude to PLANE coordinates

Converts geographic coordinates directly to PLANE coordinates (unfolded
icosahedron). Combines forward Snyder projection with the PLANE
transformation.

## Usage

``` r
hexify_lonlat_to_plane(lon, lat)
```

## Arguments

- lon:

  Longitude in degrees (-180 to 180)

- lat:

  Latitude in degrees (-90 to 90)

## Value

Data frame with columns:

- plane_x:

  X coordinate in PLANE space (range ~0 to 5.5)

- plane_y:

  Y coordinate in PLANE space (range ~0 to 1.73)

## Details

Equivalent to dggridR's dgGEO_to_PLANE().

## See also

[`hexify_cell_to_plane`](https://gcol33.github.io/hexify/reference/hexify_cell_to_plane.md)
for cell ID conversion,
[`hexify_icosa_tri_to_plane`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_plane.md)
for triangle conversion

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
[`hexify_lonlat_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_quad_ij.md),
[`hexify_quad_ij_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md),
[`hexify_quad_ij_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_icosa_tri.md),
[`hexify_quad_ij_to_xy()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_xy.md),
[`hexify_quad_xy_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md),
[`hexify_quad_xy_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_icosa_tri.md),
[`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot world cities in PLANE coordinates
cities <- data.frame(
  lon = c(2.35, -74.00, 139.69, 151.21),
  lat = c(48.86, 40.71, 35.69, -33.87)
)
plane <- hexify_lonlat_to_plane(cities$lon, cities$lat)
plot(plane$plane_x, plane$plane_y)
} # }
```
