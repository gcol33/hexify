# Convert Icosa Triangle coordinates to PLANE coordinates

Transforms icosahedral triangle coordinates to the 2D PLANE
representation (unfolded icosahedron). Each triangle is rotated and
translated to its position in the unfolded layout.

## Usage

``` r
hexify_icosa_tri_to_plane(
  icosa_triangle_face,
  icosa_triangle_x,
  icosa_triangle_y
)
```

## Arguments

- icosa_triangle_face:

  Triangle face number (0-19), integer or vector

- icosa_triangle_x:

  X coordinate on triangle face

- icosa_triangle_y:

  Y coordinate on triangle face

## Value

Data frame with columns:

- plane_x:

  X coordinate in PLANE space (range ~0 to 5.5)

- plane_y:

  Y coordinate in PLANE space (range ~0 to 1.73)

## Details

Equivalent to dggridR's dgPROJTRI_to_PLANE().

The PLANE layout arranges all 20 icosahedral faces into a roughly
rectangular region. Faces 0-4 and 5-9 form the upper row, while faces
10-14 and 15-19 form the lower row. Adjacent faces share edges in this
representation.

## See also

[`hexify_cell_to_plane`](https://gcol33.github.io/hexify/reference/hexify_cell_to_plane.md)
for direct cell ID conversion,
[`hexify_lonlat_to_plane`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_plane.md)
for lon/lat to PLANE

## Examples

``` r
if (FALSE) { # \dontrun{
# Get PLANE coordinates from triangle coordinates
fwd <- hexify_forward(lon = 2.35, lat = 48.86)
plane <- hexify_icosa_tri_to_plane(
  icosa_triangle_face = fwd$face,
  icosa_triangle_x = fwd$icosa_triangle_x,
  icosa_triangle_y = fwd$icosa_triangle_y
)
print(plane)
} # }
```
