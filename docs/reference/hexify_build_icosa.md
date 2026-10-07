# Set the default icosahedron orientation

Sets the orientation read by the functions that take no grid: the
projection functions in this family and the low-level conversions that
take a resolution and an aperture. The standard ISEA orientation (vertex
0 at 11.25E, 58.28N, azimuth 0) is the default until this is called. A
grid carries its own orientation (see the `orientation` argument of
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)),
which this does not change.

## Usage

``` r
hexify_build_icosa(
  vert0_lon = ISEA_VERT0_LON_DEG,
  vert0_lat = ISEA_VERT0_LAT_DEG,
  azimuth = ISEA_AZIMUTH_DEG
)
```

## Arguments

- vert0_lon:

  Vertex 0 longitude in degrees (default ISEA_VERT0_LON_DEG)

- vert0_lat:

  Vertex 0 latitude in degrees (default ISEA_VERT0_LAT_DEG)

- azimuth:

  Azimuth rotation in degrees (default ISEA_AZIMUTH_DEG)

## Value

Invisible NULL. Called for side effect.

## See also

Other projection:
[`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md),
[`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md),
[`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)

## Examples

``` r
# Use standard ISEA3H orientation
hexify_build_icosa()

# Custom orientation
hexify_build_icosa(vert0_lon = 0, vert0_lat = 90, azimuth = 0)
```
