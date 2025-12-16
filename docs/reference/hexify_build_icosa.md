# Initialize icosahedron geometry

Sets up the icosahedron state for ISEA projection. Uses standard ISEA3H
orientation by default (vertex 0 at 11.25E, 58.28N).

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

## Details

The icosahedron is initialized lazily at the C++ level when first
needed. Manual call is only required for non-standard orientations.

## See also

Other projection:
[`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md),
[`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md),
[`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md),
[`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md),
[`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md),
[`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md),
[`hexify_set_precision()`](https://gcol33.github.io/hexify/reference/hexify_set_precision.md),
[`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md),
[`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)

## Examples

``` r
# Use standard ISEA3H orientation
hexify_build_icosa()

# Custom orientation
hexify_build_icosa(vert0_lon = 0, vert0_lat = 90, azimuth = 0)
```
