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

## Examples

``` r
# Use standard ISEA3H orientation
hexify_build_icosa()

# Custom orientation
hexify_build_icosa(vert0_lon = 0, vert0_lat = 90, azimuth = 0)
```
