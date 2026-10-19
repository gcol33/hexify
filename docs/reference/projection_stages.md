# Points Through the Steps of Snyder's Projection

Snyder's equal-area projection of a face is Lambert's azimuthal
equal-area projection about the face's centre, adjusted to fill the
face's plane triangle. Lambert's projection has a construction with
compass and ruler, and this function returns each point at every step of
it, on one face.

## Usage

``` r
projection_stages(x, lon, lat, face = NULL)
```

## Arguments

- x:

  A HexGridInfo object from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)
  with Snyder's projection (`projection = "isea"`).

- lon, lat:

  Longitudes and latitudes in degrees (geodetic on the grid's ellipsoid,
  if it has one: the construction starts from the point's place on the
  authalic sphere).

- face:

  `NULL` to take each point onto the face it lies on, or face numbers
  from 0, one or one per point; the face's formulas are then applied as
  they stand to points off it.

## Value

A data frame with one row per point: `lon`, `lat`, `face` (from 0);
`arc`, \\z\\ in degrees; `az` (\\Az\\) and `az_prime` (\\Az'\\) in
degrees from the face's first vertex, clockwise seen from outside the
sphere; Snyder's factor `f`; the plane coordinates `sphere_u`,
`sphere_v`, `lambert_u`, `lambert_v`, `nudge_u`, `nudge_v`, `snyder_u`,
`snyder_v`; and `tx`, `ty`, Snyder's point in triangle coordinates (a
face edge is 1). The attribute `"r1"` holds \\R'\\.

## Details

With \\T\\ the face centre on the sphere and \\S\\ a point at arc \\z\\
from it and azimuth \\Az\\:

- sphere:

  \\S\\ itself.

- lambert:

  \\S\\ swings about \\T\\, in the plane of \\S\\ and the normal at
  \\T\\, down onto the tangent plane at \\T\\, keeping its distance
  \\\|TS\| = 2 \sin(z / 2)\\: the Lambert point \\L\\, at azimuth
  \\Az\\. This map is equal-area and takes the spherical face to a
  triangle with curved edges.

- nudge:

  Within the tangent plane, Snyder's adjustment moves \\L\\ to azimuth
  \\Az'\\ and radius \\2 f \sin(z / 2)\\, which fill each of the face's
  six right triangles in the plane with area in proportion (Snyder 1992,
  eqs. 6-11).

- snyder:

  The tangent plane scales by \\R'\\ about the sphere's centre onto the
  face plane, at distance \\R'\\ from it; the point lands on Snyder's
  point \\P\\, at radius \\2 R' f \sin(z / 2)\\, and the face on the
  plane triangle of the face's area.

Plane coordinates (`u`, `v`) run from the face centre along the x and y
axes of the face's plane triangle, the axes of `tx` and `ty`, in units
of the sphere's radius; an azimuth \\a\\ lies along \\(\sin a, \cos
a)\\. Seen from the sphere's centre the planes are parallel, and the
point's height along the face's normal is \\\cos z\\ for \\S\\, 1 on the
tangent plane and \\R'\\ on the face plane. The sphere point's `u` and
`v` are its projection onto the planes along the normal.

The scale factors at every step come from
[`projection_distortion`](https://gillescolling.com/hexify/reference/projection_distortion.md)
with its `stage` argument; they and these points are computed by the
same code as the forward projection, and `tx`, `ty` are those
[`hexify_forward`](https://gillescolling.com/hexify/reference/hexify_forward.md)
returns.

## References

Snyder, J. P. (1992). An equal-area map projection for polyhedral
globes. *Cartographica* 29(1), 10-21.
[doi:10.3138/27H7-8K88-4882-1752](https://doi.org/10.3138/27H7-8K88-4882-1752)

## See also

[`projection_distortion`](https://gillescolling.com/hexify/reference/projection_distortion.md);
[`hex_globe`](https://gillescolling.com/hexify/reference/hex_globe.md)
with `lambert` to follow the construction on a globe

## Examples

``` r
g <- hex_grid(resolution = 3, aperture = 3)
s <- projection_stages(g, c(0, 10, 16.37), c(40, 45, 48.21))
s

# The Lambert point keeps its distance from the face centre
cbind(lambert = sqrt(s$lambert_u^2 + s$lambert_v^2),
      chord = 2 * sin(s$arc * pi / 360))
```
