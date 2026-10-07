# Distortion of a grid's face projection

Tissot's indicatrix of the projection that maps the sphere onto the
faces of the grid's solid (Snyder's equal-area projection, or Fuller's):
a small circle on the sphere maps to an ellipse on the face, with
semi-axes `a` and `b` times the circle's radius.

## Usage

``` r
projection_distortion(x, lon, lat)
```

## Arguments

- x:

  A HexGridInfo object from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)
  with an ISEA grid. H3 is built on a gnomonic projection of its own and
  is not covered.

- lon, lat:

  Longitudes and latitudes in degrees.

## Value

A data frame with one row per point: `lon`, `lat`, `face` (from 0), the
scale factors `a` (largest) and `b` (smallest), the maximum angular
deformation `angular` \\2 \arcsin((a - b) / (a + b))\\ in degrees, the
areal scale `areal` \\a b\\, and `angle`, the direction of the longer
axis on the face plane in degrees anticlockwise from the face's x axis.

## Details

The scale factors are the singular values of the projection's
derivative, computed exactly by forward-mode automatic differentiation
of the projection code. Face-plane lengths are measured on the plane
triangle that has the face's area, so the areal scale `a * b` is 1
everywhere under Snyder's projection, and its mean over a face is 1
under Fuller's.

Snyder's projection bends along the arcs from each face's centre to its
corners, where its derivative jumps; there and on face edges the result
is the derivative on one side. Its angular deformation is largest at the
face centre, approached along those arcs: 17.27 degrees, with scale
factors 1.163 and 0.860 (Snyder 1992, Table 1). At a face centre the
projection has a different derivative along each direction, and the
result is the limit along the azimuth the point's coordinates give.

## References

Snyder, J. P. (1992). An equal-area map projection for polyhedral
globes. *Cartographica* 29(1), 10-21.
[doi:10.3138/27H7-8K88-4882-1752](https://doi.org/10.3138/27H7-8K88-4882-1752)

van Leeuwen, D. and Strebe, D. (2006). A "slice-and-dice" approach to
area equivalence in polyhedral map projections. *Cartography and
Geographic Information Science* 33(4), 269-286.
[doi:10.1559/152304006779500687](https://doi.org/10.1559/152304006779500687)

## See also

[`plot`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)
with `distortion` and `tissot` to map it

## Examples

``` r
g <- hex_grid(resolution = 3, aperture = 3)
projection_distortion(g, c(0, 16.37), c(0, 48.21))

gf <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
d <- projection_distortion(gf, runif(1000, -180, 180),
                           asin(runif(1000, -1, 1)) * 180 / pi)
range(d$areal)
```
