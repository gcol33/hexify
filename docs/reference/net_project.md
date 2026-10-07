# Place points on a net layout

Projects longitude/latitude points onto their face of the grid's solid
and places them where the layout puts that part of the face. A point on
a piece the layout shows more than once, or on a cut between two pieces,
has one row per place.

## Usage

``` r
net_project(layout, lon, lat)
```

## Arguments

- layout:

  A `hexify_net` object from
  [`net_layout`](https://gillescolling.com/hexify/reference/net_layout.md)

- lon, lat:

  Longitudes and latitudes in degrees

## Value

A data frame with columns `point` (position in `lon`), `piece` (position
in `layout$pieces`), `x` and `y`, in units of a face edge.

## Examples

``` r
grid <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
                 orientation = "gosper")
world <- net_layout(grid, "gosper")
net_project(world, c(16.37, -74.0), c(48.21, 40.71))
```
