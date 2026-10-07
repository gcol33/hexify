# Plot a grid on the sphere, on its solid, or on the unfolded net

Draws the cells of a grid in 3D, seen from above a chosen point: on the
sphere, or on the flat faces of the solid the grid is built on. All
surfaces take their cell boundaries from the same points, so a cell on
the solid is the cell on the sphere folded flat.

## Usage

``` r
# S4 method for class 'HexGridInfo,missing'
plot(
  x,
  y,
  surface = c("sphere", "solid", "net"),
  center = c(lon = 15, lat = 32),
  projection = c("orthographic", "perspective"),
  distance = NULL,
  tilt = 0,
  rotation = 0,
  fov = NULL,
  cells = NULL,
  land = TRUE,
  ocean_fill = "#F4F6F8",
  land_fill = "#C9D0D6",
  land_border = "#8C969F",
  land_lwd = 0.4,
  grid_border = "#0072B2",
  grid_lwd = 0.8,
  face_edges = NULL,
  edge_col = "black",
  edge_lwd = 1.1,
  layout = "plane",
  seams = FALSE,
  seam_col = "#D55E00",
  seam_lwd = 1.6,
  graticule = FALSE,
  graticule_col = "#9AA3AB",
  graticule_lwd = 0.5,
  distortion = c("none", "angular", "areal"),
  tissot = FALSE,
  tissot_col = "#2E3439",
  tissot_lwd = 0.8,
  parents = NULL,
  parent_col = "#1B2A38",
  parent_lwd = 1.8,
  area_legend = FALSE,
  tabs = FALSE,
  step = 0.01,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A HexGridInfo object from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)

- y:

  Ignored

- surface:

  `"sphere"`, `"solid"` or `"net"`. `"solid"` draws the flat faces of
  the grid's polyhedron (the icosahedron, or the octahedron of a grid
  built on it). The solid and the net need an ISEA grid; an H3 grid is
  drawn on the sphere.

- center:

  Point the view looks down on: a preset name from
  [`globe_centers`](https://gillescolling.com/hexify/reference/globe_centers.md)
  or `c(lon, lat)`.

- projection:

  `"orthographic"` or `"perspective"`.

- distance:

  Distance of the untilted perspective camera from the sphere's centre,
  in sphere radii; greater than 1. `NULL` uses 3. A geostationary
  satellite sits at about 6.6 Earth radii.

- tilt:

  Degrees the perspective camera swings about the point at `center`,
  between -90 and 90, keeping its range of `distance - 1` radii to that
  point. Positive values move the camera towards the bottom of the view,
  so it looks towards the top.

- rotation:

  Degrees the view turns about the line of sight; positive values turn
  north clockwise.

- fov:

  Field of view of the perspective camera in degrees, across the square
  the plot shows, centred on `center`; below 170. `NULL` frames the
  whole visible sphere instead.

- cells:

  Cell IDs to draw. `NULL` draws every cell of the grid.

- land:

  `TRUE` for the built-in
  [`hexify_world`](https://gillescolling.com/hexify/reference/hexify_world.md),
  `FALSE` for none, or an sf/sfc object of polygons.

- ocean_fill:

  Fill colour of the surface.

- land_fill:

  Fill colour of land; `NA` leaves land unfilled.

- land_border:

  Colour of country outlines; `NA` draws none.

- land_lwd:

  Line width of country outlines.

- grid_border:

  Colour of cell boundaries.

- grid_lwd:

  Line width of cell boundaries.

- face_edges:

  Draw the edges of the solid's faces. `NULL` draws them for an ISEA
  grid; H3 is built on an icosahedron of its own, so its grid takes
  none.

- edge_col:

  Colour of face edges.

- edge_lwd:

  Line width of face edges.

- layout:

  For `surface = "net"`, the layout of the faces: a name
  [`net_layout`](https://gillescolling.com/hexify/reference/net_layout.md)
  takes, or a `hexify_net` object it built for this grid.

- seams:

  For `surface = "net"`, draw the cuts of the net: the sides of its
  pieces across which the map is not continuous, including its outline.

- seam_col:

  Colour of seams.

- seam_lwd:

  Line width of seams.

- graticule:

  `FALSE` for none, `TRUE` for meridians and parallels every 30 degrees,
  or a number of degrees between them.

- graticule_col:

  Colour of the graticule.

- graticule_lwd:

  Line width of the graticule.

- distortion:

  `"none"`, `"angular"` for the maximum angular deformation of the face
  projection, or `"areal"` for its areal scale (1 everywhere on Snyder's
  equal-area projection). Needs an ISEA grid; a colour bar is drawn
  along the bottom.

- tissot:

  `FALSE` for none, `TRUE` for Tissot ellipses at points about 30
  degrees apart, or a number of degrees between them. Each ellipse is
  the image of a circle a fifth of that spacing in radius. Drawn on the
  net and the solid.

- tissot_col:

  Colour of Tissot ellipses.

- tissot_lwd:

  Line width of Tissot ellipses.

- parents:

  Resolutions up whose cell outlines are drawn over the cells, such as
  `1` or `1:2`; `NULL` for none.

- parent_col:

  Colour of parent outlines, recycled over `parents`.

- parent_lwd:

  Line width of parent outlines, recycled over `parents`.

- area_legend:

  For `surface = "net"` of a grid on Snyder's equal-area projection,
  draw one hexagonal cell of each drawn resolution at the scale of the
  map, labelled with its area.

- tabs:

  For `surface = "net"`, draw glue tabs on the seams.

- step:

  Spacing of the points along a cell boundary, as a fraction of a face
  edge.

- main:

  Plot title.

- ...:

  Ignored

## Value

The grid, invisibly

## Details

`surface = "net"` lays the faces out flat, by default in the PLANE
layout of DGGRID
([`hexify_cell_to_plane`](https://gillescolling.com/hexify/reference/hexify_cell_to_plane.md)
gives cell centres in the same coordinates); `layout` picks another,
such as Van de Sande's Gosper World on the octahedron (see
[`net_layout`](https://gillescolling.com/hexify/reference/net_layout.md)).
A cell on a cut of the net is drawn in parts, one on each side. The net
is drawn without a camera, so `center`, `projection`, `distance`,
`tilt`, `rotation` and `fov` apply to the other two surfaces.

`distortion` colours the surface by the distortion of the face
projection that maps the sphere onto the faces (see
[`projection_distortion`](https://gillescolling.com/hexify/reference/projection_distortion.md)),
and `tissot` draws Tissot's indicatrix on the faces: the ellipse each
small circle of the sphere maps to. With `graticule` on the net,
Snyder's projection shows small bends in the meridians and parallels
where they cross the arcs from each face's centre to its corners.

`parents` draws the outlines of the cells some resolutions up over the
grid's own cells, to show the hierarchy. On a net of an equal-area grid,
`area_legend` draws one cell at the scale of the map, with its area, and
the same for each resolution `parents` draws. `tabs` adds glue tabs
along the seams of a net, one on each pair of edges glued together, so a
printed net (for example through
[`pdf`](https://rdrr.io/r/grDevices/pdf.html)) folds into the solid
(Carr et al. 1997).

The view is an orthographic projection by default: parallel lines of
sight, so the whole near half of the sphere shows.
`projection = "perspective"` places a camera `distance` sphere radii
from the sphere's centre, above `center`; it sees a smaller cap, and
nearer cells appear larger. `tilt` swings that camera about the point at
`center`, which it keeps looking at, so the surface is seen at a slant
with the horizon towards the top; by default the plot frames the whole
visible sphere. `fov` sets how wide the camera sees, and `rotation`
turns the view about the line of sight in either projection.

## References

Carr, D. B., Kahn, R., Sahr, K., Olsen, A. R. (1997). ISEA discrete
global grids. Statistical Computing & Graphics Newsletter 8(2/3): 31-39.

## See also

[`grid_global`](https://gillescolling.com/hexify/reference/grid_global.md)
for the cells as sf polygons,
[`projection_distortion`](https://gillescolling.com/hexify/reference/projection_distortion.md),
[`net_cells`](https://gillescolling.com/hexify/reference/net_cells.md)

## Examples

``` r
grid <- hex_grid(resolution = 3, aperture = 3)
plot(grid)
plot(grid, surface = "solid")
plot(grid, surface = "solid", land = FALSE, center = "pacific")
plot(grid, projection = "perspective", distance = 2, center = "europe")
plot(grid, surface = "solid", projection = "perspective",
     distance = 2.5, tilt = 25, rotation = 30)
plot(grid, surface = "net")
octa <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
                 orientation = "gosper")
plot(octa, surface = "net", layout = "gosper", seams = TRUE, graticule = TRUE)

# Distortion of the face projections, with Tissot's indicatrix
plot(grid, surface = "net", land = FALSE, distortion = "angular",
     tissot = TRUE, graticule = 15, cells = numeric(0))
fuller <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
plot(fuller, distortion = "areal", land = FALSE, cells = numeric(0))

# Two resolutions with an area legend, and a net to print and fold
fine <- hex_grid(resolution = 4, aperture = 3)
plot(fine, surface = "net", land = FALSE, parents = 1, area_legend = TRUE)
plot(grid, surface = "net", layout = "land", land = FALSE, seams = TRUE,
     tabs = TRUE)
```
