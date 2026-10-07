# Interactive globe of a grid

Draws the cells of a grid on a globe that turns under the mouse,
rendered on the graphics card through WebGPU. A slider folds the solid
the grid is built on into the sphere, and the shader blends the two
surfaces. Values given per cell fill the cells through a colour ramp.

## Usage

``` r
hex_globe(
  x,
  values = NULL,
  cells = NULL,
  surface = c("sphere", "solid"),
  center = c(lon = 15, lat = 32),
  projection = c("orthographic", "perspective"),
  distance = NULL,
  tilt = 0,
  rotation = 0,
  fov = NULL,
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
  palette = "viridis",
  limits = NULL,
  na_fill = NA,
  step = 0.01,
  width = NULL,
  height = NULL,
  elementId = NULL
)
```

## Arguments

- x:

  A HexGridInfo object from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)

- values:

  Numeric values, one per cell in `cells`, filling the cells. `NULL`
  draws no fill.

- cells:

  Cell IDs to draw. `NULL` draws every cell of the grid.

- surface:

  Surface shown first, `"sphere"` or `"solid"`, the flat faces of the
  grid's polyhedron; the slider folds one into the other. An H3 grid is
  drawn on the sphere only.

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

- palette:

  Colours of the ramp from low to high values, or the name of a palette
  of [`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html).

- limits:

  Values at the two ends of the ramp; `NULL` uses the range of `values`.
  Values outside are drawn in the end colours.

- na_fill:

  Fill colour of cells whose value is `NA`; `NA` leaves them unfilled.

- step:

  Spacing of the points along a cell boundary, as a fraction of a face
  edge.

- width, height:

  Size of the widget, as CSS units or pixels.

- elementId:

  Id of the widget's HTML element.

## Value

An htmlwidget.

## Details

The cells of an ISEA grid are found per pixel on the graphics card, by
the same projection and cell numbering as
[`lonlat_to_cell`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md),
so the page holds the grid's description and the values rather than the
cells' outlines, and a grid of any resolution draws as fast as a coarse
one. Borders keep their width at every zoom and fade out where the cells
get too small to see. The pointer shows the ID and value of the cell
under it. H3 cells are drawn from their outlines.

Drag to turn the globe. Drag with Shift held to turn the view about the
line of sight and, in the perspective view, to tilt the camera. The
mouse wheel moves the perspective camera closer or zooms the
orthographic view. The view starts where the camera arguments, those of
the grid's
[`plot`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)
method, place it.

The globe needs the 'htmlwidgets' package and a viewer or browser with
WebGPU: Chrome and Edge (on Linux with some graphics cards only),
Firefox on Windows and macOS, Safari 26, and RStudio's viewer. Without
WebGPU the widget shows a notice instead. The article
<https://gillescolling.com/hexify/articles/globe.html> lists the
versions.

## See also

[`plot`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)
for the same view as a static plot

## Examples

``` r
if (requireNamespace("htmlwidgets", quietly = TRUE)) {
  grid <- hex_grid(resolution = 3, aperture = 3)
  hex_globe(grid)
  hex_globe(grid, surface = "solid", center = "pacific")

  cells <- seq_len(n_cells(grid))
  centres <- cell_to_lonlat(cells, grid)
  hex_globe(grid, values = centres$lat_deg, palette = "Blue-Red 3")
}
```
