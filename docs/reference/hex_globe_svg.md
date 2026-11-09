# Save a globe as an SVG file

Draws a globe from
[`hex_globe`](https://gillescolling.com/hexify/reference/hex_globe.md)
as an SVG drawing, in the view the widget opens with: the ocean, the
land, the cells and every line as vector paths, placed and projected as
the widget's shader places them. The drawing is made on the CPU by the
'hexglobe' package and needs no graphics card or browser.

## Usage

``` r
hex_globe_svg(widget, file, width = 800, height = 800, background = NULL)
```

## Arguments

- widget:

  A globe from
  [`hex_globe`](https://gillescolling.com/hexify/reference/hex_globe.md).

- file:

  Path of the SVG file to write.

- width, height:

  Size of the drawing in pixels.

- background:

  Colour behind the globe, or `NULL` for none.

## Value

`file`, invisibly.

## Details

The cells of an ISEA grid, which the widget finds per pixel, are drawn
as one path per colour and their borders as polylines, so the file grows
with the number of cells drawn. Borders fade out where the cells shrink
below a few pixels, as in the widget. Values blended across cells
(`smooth = TRUE`) have no vector form and are refused;
[`hex_globe_png`](https://gillescolling.com/hexify/reference/hex_globe_png.md)
draws them.

## See also

[`hex_globe_png`](https://gillescolling.com/hexify/reference/hex_globe_png.md)
for the same view as a PNG image

## Examples

``` r
if (requireNamespace("hexglobe", quietly = TRUE)) {
  grid <- hex_grid(resolution = 2, aperture = 3)
  globe <- hex_globe(grid, land = FALSE)
  hex_globe_svg(globe, tempfile(fileext = ".svg"), width = 300, height = 300)
}
```
