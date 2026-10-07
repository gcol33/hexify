# Save a globe as a PNG image

Draws a globe made by
[`hex_globe`](https://gillescolling.com/hexify/reference/hex_globe.md)
with the same WebGPU shader and pipelines the widget uses, and saves the
view it opens with, without its controls, as a PNG file. The image is
drawn on the graphics card and read back from it, so it is the same on a
machine without a display. Pixels off the globe are transparent.

## Usage

``` r
hex_globe_png(
  widget,
  file,
  width = 800,
  height = 800,
  scale = 1,
  timeout = 30,
  renderer = c("wgpu", "chrome")
)
```

## Arguments

- widget:

  A globe from
  [`hex_globe`](https://gillescolling.com/hexify/reference/hex_globe.md).

- file:

  Path of the PNG file to write.

- width, height:

  Size of the view in CSS pixels.

- scale:

  Image pixels per CSS pixel; 2 draws the view at twice the resolution.

- timeout:

  Seconds to wait for the globe to be drawn in Chrome.

- renderer:

  `"wgpu"` or `"chrome"`, as above.

## Value

`file`, invisibly.

## Details

`renderer = "wgpu"` draws through 'wgpu', the 'Rust' implementation of
WebGPU, in the 'hexglobe' package, with no browser.
`renderer = "chrome"` draws in headless Chrome, which needs the
'chromote' and 'htmlwidgets' packages and a Chromium browser, such as
Chrome or Edge, that
[`chromote::find_chrome()`](https://rstudio.github.io/chromote/reference/find_chrome.html)
finds; the browser is started with its GPU, which WebGPU draws on.

## See also

[`plot`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)
for a static plot that needs no graphics card

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- hex_grid(resolution = 3, aperture = 3)
globe <- hex_globe(grid, surface = "solid", center = "pacific")
hex_globe_png(globe, tempfile(fileext = ".png"), scale = 2)
} # }
```
