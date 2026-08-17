# Plot hexagonized globe

Renders a global hexagonal grid on an orthographic projection with
customizable rotation, land clipping, and styling options.

## Usage

``` r
plot_globe(
  area = 50000,
  center = "europe",
  clip_to_land = FALSE,
  land_data = NULL,
  exclude_antarctica = TRUE,
  fill = "#D4B896",
  border = "grey30",
  border_width = 0.2,
  ocean_fill = "white",
  ocean_border = "grey50",
  show_land = clip_to_land,
  land_fill = NA,
  land_border = "grey40",
  land_width = 0.3,
  use_ggplot = NULL,
  return_data = FALSE,
  aperture = 3L
)
```

## Arguments

- area:

  Cell area in km^2 (passed to
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md))

- center:

  Globe center: either a preset name (e.g., "europe") or numeric vector
  c(lon, lat). See
  [`globe_centers`](https://gillescolling.com/hexify/reference/globe_centers.md)
  for presets.

- clip_to_land:

  If TRUE, clip hexagons to land boundaries

- land_data:

  Optional sf object for land boundaries. If NULL and clip_to_land is
  TRUE, uses rnaturalearth::ne_countries()

- exclude_antarctica:

  If TRUE, exclude Antarctica from land clipping

- fill:

  Fill color for hexagons (default "#D4B896")

- border:

  Border color for hexagons (default "grey30")

- border_width:

  Border width for hexagons (default 0.2)

- ocean_fill:

  Fill color for ocean/globe background (default "white")

- ocean_border:

  Border color for globe circle (default "grey50")

- show_land:

  If TRUE, show land boundaries (default TRUE when clipping)

- land_fill:

  Fill color for land (default NA, transparent)

- land_border:

  Border color for land boundaries (default "grey40")

- land_width:

  Border width for land boundaries (default 0.3)

- use_ggplot:

  NULL = auto-detect, TRUE = force ggplot2, FALSE = force base

- return_data:

  If TRUE, return sf objects instead of plotting

- aperture:

  Grid aperture (default 3L)

## Value

If use_ggplot = TRUE: ggplot2 object (can add layers with +) If
use_ggplot = FALSE: NULL invisibly (plots directly) If return_data =
TRUE: list of sf objects (hexagons, land, ocean_circle, crs)

## Details

The function handles several technical challenges:

- Hexagons on the back side of the globe fail to transform - these are
  filtered out gracefully

- Invalid geometries after projection are repaired with st_buffer(0)

- Clipping is done in orthographic CRS to avoid topology errors

## See also

[`globe_centers`](https://gillescolling.com/hexify/reference/globe_centers.md)
for available presets,
[`grid_global`](https://gillescolling.com/hexify/reference/grid_global.md)
for generating global grids without plotting

## Examples

``` r
# \donttest{
# Get data for custom plotting (no rendering)
data <- plot_globe(area = 100000, center = "europe", return_data = TRUE)
nrow(data$hexagons)
class(data$ocean_circle)

# Basic usage - Europe-centered globe
plot_globe(area = 80000, center = "europe")
# }
```
