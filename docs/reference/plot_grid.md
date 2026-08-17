# Plot hexagonal grid clipped to a polygon boundary

A convenience function that creates a grid, clips it to a boundary
polygon, and plots the result in a single call.

## Usage

``` r
plot_grid(
  boundary,
  grid,
  crop = TRUE,
  boundary_fill = "gray95",
  boundary_border = "gray40",
  boundary_lwd = 0.5,
  grid_fill = "steelblue",
  grid_border = "steelblue",
  grid_lwd = 0.3,
  grid_alpha = 0.3,
  title = NULL,
  expand = 0.05
)
```

## Arguments

- boundary:

  An sf/sfc polygon to clip to (e.g., country boundary)

- grid:

  A HexGridInfo object from
  [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md)

- crop:

  If TRUE (default), cells are cropped to boundary. If FALSE, only
  complete hexagons within boundary are shown.

- boundary_fill:

  Fill color for the boundary polygon (default "gray95")

- boundary_border:

  Border color for boundary (default "gray40")

- boundary_lwd:

  Line width for boundary (default 0.5)

- grid_fill:

  Fill color for grid cells (default "steelblue")

- grid_border:

  Border color for grid cells (default "steelblue")

- grid_lwd:

  Line width for cell borders (default 0.3)

- grid_alpha:

  Transparency for cell fill (0-1, default 0.3)

- title:

  Plot title. If NULL (default), auto-generates title with cell area.

- expand:

  Expansion factor for plot limits (default 0.05)

## Value

A ggplot object that can be further customized

## Details

This is a convenience wrapper around
[`grid_clip()`](https://gillescolling.com/hexify/reference/grid_clip.md)
that handles the common use case of visualizing a hexagonal grid over a
geographic region.

## See also

[`grid_clip`](https://gillescolling.com/hexify/reference/grid_clip.md)
for the underlying clipping function,
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
grid specification

## Examples

``` r
# \donttest{
# Plot grid over France
france <- hexify_world[hexify_world$name == "France", ]
grid <- hex_grid(area_km2 = 2000)
plot_grid(france, grid)

# Customize colors
plot_grid(france, grid,
          grid_fill = "coral", grid_alpha = 0.5,
          boundary_fill = "lightyellow")

# Keep only complete hexagons
plot_grid(france, grid, crop = FALSE)

# Add ggplot2 customizations
library(ggplot2)
plot_grid(france, grid) +
  labs(subtitle = "ISEA3H Discrete Global Grid") +
  theme_void()
# }
```
