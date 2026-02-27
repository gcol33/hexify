# Quick plot of hexify results

Simple plotting function for hexify output using base R graphics. For
ggplot2/tmap, use hexify_to_polygons() instead.

## Usage

``` r
hexify_plot(
  data,
  aperture = 3L,
  col = "lightgray",
  border = "black",
  add = FALSE,
  ...
)
```

## Arguments

- data:

  Data frame returned by hexify()

- aperture:

  Grid aperture (default 3)

- col:

  Fill color for hexagons

- border:

  Border color for hexagons

- add:

  If TRUE, add to existing plot

- ...:

  Additional arguments passed to polygon()

## Value

NULL (invisibly). Creates a plot as side effect.

## See also

[`hexify_map`](https://gcol33.github.io/hexify/reference/hexify_map.md)
for basemap support,
[`hexify_heatmap`](https://gcol33.github.io/hexify/reference/hexify_heatmap.md)
for ggplot2-based heatmaps

Other visualization:
[`hexify_heatmap()`](https://gcol33.github.io/hexify/reference/hexify_heatmap.md),
[`hexify_map()`](https://gcol33.github.io/hexify/reference/hexify_map.md),
[`plot_world()`](https://gcol33.github.io/hexify/reference/plot_world.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(hexify)

df <- data.frame(lon = c(16.37, 2.35, -3.70), lat = c(48.21, 48.86, 40.42))
result <- hexify(df, lon = "lon", lat = "lat", area = 5000)

# Quick plot
hexify_plot(result, col = "lightblue", border = "darkblue")
points(df$lon, df$lat, pch = 19, col = "red")
} # }
```
