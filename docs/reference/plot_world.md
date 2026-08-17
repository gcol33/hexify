# Quick world map plot

Simple wrapper to plot the built-in world map.

## Usage

``` r
plot_world(fill = "gray90", border = "gray50", ...)
```

## Arguments

- fill:

  Fill color for countries

- border:

  Border color for countries

- ...:

  Additional arguments passed to plot()

## Value

NULL invisibly. Creates a plot as side effect.

## See also

Other visualization:
[`hexify_heatmap()`](https://gillescolling.com/hexify/reference/hexify_heatmap.md)

## Examples

``` r
# Quick world map
plot_world()

# Custom colors
plot_world(fill = "lightblue", border = "darkblue")
```
