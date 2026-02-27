# Globe center presets

Named list of lon/lat coordinates for common globe views. Used by
[`plot_globe`](https://gcol33.github.io/hexify/reference/plot_globe.md)
when center is specified as a string.

## Usage

``` r
globe_centers
```

## Format

Named list with elements:

- europe:

  c(10, 50) - Western/Central Europe

- north_america:

  c(-100, 45) - USA and Canada

- south_america:

  c(-60, -15) - Full continent

- africa:

  c(20, 5) - Central Africa

- asia:

  c(100, 35) - China, SE Asia, Japan

- oceania:

  c(135, -25) - Australia, NZ, Indonesia

- middle_east:

  c(45, 25) - Arabian Peninsula, Iran, Turkey

- south_asia:

  c(80, 20) - India, Pakistan, Bangladesh

- pacific:

  c(-160, -10) - Polynesia, Pacific islands

- caribbean:

  c(-70, 18) - Caribbean islands

- arctic:

  c(0, 90) - North pole view

- antarctic:

  c(0, -90) - South pole view

## Examples

``` r
globe_centers$europe
globe_centers$oceania
```
