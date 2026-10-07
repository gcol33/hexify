# Get Number of Cells

Counts cells: those a dataset occupies, or those a grid contains.

## Usage

``` r
n_cells(x)
```

## Arguments

- x:

  A HexData or HexGridInfo object

## Value

For a HexData, the number of distinct cells its rows fall in. For a
HexGridInfo, the number of cells the grid divides the body into, which
runs past integer range at fine resolutions and so comes back as a
double.

## Examples

``` r
grid <- hex_grid(area_km2 = 100000)
n_cells(grid)

df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
n_cells(hexify(df, lon = "lon", lat = "lat", grid = grid))
```
