# Summary of a grid or of gridded data

Reports what printing the object reports, as a list, so that a caller
can read the cell count, area, diagonal or column names without reaching
into slots. Printing an object prints its summary, so the two agree.

## Usage

``` r
# S4 method for class 'HexGridInfo'
summary(object, ...)

# S3 method for class 'hexify_grid_summary'
print(x, ...)

# S4 method for class 'HexData'
summary(object, ...)

# S3 method for class 'hexify_data_summary'
print(x, ...)
```

## Arguments

- object:

  A HexGridInfo or HexData object

- ...:

  Ignored

- x:

  A summary, as [`summary()`](https://rdrr.io/r/base/summary.html)
  returns one

## Value

For a HexGridInfo, a list of class `hexify_grid_summary` carrying
`grid_type`, `aperture`, `resolution`, `area_km2`, `diagonal_km`, `crs`,
`radius_km`, `earth`, `orientation` (`c(vert0_lon, vert0_lat, azimuth)`,
empty for H3), `projection` (`"isea"` or `"fuller"`, `NA` for H3),
`polyhedron` (`"icosahedron"` or `"octahedron"`, `NA` for H3) and
`n_cells`. For a HexData, a list of class `hexify_data_summary` carrying
`rows`, `columns`, `column_names`, `n_cells`, `type`, the `grid` summary
and a `preview` of the first rows. The print methods return their input
invisibly.

## Examples

``` r
grid <- hex_grid(area_km2 = 100000)
summary(grid)$n_cells

df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
summary(hexify(df, lon = "lon", lat = "lat", grid = grid))$column_names
```
