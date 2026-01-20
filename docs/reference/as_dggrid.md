# Convert hexify grid to 'dggridR'-compatible grid object

Creates a 'dggridR'-compatible grid specification from a hexify_grid
object. The resulting object can be used with 'dggridR' functions that
accept a dggs object.

## Usage

``` r
as_dggrid(grid)
```

## Arguments

- grid:

  A hexify_grid object from hexify_grid()

## Value

A list with 'dggridR'-compatible fields:

- pole_lon_deg:

  Longitude of grid pole (default 11.25)

- pole_lat_deg:

  Latitude of grid pole (default 58.28252559)

- azimuth_deg:

  Grid azimuth rotation (default 0)

- aperture:

  Grid aperture (3, 4, or 7)

- res:

  Resolution level

- topology:

  Grid topology ("HEXAGON")

- projection:

  Map projection ('ISEA')

- precision:

  Output decimal precision (default 7)

## See also

Other 'dggridR' compatibility:
[`dggrid_43h_sequence()`](https://gcol33.github.io/hexify/reference/dggrid_43h_sequence.md),
[`dggrid_is_compatible()`](https://gcol33.github.io/hexify/reference/dggrid_is_compatible.md),
[`from_dggrid()`](https://gcol33.github.io/hexify/reference/from_dggrid.md)
