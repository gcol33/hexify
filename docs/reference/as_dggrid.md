# Convert hexify Grid to 'dggridR'-compatible Grid Object

Creates a 'dggridR'-compatible grid specification from a hexify_grid
object. The resulting object can be used with 'dggridR' functions that
accept a dggs object.

## Usage

``` r
as_dggrid(grid)

# S3 method for class 'dggs'
format(x, ...)

# S3 method for class 'dggs'
print(x, ...)
```

## Arguments

- grid:

  A grid from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md),
  or a legacy hexify_grid

- x:

  A dggs object

- ...:

  Ignored

## Value

A list with 'dggridR'-compatible fields:

- pole_lon_deg:

  Longitude of icosahedron vertex 0 (standard 11.25)

- pole_lat_deg:

  Latitude of icosahedron vertex 0 (standard 58.282525588538995)

- azimuth_deg:

  Azimuth of vertex 1 seen from vertex 0 (standard 0)

- aperture:

  Grid aperture (3, 4, or 7)

- res:

  Resolution level

- topology:

  Grid topology ("HEXAGON")

- projection:

  Face projection ('ISEA' or 'FULLER')

- precision:

  Output decimal precision (default 7)

It prints as one line naming the grid in DGGRID's terms.

## Details

A dggs carries no body radius, so a grid built on another body cannot be
expressed as one: 'dggridR' reads every dggs on Earth's radius, and its
areas, spacings and CLS come back on Earth. Converting such a grid
warns.

## See also

Other 'dggridR' compatibility:
[`dggrid_43h_sequence()`](https://gillescolling.com/hexify/reference/dggrid_43h_sequence.md),
[`dggrid_is_compatible()`](https://gillescolling.com/hexify/reference/dggrid_is_compatible.md),
[`from_dggrid()`](https://gillescolling.com/hexify/reference/from_dggrid.md)

## Examples

``` r
dggs <- as_dggrid(hex_grid(resolution = 5, aperture = 3))
dggs
dggrid_is_compatible(dggs)
from_dggrid(dggs)
```
