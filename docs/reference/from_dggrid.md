# Convert 'dggridR' grid object to hexify_grid

Creates a hexify_grid object from a 'dggridR' dggs object. This allows
using hexify functions with grids created by 'dggridR' dgconstruct().

## Usage

``` r
from_dggrid(dggs, radius_km = EARTH_RADIUS_KM)
```

## Arguments

- dggs:

  A 'dggridR' grid object from dgconstruct()

- radius_km:

  Radius the grid is sized on, in km or as a body name. A dggs carries
  no radius, so the default is Earth's; pass this when the dggs
  describes a grid on another body.

## Value

A hexify_grid object

## Details

The 'ISEA' and 'FULLER' projections with HEXAGON topology are supported.
Other configurations will generate warnings.

A dggs has no field for the body it is sized on, so the resolution is
all that carries over and the result is an Earth grid unless `radius_km`
says otherwise. Cell area follows from the radius, so it is read from
the radius given rather than from the dggs.

The orientation carries over: `pole_lon_deg`, `pole_lat_deg` and
`azimuth_deg` place the icosahedron as DGGRID's `dggs_vert0_lon`,
`dggs_vert0_lat` and `dggs_vert0_azimuth` do, and a field left out takes
its standard ISEA value.

The function validates that the 'dggridR' grid uses compatible settings:

- Projection must be 'ISEA' or 'FULLER'

- Topology must be "HEXAGON" (DIAMOND, TRIANGLE not supported)

- Aperture must be 3, 4, or 7

## See also

Other 'dggridR' compatibility:
[`as_dggrid()`](https://gillescolling.com/hexify/reference/as_dggrid.md),
[`dggrid_43h_sequence()`](https://gillescolling.com/hexify/reference/dggrid_43h_sequence.md),
[`dggrid_is_compatible()`](https://gillescolling.com/hexify/reference/dggrid_is_compatible.md)
