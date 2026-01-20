# Convert 'dggridR' grid object to hexify_grid

Creates a hexify_grid object from a 'dggridR' dggs object. This allows
using hexify functions with grids created by 'dggridR' dgconstruct().

## Usage

``` r
from_dggrid(dggs)
```

## Arguments

- dggs:

  A 'dggridR' grid object from dgconstruct()

## Value

A hexify_grid object

## Details

Only 'ISEA' projection with HEXAGON topology is fully supported. Other
configurations will generate warnings.

The function validates that the 'dggridR' grid uses compatible settings:

- Projection must be 'ISEA' (FULLER not supported)

- Topology must be "HEXAGON" (DIAMOND, TRIANGLE not supported)

- Aperture must be 3, 4, or 7

## See also

Other 'dggridR' compatibility:
[`as_dggrid()`](https://gcol33.github.io/hexify/reference/as_dggrid.md),
[`dggrid_43h_sequence()`](https://gcol33.github.io/hexify/reference/dggrid_43h_sequence.md),
[`dggrid_is_compatible()`](https://gcol33.github.io/hexify/reference/dggrid_is_compatible.md)
