# Convert dggridR grid object to hexify_grid

Creates a hexify_grid object from a dggridR dggs object. This allows
using hexify functions with grids created by dggridR's dgconstruct().

## Usage

``` r
from_dggrid(dggs)
```

## Arguments

- dggs:

  A dggridR grid object from dgconstruct()

## Value

A hexify_grid object

## Details

Only ISEA projection with HEXAGON topology is fully supported. Other
configurations will generate warnings.

The function validates that the dggridR grid uses compatible settings:

- Projection must be "ISEA" (FULLER not supported)

- Topology must be "HEXAGON" (DIAMOND, TRIANGLE not supported)

- Aperture must be 3, 4, or 7

## Examples

``` r
if (FALSE) { # \dontrun{
library(dggridR)
library(hexify)

# Create dggridR grid
dggs <- dgconstruct(area = 1000, metric = TRUE)

# Convert to hexify format
hgrid <- from_dggrid(dggs)

# Use with hexify functions
result <- hexify_grid_to_cell(hgrid, lon = 0, lat = 45)
} # }
```
