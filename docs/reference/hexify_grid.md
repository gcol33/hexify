# Create a Legacy Grid Specification

Builds the list-based grid specification of earlier hexify versions,
kept for code written against them and for 'dggridR' interchange. New
code uses
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md),
which returns a
[`HexGridInfo-class`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
object and covers every grid hexify supports. A hexify_grid prints as
one line.

## Usage

``` r
hexify_grid(
  area,
  topology = "HEXAGON",
  metric = TRUE,
  resround = "nearest",
  aperture = 3,
  projection = "ISEA",
  radius_km = EARTH_RADIUS_KM
)

# S3 method for class 'hexify_grid'
format(x, ...)

# S3 method for class 'hexify_grid'
print(x, ...)
```

## Arguments

- area:

  Target cell area in km^2 (if metric=TRUE) or area code

- topology:

  Grid topology (only "HEXAGON" supported)

- metric:

  Whether area is in metric units (km^2)

- resround:

  How to round resolution ("nearest", "up", "down")

- aperture:

  Aperture sequence (3, 4, or 7)

- projection:

  Face projection, DGGRID's `dggs_proj`: 'ISEA' (default) or 'FULLER'

- radius_km:

  Radius of the body the grid covers, in kilometers, or a body name such
  as "mars" (default Earth). See
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md).

- x:

  A hexify_grid object

- ...:

  Ignored

## Value

A hexify_grid object containing:

- area:

  Target cell area

- resolution:

  Calculated resolution level

- aperture:

  Grid aperture (3, 4, or 7)

- topology:

  Grid topology ("HEXAGON")

- projection:

  Face projection ("ISEA" or "FULLER")

- radius_km:

  Radius of the body, in kilometers

- index_type:

  Index encoding type ("z3", "z7", or "zorder")

## See also

[`hexify`](https://gillescolling.com/hexify/reference/hexify.md) for the
main user function,
[`hexify_grid_to_cell`](https://gillescolling.com/hexify/reference/hexify_grid_to_cell.md)
for coordinate conversion

Other hexify main:
[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md)

## Examples

``` r
# Create a grid with ~1000 km^2 cells
grid <- hexify_grid(area = 1000, aperture = 3)
print(grid)

# Create a finer resolution grid (~100 km^2 cells)
fine_grid <- hexify_grid(area = 100, aperture = 3, resround = "up")
```
