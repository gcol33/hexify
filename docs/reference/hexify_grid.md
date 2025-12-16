# Create a hexagonal grid specification

Creates a discrete global grid system (DGGS) object with hexagonal cells
at a specified resolution. This is the main constructor for hexify
grids.

`hexify_construct()` was renamed to `hexify_grid()` for clarity. Use
`hexify_grid()` instead.

## Usage

``` r
hexify_grid(
  area,
  topology = "HEXAGON",
  metric = TRUE,
  resround = "nearest",
  aperture = 3,
  projection = "ISEA"
)

hexify_construct(
  area,
  topology = "HEXAGON",
  metric = TRUE,
  resround = "nearest",
  aperture = 3,
  projection = "ISEA"
)
```

## Arguments

- area:

  Target cell area in km² (if metric=TRUE) or area code

- topology:

  Grid topology (only "HEXAGON" supported)

- metric:

  Whether area is in metric units (km²)

- resround:

  How to round resolution ("nearest", "up", "down")

- aperture:

  Aperture sequence (3, 4, or 7)

- projection:

  Projection type (only "ISEA" supported currently)

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

  Map projection ("ISEA")

- index_type:

  Index encoding type ("z3", "z7", or "zorder")

## See also

[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) for the
main user function,
[`hexify_grid_to_cell`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md)
for coordinate conversion

Other hexify main:
[`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a grid with ~1000 km² cells
grid <- hexify_grid(area = 1000, aperture = 3)
print(grid)

# Create a finer resolution grid (~100 km² cells)
fine_grid <- hexify_grid(area = 100, aperture = 3, resround = "up")
} # }
```
