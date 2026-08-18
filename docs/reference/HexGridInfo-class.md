# HexGridInfo Class

An S4 class representing a hexagonal grid specification. Stores all
parameters needed for grid operations.

## Details

Create HexGridInfo objects using the
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)
constructor function. Do not use `new("HexGridInfo", ...)` directly.

The aperture can be "3", "4", "7" for grids that refine by one aperture
at every level; a family name such as "4/3" or "4/7", which refines by
the first aperture for the first floor(resolution / 2) levels and by the
second for the rest; or one aperture per level, "4/4/7/3".

For H3 grids, the aperture is fixed at "7" and resolution ranges from 0
to 15.

## Slots

- `aperture`:

  Character. Grid aperture: "3", "4", "7", a mixed family such as "4/3"
  or "4/7", or one aperture per resolution level ("4/4/7/3").

- `resolution`:

  Integer. Grid resolution level (0-30 for ISEA, 0-15 for H3).

- `area_km2`:

  Numeric. Cell area in square kilometers.

- `diagonal_km`:

  Numeric. Cell diagonal (long diagonal) in kilometers.

- `crs`:

  Integer or character. Coordinate reference system: an EPSG code, or a
  'PROJ' or 'WKT' string. Defaults to 'WGS84' on Earth, and to a longlat
  CRS on the sphere of `radius_km` on any other body.

- `grid_type`:

  Character. Grid system: "isea" (default) or "h3".

- `radius_km`:

  Numeric. Radius of the body the grid covers, in kilometers. `NA` reads
  as Earth's mean radius.

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
the constructor function,
[`HexData-class`](https://gillescolling.com/hexify/reference/HexData-class.md)
for hexified data objects
