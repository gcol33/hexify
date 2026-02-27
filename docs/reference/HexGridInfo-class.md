# HexGridInfo Class

An S4 class representing a hexagonal grid specification. Stores all
parameters needed for grid operations.

## Details

Create HexGridInfo objects using the
[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md)
constructor function. Do not use `new("HexGridInfo", ...)` directly.

The aperture can be "3", "4", "7" for standard grids, or "4/3" for mixed
aperture grids that start with aperture 4 and switch to aperture 3.

For H3 grids, the aperture is fixed at "7" and resolution ranges from 0
to 15.

## Slots

- `aperture`:

  Character. Grid aperture: "3", "4", "7", or "4/3" for mixed.

- `resolution`:

  Integer. Grid resolution level (0-30 for ISEA, 0-15 for H3).

- `area_km2`:

  Numeric. Cell area in square kilometers.

- `diagonal_km`:

  Numeric. Cell diagonal (long diagonal) in kilometers.

- `crs`:

  Integer. Coordinate reference system (default 4326 = 'WGS84').

- `grid_type`:

  Character. Grid system: "isea" (default) or "h3".

## See also

[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
the constructor function,
[`HexData-class`](https://gcol33.github.io/hexify/reference/HexData-class.md)
for hexified data objects
