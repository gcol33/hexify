# HexGrid Class

An S4 class representing a hexagonal grid specification. Stores all
parameters needed for grid operations.

## Usage

``` r
# S4 method for class 'HexGrid'
x$name

# S4 method for class 'HexGrid'
names(x)

# S4 method for class 'HexGrid'
show(object)

# S4 method for class 'HexGrid'
as.list(x, ...)
```

## Arguments

- x:

  HexGrid object

## Details

Create HexGrid objects using the
[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md)
constructor function. Do not use `new("HexGrid", ...)` directly.

The aperture can be "3", "4", "7" for standard grids, or "4/3" for mixed
aperture grids that start with aperture 4 and switch to aperture 3.

## Functions

- `$`: Get aperture value

- `names(HexGrid)`: Get slot names

- `show(HexGrid)`: Print summary

- `as.list(HexGrid)`: Convert to list

## Slots

- `aperture`:

  Character. Grid aperture: "3", "4", "7", or "4/3" for mixed.

- `resolution`:

  Integer. Grid resolution level (0-30).

- `area_km2`:

  Numeric. Cell area in square kilometers.

- `diagonal_km`:

  Numeric. Cell diagonal (long diagonal) in kilometers.

- `crs`:

  Integer. Coordinate reference system (default 4326 = WGS84).

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md) for
the constructor function,
[`HexData-class`](https://gillescolling.com/hexify/reference/HexData-class.md)
for hexified data objects
