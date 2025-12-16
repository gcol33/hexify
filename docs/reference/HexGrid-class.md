# HexGrid Class

An S4 class representing a hexagonal grid specification. Stores all
parameters needed for grid operations so they don't need to be repeated
in downstream function calls.

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

- name:

  Slot name for \$ access

- object:

  A HexGrid object (for show method)

- ...:

  Additional arguments (ignored)

## Details

Create HexGrid objects using the
[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md)
constructor function. Do not use `new("HexGrid", ...)` directly.

## Functions

- `$`: Get aperture value

- `names(HexGrid)`: Get slot names

- `show(HexGrid)`: Print summary

- `as.list(HexGrid)`: Convert to list

## Slots

- `aperture`:

  Integer. Grid aperture (3, 4, or 7).

- `resolution`:

  Integer. Grid resolution level (0-30).

- `area_km2`:

  Numeric. Target cell area in square kilometers.

- `grid_system`:

  Character. Grid system identifier (default "ISEA").

- `topology`:

  Character. Grid topology (default "H" for hexagon).

- `index_type`:

  Character. Index encoding type ("integer" or "character").

- `crs_input`:

  Integer. Input coordinate reference system (default 4326).

- `crs_work`:

  Integer. Working CRS for internal operations.

- `mixed_aperture`:

  Logical. Whether mixed aperture (4/3) is used.

- `mixed_aperture_level`:

  Integer. Level at which aperture switches (for 4/3).

- `meta`:

  List. Additional metadata for future extensions.

## See also

[`hex_grid`](https://gcol33.github.io/hexify/reference/hex_grid.md) for
the constructor function,
[`HexData-class`](https://gcol33.github.io/hexify/reference/HexData-class.md)
for hexified data objects
