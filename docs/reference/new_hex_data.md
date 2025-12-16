# Create a HexData Object (Internal)

Internal constructor for HexData objects. Users should use
[`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)
instead.

## Usage

``` r
new_hex_data(data, grid, mapping = list(), kind = "unknown", meta = list())
```

## Arguments

- data:

  Data frame or sf object with cell assignments

- grid:

  HexGrid object

- mapping:

  List of column name mappings

- kind:

  Type of data: "points", "cells", or "unknown"

- meta:

  Additional metadata

## Value

A HexData object
