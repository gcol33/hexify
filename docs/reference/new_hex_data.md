# Create a HexData Object (Internal)

Internal constructor for HexData objects. Users should use
[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md)
instead.

## Usage

``` r
new_hex_data(data, grid, cell_id, cell_center)
```

## Arguments

- data:

  Data frame or sf object (original user data, untouched)

- grid:

  HexGridInfo object

- cell_id:

  Cell IDs for each row: integer64, or whole numbers below 2^53

- cell_center:

  Matrix with columns lon, lat for cell centers

## Value

A HexData object
