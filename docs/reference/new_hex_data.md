# Create a HexData Object (Internal)

Internal constructor for HexData objects. Users should use
[`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md)
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

  Numeric vector of cell IDs for each row

- cell_center:

  Matrix with columns lon, lat for cell centers

## Value

A HexData object
