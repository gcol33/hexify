# Convert longitude/latitude to hexagonal cell hierarchical index

Converts geographic coordinates (longitude, latitude) to hexagonal cell
hierarchical index strings. These strings encode the face, resolution,
and cell location in a Z-order (Morton code) format.

## Usage

``` r
hexify_lonlat_to_h_index(grid, lon, lat)
```

## Arguments

- grid:

  Grid specification from hexify_grid()

- lon:

  Longitude vector in degrees (numeric, -180 to 180)

- lat:

  Latitude vector in degrees (numeric, -90 to 90)

## Value

Data frame with columns:

- h_index:

  Hierarchical index (character string)

- face:

  Icosahedron face number (integer, 0-19)

## Details

Most users should use
[`hexify_lonlat_to_cell`](https://gillescolling.com/hexify/reference/hexify_lonlat_to_cell.md)
or
[`hexify_grid_to_cell`](https://gillescolling.com/hexify/reference/hexify_grid_to_cell.md)
which return DGGRID-compatible integer cell IDs.

This function returns hierarchical index strings useful for:

- Understanding the cell's position in the hierarchy

- Prefix-based spatial queries

- Parent/child cell relationships
