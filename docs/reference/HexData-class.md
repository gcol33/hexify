# HexData Class

An S4 class representing hexified data. Contains the original user data
plus cell assignments from the hexification process.

## Details

HexData objects are created by
[`hexify`](https://gillescolling.com/hexify/reference/hexify.md). The
original data is preserved in the `data` slot, while cell assignments
are stored separately in `cell_id` and `cell_center`.

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
get a combined data frame with cell columns.

## Slots

- `data`:

  Data frame or sf object. The original user data (untouched).

- `grid`:

  HexGridInfo object. The grid specification used.

- `cell_id`:

  Cell IDs for each row of data. Numeric for ISEA grids, character for
  H3 grids.

- `cell_center`:

  Matrix. Two-column matrix (lon, lat) of cell centers.

## See also

[`hexify`](https://gillescolling.com/hexify/reference/hexify.md) for
creating HexData objects,
[`HexGridInfo-class`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
for grid specifications
