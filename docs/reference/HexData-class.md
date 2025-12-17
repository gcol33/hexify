# HexData Class

An S4 class representing hexified data. Contains the original user data
plus cell assignments from the hexification process.

## Usage

``` r
# S4 method for class 'HexData'
grid_info(x)

# S4 method for class 'HexData'
cells(x)

# S4 method for class 'HexData'
n_cells(x)

# S4 method for class 'HexData'
nrow(x)

# S4 method for class 'HexData'
ncol(x)

# S4 method for class 'HexData'
dim(x)

# S4 method for class 'HexData'
names(x)

# S4 method for class 'HexData'
x$name

# S4 method for class 'HexData'
x$name <- value

# S4 method for class 'HexData'
x[i, j, ..., drop = FALSE]

# S4 method for class 'HexData'
x[[i]]

# S4 method for class 'HexData,ANY,missing'
x[[i, j]] <- value

# S4 method for class 'HexData'
show(object)

# S4 method for class 'HexData'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S4 method for class 'HexData'
as.list(x, ...)
```

## Details

HexData objects are created by
[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md). The
original data is preserved in the `data` slot, while cell assignments
are stored separately in `cell_id` and `cell_center`.

Use [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
get a combined data frame with cell columns.

## Functions

- `grid_info(HexData)`: Extract grid specification

- `cells(HexData)`: Extract unique cell IDs

- `n_cells(HexData)`: Count unique cells

- `nrow(HexData)`: Get number of rows

- `ncol(HexData)`: Get number of columns (includes virtual cell columns)

- `dim(HexData)`: Get dimensions

- `names(HexData)`: Get column names (includes virtual cell columns)

- `$`: Access columns via \$ (includes virtual cell columns)

- `` `$`(HexData) <- value ``: Set columns via \$\<-

- `[`: Subset rows/columns

- `[[`: Subset single column (includes virtual cell columns)

- `` `[[`(x = HexData, i = ANY, j = missing) <- value ``: Set single
  column

- `show(HexData)`: Print summary

- `as.data.frame(HexData)`: Convert to data.frame (includes cell
  columns)

- `as.list(HexData)`: Convert to list

## Slots

- `data`:

  Data frame or sf object. The original user data (untouched).

- `grid`:

  HexGridInfo object. The grid specification used.

- `cell_id`:

  Numeric vector. Cell IDs for each row of data.

- `cell_center`:

  Matrix. Two-column matrix (lon, lat) of cell centers.

## See also

[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) for
creating HexData objects,
[`HexGridInfo-class`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
for grid specifications
