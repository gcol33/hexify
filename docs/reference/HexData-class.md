# HexData Class

An S4 class representing hexified data. Wraps the user's data with a
reference to the grid specification used, enabling downstream operations
without repeated parameter specification.

## Usage

``` r
# S4 method for class 'HexData'
grid(x)

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

## Arguments

- x:

  A HexData object

- name:

  Column name for \$ access

- value:

  Value to assign

- i:

  Row indices or logical vector for subsetting

- j:

  Column indices, names, or logical vector for subsetting

- ...:

  Additional arguments passed to underlying methods

- drop:

  Whether to drop dimensions when subsetting (default FALSE)

- object:

  A HexData object (for show/as.data.frame methods)

- row.names:

  Row names for as.data.frame conversion

- optional:

  Logical; if TRUE, row.names may be omitted

## Details

HexData objects are created by
[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) and can
be used with standard R functions. Use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) to
extract the underlying data as a plain data frame.

## Functions

- `grid(HexData)`: Extract grid specification

- `cells(HexData)`: Extract unique cell IDs

- `n_cells(HexData)`: Count unique cells

- `nrow(HexData)`: Get number of rows

- `ncol(HexData)`: Get number of columns

- `dim(HexData)`: Get dimensions

- `names(HexData)`: Get column names

- `$`: Access columns via \$

- `` `$`(HexData) <- value ``: Set columns via \$\<-

- `[`: Subset rows/columns

- `[[`: Subset single column

- `` `[[`(x = HexData, i = ANY, j = missing) <- value ``: Set single
  column

- `show(HexData)`: Print summary

- `as.data.frame(HexData)`: Convert to data.frame

- `as.list(HexData)`: Convert to list

## Slots

- `data`:

  Data frame or sf object. The underlying data with cell assignments.

- `grid`:

  HexGrid object. The grid specification used for hexification.

- `mapping`:

  List. Column name mappings (lon, lat, geometry columns used).

- `kind`:

  Character. Data type: "points", "cells", or "unknown".

- `meta`:

  List. Additional metadata (e.g., cached polygons).

## Compatibility

HexData objects preserve the structure of the underlying data:

- If input was data.frame, output data slot is data.frame

- If input was sf, output data slot is sf (geometry preserved)

- Subsetting operations work transparently via `[` method

## See also

[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) for
creating HexData objects,
[`HexGrid-class`](https://gcol33.github.io/hexify/reference/HexGrid-class.md)
for grid specifications
