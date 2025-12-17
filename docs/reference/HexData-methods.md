# HexData S4 Methods

S4 methods for HexData objects. These provide standard R operations for
accessing data, subsetting, and conversion.

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

## Arguments

- x:

  HexData object

- name:

  Column name

- value:

  Replacement value

- i, j:

  Row/column indices

- ...:

  Additional arguments

- drop:

  Logical, whether to drop dimensions

- object:

  HexData object (for show)

- row.names:

  Optional row names

- optional:

  Logical (ignored)
