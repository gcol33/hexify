# HexGridInfo S4 Methods

S4 methods for HexGridInfo objects. These provide standard R operations
like `$`, [`names()`](https://rdrr.io/r/base/names.html), `show()`, and
[`as.list()`](https://rdrr.io/r/base/list.html).

## Usage

``` r
# S4 method for class 'HexGridInfo'
x$name

# S4 method for class 'HexGridInfo'
names(x)

# S4 method for class 'HexGridInfo'
show(object)

# S4 method for class 'HexGridInfo'
as.list(x, ...)
```

## Arguments

- x:

  HexGridInfo object

- name:

  Slot name

- object:

  HexGridInfo object (for show)

- ...:

  Additional arguments

## Value

- `$`: The value of the requested slot

- `names`: Character vector of slot names

- `show`: The object, invisibly (called for side effect of printing)

- `as.list`: A named list of slot values
