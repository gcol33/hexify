# Decode a cell index to face, i, j, and resolution

Internal function to decode a cell index string into its constituent
components: face number, grid coordinates (i, j), and resolution level.

## Usage

``` r
index_to_cell_internal(index, aperture, index_type)
```

## Arguments

- index:

  Cell index string

- aperture:

  Grid aperture (3, 4, or 7)

- index_type:

  Index encoding type ("z3", "z7", or "zorder")

## Value

List with components:

- face:

  Face number (integer)

- i:

  Grid coordinate i (integer)

- j:

  Grid coordinate j (integer)

- resolution:

  Resolution level (integer)

## Examples

``` r
if (FALSE) { # \dontrun{
# This is an internal function, typically called by other functions
cell_info <- index_to_cell("0112345", aperture = 3, index_type = "z3")
print(cell_info)
} # }
```
