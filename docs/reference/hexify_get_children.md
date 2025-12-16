# Get children indices

Returns all children indices (one resolution finer).

## Usage

``` r
hexify_get_children(
  index,
  aperture = 3L,
  index_type = c("auto", "z3", "z7", "zorder")
)
```

## Arguments

- index:

  Index string

- aperture:

  Aperture (3, 4, or 7)

- index_type:

  Index encoding. Default "auto".

## Value

Character vector of child indices
