# Get parent index

Returns the parent index (one resolution coarser).

## Usage

``` r
hexify_get_parent(
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

Parent index string
