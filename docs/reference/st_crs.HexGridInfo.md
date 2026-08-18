# Coordinate reference system of a grid

Reads the CRS a grid's coordinates are in, as
[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
does for any spatial object. An Earth grid returns 'WGS84'; a grid built
on another body returns a longlat CRS on the sphere of its radius.

## Usage

``` r
# S3 method for class 'HexGridInfo'
st_crs(x, ...)

# S3 method for class 'HexData'
st_crs(x, ...)
```

## Arguments

- x:

  A HexGridInfo or HexData object

- ...:

  Passed on to
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)

## Value

An object of class `crs`, as
[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
returns: a list carrying the reference system in 'PROJ' and 'WKT' form.
It says how to read the coordinates that the grid's cells, centres and
sf exports come back in.

## Examples

``` r
st_crs(hex_grid(resolution = 5))
st_crs(hex_grid(resolution = 5, radius_km = "mars"))
```
