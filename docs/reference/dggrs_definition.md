# Definition of a grid as a Discrete Global Grid Reference System

Describes a grid in the DGGRS definition schema of OGC API - Discrete
Global Grid Systems: the hierarchy of grids (solid, refinement ratio,
refinement strategy, zone shapes, orientation), the zone identifiers and
the order of sub-zones.

## Usage

``` r
dggrs_definition(grid)
```

## Arguments

- grid:

  A HexGridInfo or HexData object.

## Value

A named list with the schema's fields `title`, `description`, `dggh`,
`zirs` and `subZoneOrder`.
`jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE)` writes it as
JSON.

## Details

The definition describes the grid as hexify builds it. The zone
identifiers are the strings
[`cell_to_index`](https://gillescolling.com/hexify/reference/cell_to_index.md)
writes, which carry their resolution and are unique across all
resolutions, and the cell IDs, which number the cells of one resolution
from 1 and are unique only together with it.

The ISEA3H and ISEA7H definitions registered with OGC
(<https://www.opengis.net/def/dggrs/OGC/1.0/ISEA7H>) are other grids:
they convert WGS84 geodetic latitude to authalic latitude before
projecting, and place the icosahedron's first vertex at 11.20 degrees E
and authalic latitude \\\arctan(\phi)\\, where hexify, like DGGRID,
reads latitude on the sphere and places it at 11.25 degrees E. A hexify
grid is therefore described without an OGC URI.

## References

Open Geospatial Consortium. OGC API - Discrete Global Grid Systems -
Part 1: Core. OGC 21-038r1.
<https://docs.ogc.org/is/21-038r1/21-038r1.html>

## See also

[`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md),
[`cell_to_index`](https://gillescolling.com/hexify/reference/cell_to_index.md)

## Examples

``` r
def <- dggrs_definition(hex_grid(resolution = 5, aperture = 7))
def$dggh$definition$refinementRatio
def$zirs$textZIRS$description
```
