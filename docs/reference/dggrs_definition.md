# Definition of a Grid as a Discrete Global Grid Reference System

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

The ISEA3H, ISEA7H, IVEA3H and IVEA7H definitions registered with OGC
(<https://www.opengis.net/def/dggrs/OGC/1.0/ISEA7H>) convert WGS84
geodetic latitude to authalic latitude before projecting, and place the
icosahedron's first vertex at 11.20 degrees E and authalic latitude
\\\arctan(\phi)\\:
`hex_grid(aperture = 3, ellipsoid = "WGS84", orientation = "ogc")`, its
aperture-7 twin and the same with `projection = "ivea"`. For such a grid
the definition carries OGC's zone identifiers
(`cell_to_index(form = "textZIRS")`) and links to the OGC definition,
which it matches in its hierarchy and zone identifiers; it carries no
OGC URI, since hexify lists sub-zones in ascending order of cell ID
rather than in OGC's scanlines. On a grid with an ellipsoid the
orientation's latitude is geodetic, as OGC states it.

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
