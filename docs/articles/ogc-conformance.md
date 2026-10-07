# OGC Topic 21 Conformance

The Open Geospatial Consortium (OGC) describes discrete global grid
systems (DGGS) in its Abstract Specification Topic 21 (OGC 20-040r3,
consistent with ISO 19170-1). Topic 21 defines a DGGS as a hierarchy of
grids that each cover the globe once, zonal identifiers that name the
cells (zones), and functions to assign data to zones, query them and
exchange results. Its 29 requirements fall into three conformance
classes: common spatio-temporal classes (Requirements 1-5), the DGGS
core (6-19) and the Equal-Area Earth Reference System, EAERS (20-29).
This page sets hexify against each of them, gives the measurements a few
of them ask for, and shows the functions that answer Topic 21’s zone
queries.

## Earth model and area error

Topic 21 lets a DGGS sit on a sphere or an ellipsoid; the misfit of a
sphere to the Earth goes into an area error budget, which Requirement 28
caps at 1%. hexify builds an Earth grid on a sphere whose area is that
of the WGS84 ellipsoid, and reads geodetic latitude as latitude on that
sphere. On the sphere, Snyder’s projection makes every hexagon of a
resolution the same area (pentagons 5/6 of it). On the ellipsoid, a
small region at geodetic latitude $`\varphi`$ then has $`MN/R_q^2`$
times its spherical area, with $`M`$ and $`N`$ the meridional and
prime-vertical radii of curvature and $`R_q`$ the authalic radius:

``` r

a <- 6378137; f <- 1 / 298.257223563; e2 <- f * (2 - f); e <- sqrt(e2)
q <- function(phi) (1 - e2) * (sin(phi) / (1 - e2 * sin(phi)^2) -
                                 log((1 - e * sin(phi)) / (1 + e * sin(phi))) / (2 * e))
rq2 <- a^2 * q(pi / 2) / 2
area_scale <- function(lat_deg) {
  s2 <- sin(lat_deg * pi / 180)^2
  a^2 * (1 - e2) / (1 - e2 * s2)^2 / rq2
}
round(100 * (area_scale(c(0, 30, 45, 60, 90)) - 1), 3)
#> [1] -0.447 -0.113  0.223  0.560  0.899
```

A cell’s area on WGS84 therefore lies between 0.45% below and 0.90%
above the sphere’s, by latitude alone: inside the 1% budget. The cells
of every resolution still tile the ellipsoid exactly, since the sphere
and the ellipsoid have the same area.
`paper/bench/bench_ogc_conformance.R` reads the same factor at the area
centroid of every cell of ISEA3H, ISEA4H and ISEA7H grids, with the same
extremes. The DGGRS definitions OGC registers for ISEA3H and ISEA7H
instead convert geodetic latitude to authalic latitude, which makes
cells equal-area on the ellipsoid itself (issue \#104).

Fuller’s projection (`projection = "fuller"`) is not equal-area, so a
Fuller grid falls outside the EAERS class (Requirements 22 and 29), as
does H3.

## Requirements

“Has” means hexify meets the requirement, “partial” that it meets part
of it, “lacks” that it does not. Requirements 2, 3, 5 and 20 ask that a
system’s data model follow Topic 21’s UML diagrams; hexify is an R
library whose objects were not built from them, and they are not
assessed here.

| Req | Asks for | hexify | How |
|----|----|----|----|
| 1 | Temporal geometry and topology | lacks | Zones are spatial only |
| 2 | Zonal geometry and topology classes | not assessed | UML data model |
| 3 | Spatial location classes | not assessed | UML data model |
| 4 | Temporal reference system by period identifiers | lacks | Zones are spatial only |
| 5 | Reference system by zonal identifiers | not assessed | UML data model |
| 6 | Core reference system data model | has | [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md) holds solid, aperture, resolution, orientation, projection, radius and CRS; [`dggrs_definition()`](https://gillescolling.com/hexify/reference/dggrs_definition.md) adds refinement strategy and grid constraints |
| 7 | A CRS, with coordinate epoch | partial | The grid’s `crs`, WGS84 on Earth; no coordinate epoch |
| 8 | A global domain and its dimension | has | The whole sphere of `radius_km`, two-dimensional |
| 9 | Level 0 covers the domain | has | The twelve vertex cells of the icosahedron; [`n_cells()`](https://gillescolling.com/hexify/reference/n_cells.md) counts every level |
| 10 | Every location in exactly one level-0 zone | has | Point assignment is a function; every cell corner is shared by exactly three cells |
| 11 | Simple zone geometry | has | Hexagons and twelve pentagons (squares on the octahedron) |
| 12 | A position inside each zone | has | [`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md), the cell centre |
| 13 | Globally unique zonal identifiers | has | [`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md) strings carry their resolution (see below) |
| 14 | A grid is all zones of one level | has | One [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md) is one level |
| 15 | Grids ordered by refinement | has | Resolutions 0-30 (ISEA), 0-15 (H3) |
| 16 | Quantization functions | partial | [`hexify()`](https://gillescolling.com/hexify/reference/hexify.md), [`lonlat_to_cell()`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md), [`hex_summarize()`](https://gillescolling.com/hexify/reference/hex_summarize.md), [`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md), [`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md), [`hex_compact()`](https://gillescolling.com/hexify/reference/hex_compact.md); no tile roles |
| 17 | Zone queries | partial | [`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md) (one parent or all overlapping), [`get_children()`](https://gillescolling.com/hexify/reference/get_children.md), [`get_siblings()`](https://gillescolling.com/hexify/reference/get_siblings.md), [`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md) (disk or ring), [`hex_distance()`](https://gillescolling.com/hexify/reference/hex_distance.md); no DE-9IM predicates on zone sets |
| 18 | Read and run external queries | lacks | Queries are R function calls |
| 19 | Translate and broadcast results | partial | [`cell_to_sf()`](https://gillescolling.com/hexify/reference/cell_to_sf.md) (any GDAL format), [`as_dggrid()`](https://gillescolling.com/hexify/reference/as_dggrid.md), [`from_dggrid()`](https://gillescolling.com/hexify/reference/from_dggrid.md), [`import_h3()`](https://gillescolling.com/hexify/reference/import_h3.md), [`dggrs_definition()`](https://gillescolling.com/hexify/reference/dggrs_definition.md); no transport |
| 20 | EAERS data model | not assessed | UML data model |
| 21 | Domain is the Earth model’s surface; equal-sized cells | has | Sphere of WGS84 area; error above |
| 22 | Initial equal-area tessellation from a base polyhedron | has (Snyder) | Fuller is not equal-area |
| 23 | Progressively smaller cells | has | Apertures 3, 4, 7 and any sequence of them |
| 24 | Limit on refinement tied to the area budget | partial | Limit set by 64-bit cell IDs (resolution 30, 29, 21), not by the budget |
| 25 | Complete and unique at every level | has | As 9 and 10, at every resolution |
| 26 | EAERS cells are simple polygons | has | As 11 |
| 27 | Position is the area centroid | partial | The cell centre, within a small offset of the centroid (below) |
| 28 | Area error budget of 1% or less | has | -0.45% to +0.90% on WGS84 (above) |
| 29 | Equal area within the budget | has (Snyder) | Exact on the sphere; Fuller and H3 are not equal-area |

## Zonal identifiers

A cell ID numbers the cells of one resolution from 1, as DGGRID’s SEQNUM
does, so it names a cell only together with its grid. The strings
[`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)
writes carry the resolution in their length, so they name one cell among
all resolutions and serve as the zonal identifier of Requirement 13: Z7
(IGEO7) for aperture 7, Z3 for aperture 3, Z-order for aperture 4, and a
two-digit-per-level string for mixed sequences. Each is hierarchy-based,
one of the indexing methods Topic 21 lists (8.2.4.3).

``` r

g7 <- hex_grid(resolution = 5, aperture = 7)
cell <- lonlat_to_cell(16.37, 48.21, g7)
cell
#> integer64
#> [1] 35965
cell_to_index(cell, g7)
#> [1] "0003246"
cell_to_index(get_parent(cell, g7), hex_grid(resolution = 4, aperture = 7))
#> [1] "000324"
```

## Zone queries

Topic 21’s `parent()` returns every coarser zone a zone overlaps, or
with `inheritID = true` only the one its identifier descends from.
hexify’s default parent is the coarser cell holding the cell’s centre,
the `inheritID = true` answer; `overlapping = TRUE` returns all of them.
The hexagons of successive resolutions do not nest: an aperture-3 cell
on a parent’s corner lies a third in each of three parents, an
aperture-4 cell on a parent’s edge half in each of two, and an
aperture-7 cell off the parent’s centre 11/12 in its parent and 1/12 in
one neighbour.

``` r

g3 <- hex_grid(resolution = 6, aperture = 3)
near <- get_neighbors(lonlat_to_cell(16.37, 48.21, g3), g3, include_self = TRUE)[[1]]
get_parent(near, g3)
#> integer64
#> [1] 532 532 532 523 541 747 747
ov <- get_parent(near, g3, overlapping = TRUE)
lengths(ov)
#> [1] 1 3 3 3 3 3 3
ov[lengths(ov) == 3][1]
#> [[1]]
#> integer64
#> [1] 532 514 523
```

Siblings are the other children of a cell’s parent, and a hollow ring
the cells exactly `k` steps away:

``` r

get_siblings(cell, g7)
#> [[1]]
#> integer64
#> [1] 35916 36015 36064 36112 36162 36211
get_neighbors(cell, g7, k = 2, ring = TRUE)
#> [[1]]
#> integer64
#>  [1] 35671 35720 35769 35768 35915 35868 36015 36063 36161 36162 36211 36259
```

Topic 21 gives the refinement ratios of a DGGS as a list that recurs
when it is shorter than the hierarchy.
[`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md)
reads a short aperture list the same way:

``` r

hex_grid(resolution = 5, aperture = c(4, 3))@aperture
#> [1] "4,3,4,3,4"
```

## Centroid

Requirement 27 asks that a zone’s position be its centroid, “computed as
the geodesic centre of surface area”.
[`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
gives the inverse projection of the centre of the planar hexagon, which
is the area centroid where a cell is symmetric on the sphere and lies
close to it elsewhere. `paper/bench/bench_ogc_conformance.R` measures
the offset: the area centroid of a region of the sphere points along
half the integral of $`\mathbf{x} \times d\mathbf{x}`$ around its
boundary.

| Grid | Resolution | Centre spacing (km) | Median offset (km) | Largest offset (km) | Largest / spacing |
|----|----|----|----|----|----|
| ISEA3H | 3 | 1480 | 4.41 | 72.4 | 0.049 |
| ISEA3H | 6 | 286 | 0.29 | 16.4 | 0.057 |
| ISEA3H | 9 | 55.0 | 0.0088 | 3.30 | 0.060 |
| ISEA3H | 12 | 10.6 | 0.00032 | 0.627 | 0.059 |
| ISEA4H | 4 | 482 | 0.93 | 26.8 | 0.056 |
| ISEA4H | 8 | 30.2 | 0.0026 | 1.79 | 0.059 |
| ISEA7H | 3 | 417 | 0.51 | 22.3 | 0.054 |
| ISEA7H | 5 | 59.5 | 0.010 | 3.54 | 0.059 |
| ISEA7H | 7 | 8.51 | 0.00021 | 0.510 | 0.060 |

Grids up to a few thousand cells are measured whole; finer ones on the
cells holding 20,000 uniform points. The median offset falls by about
two orders of magnitude every time the cells shrink tenfold, as expected
where the projection varies smoothly across a cell. The largest offset
stays near 6% of the distance between neighbouring centres at every
resolution. For data that need the area centroid as the zone’s position,
the centre is therefore an approximation, within the stated fraction of
a cell.

## Binning and the hierarchy

Because the hexagons do not nest, a point’s cell at a fine resolution
and its parent do not always agree with the point’s cell at the parent’s
resolution. Griffin (2026) measures this on 10,000 points. On the
lattice the rate one level up follows from the shares above: aperture 3,
$`2/3 \cdot 2/3 = 4/9`$; aperture 4, $`3/4 \cdot 1/2 = 3/8`$; aperture
7, $`6/7 \cdot 1/12 = 1/14`$. Snyder’s projection is equal-area, so the
sphere gives the same rates up to the cells beside the twelve vertices.
`paper/bench/bench_hierarchy_commutation.R` measures it on 10,000
uniform points:

| Grid                   | Lattice, one level | 1 level | 2     | 3     | 4     | 5     | 6     |
|------------------------|--------------------|---------|-------|-------|-------|-------|-------|
| ISEA3H, resolution 12  | 4/9 = 0.444        | 0.452   | 0.473 | 0.448 | 0.604 | 0.455 | 0.651 |
| ISEA4H, resolution 10  | 3/8 = 0.375        | 0.371   | 0.528 | 0.599 | 0.632 | 0.654 | 0.651 |
| ISEA7H, resolution 8   | 1/14 = 0.071       | 0.070   | 0.061 | 0.066 | 0.067 | 0.065 | 0.066 |
| FULLER7H, resolution 8 | 1/14 = 0.071       | 0.077   | 0.064 | 0.068 | 0.066 | 0.067 | 0.065 |
| H3, resolution 8       | 1/14 = 0.071       | 0.074   | 0.066 | 0.068 | 0.066 | 0.066 | 0.068 |

The standard error of each rate is about 0.005 (0.0025 for aperture 7).
One level up, each equal-area grid gives its lattice rate within two
standard errors. Further up, aperture 7 stays near 6.5%, as Griffin
(2026) reports for H3; apertures 3 and 4 rise to about two thirds. One
level up, the parent from
[`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md)
is the cell a point falls in at the coarser resolution for 93% of points
at aperture 7, 55% at aperture 3 and 63% at aperture 4; an aggregation
that must follow the coarser grid exactly bins the points at that
resolution directly.

## Area and shape against other grids

Kmoch et al. (2022) compared the normalized area (a cell’s area over the
mean) and the isoperimetric quotient $`4\pi A/p^2`$ of ten open-source
DGGS, reading each cell in a Lambert azimuthal equal-area plane on WGS84
centred on the cell.
[`cell_metrics()`](https://gillescolling.com/hexify/reference/cell_metrics.md)
returns both, measured on the sphere, as `normalized_area` and `ipq`.
`paper/bench/bench_kmoch_distortion.R` runs their pipeline on hexify’s
cells and sets the result beside their tables:

| Grid | Resolution | Cells | SD of normalized area, this pipeline | Kmoch et al. | Mean ipq, this pipeline | Kmoch et al. | SD of normalized area on the sphere |
|----|----|----|----|----|----|----|----|
| ISEA3H | 6 | 7,292 | 0.0045 | 0.0040 | 0.8964 | 0.8964 | 0.0068 |
| FULLER3H | 6 | 7,292 | 0.0291 | 0.0287 | 0.9015 | 0.9015 | 0.0303 |
| ISEA4H | 7 | 163,842 | 0.0040 |  | 0.8964 |  | 0.0014 |
| ISEA43H | 8 | 207,362 | 0.0040 |  | 0.8964 |  | 0.0013 |
| ISEA7H | 4 | 24,012 | 0.0040 | 0.0040 | 0.8964 | 0.8964 | 0.0037 |
| ISEA7H | 5 | 168,072 | 0.0048 | 0.0043 | 0.8968 | 0.8969 | 0.0014 |
| FULLER7H | 5 | 168,072 | 0.0295 | 0.0290 | 0.9010 | 0.9010 | 0.0292 |
| H3 | 4 | 288,122 | 0.1334 | 0.1306 | 0.9047 | 0.9047 | 0.1334 |

The pipeline columns leave out cells crossing the antimeridian and the
twelve pentagons, whose walls bend where they cross face edges, so that
their corners joined by straight lines miss part of them; Kmoch et
al. left out cells whose geometry their libraries drew invalid, which
removes a few more. The last column is
[`cell_metrics()`](https://gillescolling.com/hexify/reference/cell_metrics.md)
on every cell, pentagons included, with areas on the sphere.

Cell by cell, the corners hexify draws are DGGRID’s for every ISEA3H,
FULLER3H, ISEA7H and FULLER7H cell Kmoch et al. kept, and H3’s for every
H3 cell, except the aperture-7 cells at odd resolutions whose DGGRID
corners are misplaced (134 cells at resolution 3, 975 at resolution 5).
On the rest, area and isoperimetric quotient agree with theirs to 0.03%
and 2e-5.

Two readings follow. On the sphere, hexify’s ISEA hexagons of one
resolution have exactly one area, and the spread in the last column
comes from the pentagons alone; the spread of about 0.004 that both
pipelines find for every ISEA grid on WGS84 is the latitude factor of
the first section, the same for apertures 3, 4, 7 and 4/3. Aperture 3
and 4 hexagons, which Kmoch et al. left out as not congruent, have the
same area and shape spread as aperture 7.

## OGC API - DGGS

OGC API - Discrete Global Grid Systems (OGC 21-038r1) serves data
organised by a Discrete Global Grid Reference System (DGGRS) and
describes each DGGRS in a JSON schema.
[`dggrs_definition()`](https://gillescolling.com/hexify/reference/dggrs_definition.md)
writes a hexify grid in that schema:

``` r

def <- dggrs_definition(g7)
str(def$dggh$definition)
#> List of 8
#>  $ spatialDimensions : int 2
#>  $ temporalDimensions: int 0
#>  $ crs               : chr "EPSG:4326"
#>  $ basePolyhedron    : chr "icosahedron"
#>  $ refinementRatio   : int 7
#>  $ refinementStrategy: chr [1:2] "centredChildCell" "nodeSharingChildCell"
#>  $ constraints       :List of 1
#>   ..$ cellEqualSized: logi TRUE
#>  $ zoneTypes         : chr [1:2] "hexagon" "pentagon"
def$zirs$textZIRS$description
#> [1] "Z7 (IGEO7): two digits naming the resolution-0 cell (00-11) the zone descends from, then one digit 0-6 per resolution naming the child of the parent: 0 the centred child, 1-6 the six around it; the string DGGRID writes. Dropping the last digit gives the parent's identifier."
```

The ISEA3H and ISEA7H that OGC registers are not hexify’s grids of those
names: they convert geodetic latitude to authalic latitude and place the
icosahedron’s first vertex at 11.20 degrees E and authalic latitude
$`\arctan(\phi)`$, where hexify, like DGGRID, reads latitude on the
sphere and places it at 11.25 degrees E. Their zone identifiers
therefore name other zones, and writing them waits on the
authalic-latitude option (issue \#104).

## References

Gibb, R. (ed.) (2021). Topic 21: Discrete Global Grid Systems, Part 1:
Core Reference System and Operations and Equal Area Earth Reference
System. OGC Abstract Specification 20-040r3.
<https://docs.ogc.org/as/20-040r3/20-040r3.html>

Griffin, B. (2026). Hex9: A Quasi-Authalic, Quasi-Continuous Hexagonal
DGGS on the Reference Ellipsoid. arXiv:2608.00022.

Kmoch, A., Vasilyev, I., Virro, H. and Uuemaa, E. (2022). Area and shape
distortions in open-source discrete global grid systems. Big Earth Data
6(3), 256-275. <https://doi.org/10.1080/20964471.2022.2094926>

Open Geospatial Consortium (2025). OGC API - Discrete Global Grid
Systems - Part 1: Core. OGC 21-038r1.
<https://docs.ogc.org/is/21-038r1/21-038r1.html>
