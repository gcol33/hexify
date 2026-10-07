# Lay out the faces of a grid's solid in the plane

Builds a flat map of the sphere from the faces of the solid an ISEA grid
is built on: the face projection of the grid (Snyder's equal-area
projection, or Fuller's) maps the sphere onto the faces, and the layout
cuts the faces into pieces and places each piece in the plane by a
rotation and a translation (the rhombic layout also shears). Draw a grid
on it with
[`plot`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)`(grid, surface = "net", layout = )`,
and place points on it with
[`net_project`](https://gillescolling.com/hexify/reference/net_project.md).

## Usage

``` r
net_layout(
  x,
  layout = c("plane", "gosper", "gosper_flower", "gosper_land", "rhombic", "land"),
  centre = NULL,
  mirror = FALSE,
  land = TRUE
)
```

## Arguments

- x:

  A HexGridInfo object from
  [`hex_grid`](https://gillescolling.com/hexify/reference/hex_grid.md),
  built on the icosahedron or the octahedron.

- layout:

  Name of the layout; see Details.

- centre:

  For the Gosper layouts, the point `c(lon, lat)` whose tile lies in the
  middle, at the origin. The map is turned so the hexagons have
  horizontal top and bottom sides, with north at that point within 30
  degrees of up the page. `NULL` uses the north pole.

- mirror:

  For the Gosper layouts, `TRUE` builds the hexagons on the other four
  faces of the octahedron, the mirror image of the default tiling.

- land:

  For `layout = "land"` and `"gosper_land"`, the land the net keeps
  joined: `TRUE` for the built-in
  [`hexify_world`](https://gillescolling.com/hexify/reference/hexify_world.md),
  or an sf/sfc object of polygons.

## Value

A `hexify_net` object: a list with the layout's `name`, the grid's solid
in `polyhedron`, the `icosa` argument of the grid, and `pieces`, one
list per placed piece with its `face` (from 0), its `region` (three
corners in the face's triangle coordinates), their `labels`, the affine
map `A` and `b` placing a point t at `A t + b`, and the `tile` it
belongs to.

## Details

The layouts are:

- `"plane"`:

  The faces in the PLANE layout of DGGRID, the layout
  [`hexify_cell_to_plane`](https://gillescolling.com/hexify/reference/hexify_cell_to_plane.md)
  reads: five strips of four triangles on the icosahedron, four diamonds
  side by side on the octahedron.

- `"gosper"`:

  Van de Sande's Gosper World on the octahedron: four regular hexagons,
  each one face of the octahedron and a third of each of its three
  neighbours, the three neighbouring faces cut along the lines from
  their centres to their corners. The tile holding `centre` lies in the
  middle and the other three are attached on alternate sides.

- `"gosper_flower"`:

  The tile holding `centre` surrounded by the six tiles that border it,
  so each of the other three tiles appears twice (Van de Sande 2024,
  Fig. 5). The map is continuous across the sides of the middle tile;
  the three seams run out from its corners that are vertices of the
  octahedron.

- `"gosper_land"`:

  The four Gosper tiles joined along the three hexagon sides, out of the
  twelve, that cross the most land of `land`, so the continents stay
  whole where they can, as in Van de Sande's Fig. 1 (2024).

- `"rhombic"`:

  DGGAL's 5 x 6 rhombic space on the icosahedron (Jacovella-St-Louis et
  al. 2025): each face sheared onto half of a unit square, the ten
  diamonds a staircase of squares. Areas keep one constant. DGGAL's y
  axis points down the page, so a point at DGGAL's (x, y) is placed at
  (x, -y).

- `"land"`:

  A net cut along the solid's edges that cross the least land of `land`,
  so the faces stay joined across land where they can, in the spirit of
  Fuller's Dymaxion map, which also splits two faces. With
  `hex_grid(orientation = "dymaxion", projection = "fuller")` it is a
  Dymaxion-style map of whole faces.

On the octahedron with Snyder's projection the Gosper World is
equal-area, where Van de Sande built it from a gnomonic projection onto
a rhombic dodecahedron. His poles lie at midpoints of octahedron edges,
which is `hex_grid(polyhedron = "octahedron", orientation = "gosper")`;
Rus's four-hexagon map puts them at octahedron vertices, which is the
standard orientation of the octahedron.

## References

Van de Sande, A. (2024). Gosper World: A Hexagonal Map Using Gosper
Fractals. *Bridges 2024 Conference Proceedings*, 507-510.

Rus, J. (2017). Flowsnake Earth. *Bridges 2017 Conference Proceedings*,
237-244.

Jacovella-St-Louis, J., Padilla Ruiz, M., Moreira de Sousa, L.,
Peterson, P., Stefanakis, E. (2025). Interoperable global and local
indexing of discrete global grid systems based on the ISEA projection
for efficient storage, processing and transmission. *Abstracts of the
International Cartographic Association* 10, 126.
[doi:10.5194/ica-abs-10-126-2025](https://doi.org/10.5194/ica-abs-10-126-2025)

## See also

[`net_project`](https://gillescolling.com/hexify/reference/net_project.md)
to place points on a layout

## Examples

``` r
grid <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
                 orientation = "gosper")
world <- net_layout(grid, "gosper")
plot(grid, surface = "net", layout = world, seams = TRUE)
plot(grid, surface = "net", layout = "gosper_flower", land = TRUE)
```
