# hexify (development version)

## Breaking changes

* `plot_globe()` is removed. `plot(<grid>)` draws a grid in 3D, on the sphere
  (`surface = "sphere"`) or on the flat faces of the solid the grid is
  built on (`surface = "solid"`), with land fill and country outlines
  from `hexify_world` or any sf polygons (`land =`), a `center` to look down
  on, and the solid's face edges. It takes the grid object, so any
  aperture, mixed sequence, resolution, solid or body reaches it; an H3 grid
  is drawn on the sphere. `globe_centers` names the views as before.

* A grid built with one aperture per level stores it comma-separated:
  `hex_grid(resolution = 4, aperture = c(4, 4, 7, 3))@aperture` is
  `"4,4,7,3"`, no longer `"4/4/7/3"`; `aperture = "4,4,7,3"` is accepted and
  "/" is kept for the two-aperture families such as `"4/3"`. A two-level
  sequence `c(4, 7)` was stored as `"4/7"`, the family name, so its parents
  were cells of the family's resolution-1 grid (aperture 7) rather than of its
  own first level (aperture 4). Cell IDs are unchanged. `get_children()` on a
  per-level grid now says that the spelling names no finer level.

* The aperture-7 Z7 index on the icosahedron is now IGEO7's, string for
  string the one DGGRID writes: every cell at resolutions 0 to 7 agrees with
  DGGRID, on ISEA and FULLER and under a given orientation
  (`paper/bench/bench_dggrid_z7.R`). Previously about two cells in three
  differed: a cell whose ancestry leaves its quad was written as quad + 12 *
  seed instead of the base cell it descends from, and the digits were not
  turned into a polar base cell's frame or past the direction each pentagon
  lacks. DGGRID's encoder is bijective: the two cells #53 reported as sharing
  a string get distinct strings from it, so that report came from hexify's
  own earlier port, not from DGGRID. The axis names of digits 3 to 6 were
  also swapped in hexify, which turned digit rotations the wrong way. A
  pentagon has six children, so `hexify_get_children()` returns six strings
  under one, and `hexify_z7_canonical()` sends a string in the deleted
  subsequence to the cell's own index. Cell IDs are unchanged. Octahedral
  grids keep the quad + 6 * seed leading field, since IGEO7 is defined on the
  icosahedron (#86).

* `hexify_cell_id_to_quad_ij()` is removed. `hexify_cell_to_quad_ij()` takes
  the same arguments and returns the same `quad`, `i`, `j` data frame.

* The inverse of Snyder's projection is now in closed form (Recht 2021), so
  `hexify_set_precision()`, `hexify_get_precision()`, `hexify_set_verbose()`
  and `hexify_projection_stats()` are removed, and `hexify_inverse()` no
  longer takes `tol` or `max_iters`. The inverse agrees with the previous
  Newton iteration to within 1e-13 radians on the icosahedron, octahedron and
  tetrahedron, returns a forward-projected point to within a few times 1e-15
  radians, and is about twice as fast (#84).

## New features

* `hex_grid(polyhedron = "octahedron")` builds an ISEA-family grid on the
  octahedron: Snyder's equal-area projection with g = 54.73561032 deg, G = 45
  deg (Snyder 1992), four diamonds, and six square vertex cells, each two
  thirds of a hexagon's area. Apertures 3, 4, 7 and mixed sequences, the
  hierarchy, indices, neighbours, `cell_area()`, `plot()` and `hex_globe()`
  work as on the icosahedron. The solid is a descriptor (vertices, faces,
  Snyder's constants) from which the face pairing, the quad adjacency and
  the vertex cells are derived, so the icosahedron's grids are unchanged
  bit for bit. Fuller's projection stays icosahedron-only (#94).

* `hexify_forward()`, `hexify_inverse()`, `hexify_which_face()`,
  `hexify_face_centers()` and `hexify_lonlat_to_plane()` take
  `polyhedron = "octahedron"` or `"tetrahedron"`. The tetrahedron carries the
  projection only, and `hex_grid(polyhedron = "tetrahedron")` says why: the
  cell numbering pairs faces into diamonds that leave exactly two vertices as
  single-cell poles, a vertex where three faces meet has to be one of them,
  and all four vertices of the tetrahedron are such vertices.

* `cell_to_sf(shape = , depth = )` draws ISEA cells as Gosper islands
  (`"gosper"`) or as the outline of their descendants (`"descendants"`);
  `plot()`, `hexify_heatmap()` and `st_as_sf()` on gridded data take the same
  arguments. A Gosper island replaces each cell edge `depth` times by three
  segments `1/sqrt(7)` as long, turned by `atan(sqrt(3)/5)`, on the face plane
  of the projection, so on an equal-area grid each island has its cell's area
  and the islands of a grid tile the sphere. The descendant outline of an
  aperture-7 cell is the boundary of its descendants `depth` resolutions down,
  which nests across resolutions; the child lattice of ISEA7H turns one way at
  one resolution and back at the next, so this outline stays close to the
  hexagon (#82).

* `get_children()` on ISEA grids looks children up by `match()` rather than by
  name, so its time grows linearly with the number of cells: all 492 cells of
  an aperture-7 resolution-2 grid, three levels down, go from 123 s to 0.5 s.

* `hex_grid(projection = "fuller")` builds an ISEA-family grid on Fuller's
  projection instead of Snyder's equal-area one: DGGRID's FULLER3H,
  FULLER4H, FULLER7H and FULLER43H, and any aperture sequence or orientation.
  The forward projection follows Gray (1995) and the inverse solves Gray's
  equation (39) in Crider's (2008) form in closed form, as a cubic solved
  by Viete's trigonometric method (#96). Cell IDs, the
  hierarchy and neighbours are those of the ISEA grid; centres and corners
  agree with DGGRID's FULLER output. Fuller cells are not equal-area, so
  `area_km2` is the mean cell area and `cell_area()` returns each cell's own,
  from its solid angle. `hexify_forward()`, `hexify_forward_to_face()` and
  `hexify_inverse()` take `projection`, `as_dggrid()` and `from_dggrid()`
  carry FULLER, `dggrid_is_compatible()` accepts it, and `hex_globe()` finds
  Fuller cells per pixel as it does ISEA cells.

* `hex_globe(<grid>)` draws a grid on an interactive globe, rendered on the
  graphics card through WebGPU (an htmlwidget; needs 'htmlwidgets'). Drag turns
  the globe, Shift-drag tilts and turns the view, the wheel zooms, and a
  slider folds the icosahedron into the sphere. `values` fills the cells
  through a colour ramp (`palette`, `limits`, `na_fill`). It takes the camera
  and style arguments of `plot(<grid>)`. The cells of an ISEA grid are found
  per pixel on the graphics card, with `lonlat_to_cell()`'s projection and
  numbering, so the page carries the grid's description and one value per
  cell rather than cell outlines: at aperture 3, resolution 10 (590,492
  cells), a globe with a value on every cell builds in 0.4 s, its page is
  7 MB and a frame takes about 3 ms at 800 and at 1600 pixels. The shader's
  cells agree with `lonlat_to_cell()` to within 1e-5 radians. Borders keep
  their width at every zoom and fade out where cells shrink to a few pixels.
  The pointer shows the ID and value of the cell under it. H3 cells are drawn
  from their outlines. Without WebGPU the widget shows a notice.

* `hex_globe_png()` saves a globe as a PNG, drawn by the widget's own WebGPU
  renderer in headless Chrome (needs 'chromote').

* Land is triangulated on the faces of the icosahedron in C++, so filled land
  folds with the grid; the land area of the mesh matches `sf::st_area()` of
  `hexify_world`.

* Both surfaces read their cell boundaries from the same points: each cell
  edge is walked in the face plane, where it is straight, and cut where it
  crosses a face edge, so a cell on the icosahedron is the cell on the sphere
  folded flat.

* `plot(<grid>)` takes a viewpoint. `projection = "perspective"` places a
  camera `distance` sphere radii from the centre, `tilt` swings it about the
  point at `center`, which it keeps looking at, until the horizon shows, and
  `fov` sets its field of view; by default the frame holds the whole visible
  sphere. `rotation` turns either projection about the line of
  sight. Filled shapes are clipped at a plane in front of the camera, so any
  tilt below 90 degrees draws.

* `hex_grid(orientation = )` places an ISEA grid's icosahedron anywhere on
  the sphere (#80): `c(vert0_lon, vert0_lat, azimuth)` as DGGRID's
  `dggs_vert0_lon`, `dggs_vert0_lat` and `dggs_vert0_azimuth` take it,
  `"random"` (repeatable with `set.seed()`), or `"region"` with `region =` a
  point or an sf object, which centres the grid as DGGRID's `REGION_CENTER`
  does (on the midpoint of an icosahedron edge), or `"face"` with `region =`,
  which puts the region's centre on the centre of a face. The grid stores it in a new `orientation` slot, and every function
  taking the grid reads it, the `hex_globe()` shader included; grids saved
  before the slot existed read the standard orientation. Cell IDs, the
  hierarchy and neighbours are the same under every orientation. Under three
  given orientations and three region placements, 50,000 random points per
  grid fall in the same cells as in DGGRID for apertures 3, 4, 7, ISEA43H and
  a mixed sequence, with cell centres within 2 mm. H3 grids keep H3's own
  orientation.

* `as_dggrid()`, `from_dggrid()` and `dggrid_is_compatible()` carry the
  orientation (`pole_lon_deg`, `pole_lat_deg`, `azimuth_deg`) instead of
  warning about or rejecting a non-standard one.

* `hexify_build_icosa()` sets the orientation of the functions that take no
  grid, and no longer reaches grids; `hex_grid()` and `hexify_grid()` no
  longer reset it.

* `hex_grid(orientation = "dymaxion")` places the icosahedron in Fuller's
  Dymaxion orientation (#89; Sahr et al. 2003): vertex 0 at 5.2454W, 2.3009N,
  azimuth 7.46658, with all twelve pentagons in the ocean. With
  `projection = "fuller"` the grid lies on Fuller's Dymaxion map.

* New `cell_metrics()` and `wall_metrics()` measure a grid's cells: area,
  perimeter and compactness (White et al. 1998) per cell, and per wall
  between adjacent cells its length, the distance between the two centres and
  the cell wall midpoint ratio (Gregory et al. 2008). Without `cell_id` they
  measure every cell of the grid. ISEA walls are followed on the sphere as the
  projection curves them, for Snyder's and Fuller's projections alike; H3
  walls are read from H3's directed edges.

## Bug fixes

* `get_parent()` on a mixed aperture grid finds the parent in exact lattice
  coordinates (#85). After an aperture-3 step a child centre sits on a
  corner of three parents, and after an aperture-4 step on an edge shared by
  two; the parent of such a child was decided by floating-point rounding of
  the projected centre, so a hexagonal parent received anywhere from 1 to 7
  children. A fixed tie rule now gives every hexagonal parent inside a quad
  exactly as many children as the step's aperture. Children whose centre
  lies strictly inside a parent keep their parent. `get_children()` and
  `cell_to_index()` follow `get_parent()`, so index strings of mixed grids
  change for the tied children; regenerate stored ones. `levels > 1` steps
  up one level at a time.

* Cells in the last column of a quad at aperture 3, resolutions 25-30, were
  assigned to another quad: the quad edge dimension was computed in floating
  point and came out one too small. It is now exact.

* `cell_area()` on an ISEA grid returns each cell's area: every hexagon of a
  resolution has `S / (N - 2)` of a body of area `S` split into
  `N` cells, and each of the 12 pentagons 5/6 of that. It returned the
  mean `S / N` for every cell. The `cell_area_km2` column of a HexData
  object and of `as.data.frame()` reads `cell_area()`, so it carries the same
  areas, on the Fuller projection each cell's own.

* `hexify_world` is simplified as a coverage, so neighbouring countries share
  their simplified borders exactly. It was simplified country by country,
  which left about 480 slivers of gap between neighbours; they showed as thin
  lines wherever land was filled without country outlines. Attributes are
  unchanged.

* A cell at an icosahedral vertex now drops the corner that lies in the
  icosahedron's angular deficit, found from the faces rather than a fixed
  corner number. Under the lattice turn of a mixed sequence such as
  `c(4, 4, 7, 3)` the fixed number named the wrong corner, and those twelve
  pentagons were drawn with a corner read off the wrong face.

* On a mixed aperture sequence with an odd number of aperture-7 levels
  (`"4/7"` at odd resolutions, `c(4, 4, 7, 3)`, `c(7, 3)`, ...), the cell
  lattice is turned by atan(sqrt(3)/5) and its cells straddle the quad edges,
  so a point near an edge can have its nearest centre in the next quad.
  `lonlat_to_cell()` clamped that centre back onto the quad and read it as an
  unrelated cell: 2,500 to 4,200 boundary points per grid went to the wrong
  cell, and `grid_global()` missed a cell of `c(4, 4, 7, 3)` and of
  `c(4, 3, 7)`. The centre now moves to the quad that owns it through
  DGGRID's edge table, as pure aperture 7 already did.

* At odd aperture-3 resolutions, a point on or just across an icosahedron
  face edge could quantize to a centre two rows outside its quad, which
  `lonlat_to_cell()` read as an unrelated cell (often a pentagon). Such
  centres now go through the same edge table. Points away from face edges
  keep their cells.

* The Snyder projection constants are computed from their closed forms at
  full double precision (tan θ = 3 - sqrt(5), a face edge of
  sqrt(4π / (5 sqrt(3))), R' = edge / (sqrt(3) tan θ), vertex latitude
  atan(1/2)) instead of being carried to 8 to 10 digits, and the icosahedron's
  vertices are rotated into place on unit vectors, which places them to
  4e-16 where they were off by 1e-8. Cell assignments move only for points
  within about 1e-8 rad (6 cm) of a cell boundary: 2 of 3.8 million random
  points across 19 grids.

* Forward and inverse projection put an azimuth that lies exactly on a
  120° sector boundary, a line from a face centre through a vertex, in
  different sectors, so the inverse returned such a point on the mirrored
  side of the face. Both now read one half-open sector rule. A point on a
  face edge now gets one of the two cells beside it at every resolution.

# hexify 0.8.4

## Breaking changes

* ISEA cell polygons follow the true cell edges. An edge is straight in the
  projection plane and curved in lon/lat, so `cell_to_sf()`, `grid_rect()`,
  `grid_global()`, `hexify_cell_to_sf()` and the plot methods now split each
  edge until every piece is a straight lon/lat chord to within 0.1% of its
  length, as DGGRID's densification does. Drawn with their six (or five)
  corners alone, cells were off by up to 4% of their area at aperture 3
  resolution 8 and up to 19% at coarse resolutions; now within about 0.1%. A polygon now has
  six corners and the points between them, so a cell's vertex count is no
  longer six; the corners themselves are unchanged.

## Bug fixes

* At odd aperture-7 resolutions, `lonlat_to_cell()`, `hexify()` and the index
  functions gave a point near a cell edge to a neighbouring cell: 7 to 11% of
  random points over central Europe at resolutions 2 to 7. A point was snapped to the finer
  substrate lattice and then to that point's parent, so a cell was the union
  of seven substrate hexagons rather than a hexagon. Points now go to the cell
  whose centre is nearest, the hexagon its polygon draws, which changes the
  cell ID of those points.

* A cell straddling a quad edge was drawn with the corners beyond the edge
  read on its centre's face stretched across it, and an icosahedral-vertex
  cell with its folded edge drawn straight across the missing sector. Both
  are now drawn where `lonlat_to_cell()` places their boundary.

* `grid_rect()` on ISEA grids missed cells that intersect the bounding box, mostly
  along its edges and occasionally inside it. It now seeds from a point sampling and
  grows through neighbouring cells that meet the box until none is added, so the
  result is every cell meeting the box, independent of the sampling density (#77).
  `hexify_grid_rect()` had its own copy of the sampling and now calls `grid_rect()`.

# hexify 0.8.3

## Breaking changes

* `as_sf()` is now `sf::st_as_sf()`. sf owns the verb, so hexify registers
  methods on it for `HexData` and `HexGridInfo` and re-exports it, and an sf
  pipeline reaches a hexify object without knowing hexify's own name for the
  conversion. `st_as_sf()` on a grid specification returns its global cell set.
  The `geometry` argument follows `...`, as the generic's signature requires, so
  name it: `st_as_sf(x, geometry = "polygon")`.

* The index entry points that returned one value now return one per index:
  `hexify_index_to_cell()` and `hexify_index_to_lonlat()` a data frame with a
  row per index, `hexify_get_children()` a list with an element per index.
  Reading a single index out of them takes `[[1]]` or `[1, ]`.

## New features

* `summary()` on a grid or on gridded data returns what printing it reports, as
  a list: grid type, aperture, resolution, area, diagonal, CRS, radius and cell
  count for a grid; row and column counts, column names, cell count, storage
  type and the grid's own summary for data. Printing an object prints this
  summary, so what it shows and what it returns are the same thing.

## Bug fixes

* `get_children()` returned cell IDs past the end of the child grid on apertures
  3 and 4, and `NA`s on aperture 7, at the twelve icosahedron vertices, so
  `get_parent()` stopped on its own output there. A vertex sits at the corner of
  several quads and has an index spelling in each, and appending a digit to one
  spelling names cells the other quads hold: some digits named nothing, and the
  children spelled under another quad were never reached. Children are now read
  back from `get_parent()` over candidates drawn from both the index expansion
  and the ring of child cells around the parent, so they invert `get_parent()`
  exactly. Between them the parents of a resolution now claim every cell of the
  child grid once, and an aperture-7 pentagon has six children rather than five
  and two `NA`s.

* `hex_compact()` and `hex_uncompact()` refused every ISEA grid but aperture 7,
  and did their arithmetic on the index strings. Both read cells now, so all
  three pure apertures compact: cells are grouped by `get_parent()`, and a
  group holding as many distinct cells as the parent has children replaces
  them, which reads the six children of a cell at an icosahedron vertex rather
  than assuming seven. Uncompaction expands through `get_children()` for the
  same reason, so the children spelled under a neighbouring quad are reached.
  An index string that names no cell of the grid, which appending a digit to a
  vertex cell writes, is refused by name rather than reported as a cell number
  the caller never wrote. A full grid compacts to the twelve resolution-0 cells
  at every aperture.

* `get_parent()` on an aperture-7 grid at resolution 1 returned cell IDs past
  the end of the resolution-0 grid -- up to 83, where resolution 0 holds 12 --
  and only 42 distinct ones, so `get_children()` at resolution 0 came back in
  threes rather than sixes. A Z7 index spells its leading field as
  `quad + 12 * seed`, the seed naming the point the hierarchy walk arrives at,
  and a two-character index was read as a bare quad number. It is decoded now:
  the seed says which of the base cells meeting at the quad's corner the index
  reaches, through the same adjacency the encoder walks. Every resolution-0 cell
  sits on an icosahedron vertex and now has the six children that make up
  resolution 1's 72, at every aperture the parents of a resolution claim each
  cell of the finer one exactly once, and a full grid compacts to the twelve
  resolution-0 cells at apertures 3, 4 and 7 rather than stopping at 52.

* `hexify_lonlat_to_index()` documented itself as vectorised over `lon` and
  `lat` and took one point, and nothing that consumes an index string accepted
  the vector `cell_to_index()` returns: `hexify_index_to_cell()`,
  `hexify_index_to_lonlat()`, `hexify_get_parent()`, `hexify_get_children()`,
  `hexify_get_resolution()`, `hexify_compare_indices()` and
  `hexify_z7_canonical()` each stopped on Rcpp's own scalar-coercion message,
  which named no hexify function and no argument. Every index entry point takes
  vectors now and loops in C++, so a caller pays one call rather than one per
  element, and a missing coordinate or index carries through as a missing
  result. Arguments read in step whose lengths differ are named, with the
  function they were passed to.

* `hex_distance()` over-reported the hop count between cells of one quad, by up
  to three times, and at aperture 7 also under-reported it by half. It read the
  difference of the two cells straight off the substrate the cell IDs are packed
  on, which is the cell lattice only some of the time: at aperture 3's odd
  resolutions one substrate point in three is a cell and the six neighbours
  stand `sqrt(3)` substrate units apart, so a step was counted as three. The
  handedness of the axial formula was wrong as well. The difference is now
  divided by the lattice generator first, which reads it in cell steps, and the
  six neighbours are the six units there. Every adjacent pair of cells is one
  hop apart at apertures 3, 4 and 7, and the distances agree with breadth-first
  search over `get_neighbors()`.

* `get_neighbors()` and `hex_distance()` read a grid's aperture as a single
  number, which is `NA` for a sequence spelling, so every grid built from one --
  `"4/3"`, `"4/7"` or a per-level vector -- stopped on a coerced-aperture
  message while every other operation on it worked. One walk now serves both:
  it takes what a quad edge measures, which substrate lattice the cells sit on,
  and how a quad coordinate maps to and from the plane, and a pure aperture and
  a sequence differ only in those. A mixed grid has the twelve cells with five
  neighbours and six everywhere else that a pure one has, at every sequence and
  resolution tested, and its distances agree with breadth-first search.

* `grid_clip()`, `hex_zonal()` and `hex_extract()` stopped on a grid built for
  another body, with sf's own `st_crs(x) == st_crs(y) is not TRUE`, which names
  neither the function nor the boundary nor the radius that separates them. A
  boundary and a grid on two bodies are now reconciled in one place: two CRSs on
  the same body reproject as before, and where the bodies differ, longitude and
  latitude name the same region on either sphere, so a lon/lat boundary is read
  in the grid's own CRS and the caller is told which two CRSs met. A projected
  CRS carries lengths that do not carry over, and is refused with both named.
  `hexify()` read sf coordinates the same way, transforming them to WGS84
  whatever body the grid was on, and now reads them on the grid's own.

* The `grid_clip()` and `plot_grid()` examples asked for France and got a third
  of the globe: the world map's France carries the overseas departements, and
  both functions cover the bounding box of what they are given, which for all of
  France runs from -62 to 56 degrees of longitude. That built 36707 cells to
  draw metropolitan France. The examples crop to the metropolitan extent, which
  is what they meant, and together they now take under a second rather than
  fifty.

* `as_dggrid()` on a grid built for another body returned a dggs that read as
  an Earth grid, and said nothing. A dggs has no field for the radius, so the
  conversion warns and names the radius being dropped. `from_dggrid()` is the
  same boundary in reverse and takes a `radius_km` argument, defaulting to
  Earth, so a grid on another body round-trips.

* `hexify_heatmap()` handled an RColorBrewer palette only between three levels
  and the palette's own length. A `breaks` argument with ten or more bins passed
  fewer colours than the scale needed and stopped in ggplot2, naming ggplot2
  rather than the palette; fewer than three bins let RColorBrewer's own warning
  through. The palette is now taken at exactly the number of levels asked for,
  by interpolating over its full range past its length, and it matches
  `RColorBrewer::brewer.pal()` colour for colour within it.

* A `colors` name in neither the viridis nor the RColorBrewer family was passed
  to viridisLite, which drew the default palette and warned about its own
  argument. It is an error now, naming both families.

* `n_cells()` had no `HexGridInfo` method, so asking a grid how many cells it
  holds stopped on dispatch while printing the grid reported exactly that
  number. Both read one function now. On a `HexData` the count is still the
  number of distinct cells its rows fall in.

* `tibble::as_tibble()` on a `HexData` reached tibble's default method and
  stopped there, unable to coerce an S4 object. `as_tibble.HexData()` was
  exported rather than registered, so dispatch never found it. It is registered
  on tibble's generic now, and tibble stays in Suggests.

* `cell_to_sf()` returned an empty polygon for every cell whose boundary lay
  within two degrees of a pole, which from resolution 4 at aperture 7 and
  resolution 5 at aperture 4 is each cell nearest the poles. Corners inside
  that band were snapped onto the pole, so once a cell was narrower than the
  band all of its corners landed on one point. A pole falls inside a cell or on
  a cell edge, never on a corner, and corners now keep the position the inverse
  projection gives them.

* `grid_global()` emitted GDAL errors and returned an invalid geometry for the
  cells that cover a pole. The ring of such a cell winds once around the globe,
  which no lon/lat ring closes, and `sf::st_wrap_dateline()` joined the two
  halves it cut along a line of constant latitude, across the cell's own edges.
  The ring is now cut at the antimeridian and carried up to the pole, the way a
  polar cap is drawn in lon/lat.

* A cell split at the antimeridian keeps its shape. The cut ends were joined
  along a line of constant latitude, which is not where the cell's edge runs and
  moved as much as a tenth of the area of a cell at the seam. The ring now
  carries a corner on its own edge where it meets +/-180, and a cell covers the
  same area whether or not it is split.

* Cell rings cross the antimeridian by carrying their longitudes on
  continuously rather than by shifting the negative ones, which placed a ring
  spanning exactly half the globe on both sides of the map at once. Only cells
  reaching past +/-180 go to `sf::st_wrap_dateline()`, which cuts any segment
  half a turn wide and so also split the edge a cell runs over a pole.

# hexify 0.8.2

## Bug fixes

* `lonlat_to_cell()` on a mixed-aperture grid names a cell for a point that
  sits on a quad corner, which is where the twelve icosahedron vertices put the
  poles. The rotate/requantize chain can leave the quad's coordinate range by a
  tie-breaking unit there, and the cell index reads the pair as an unsigned
  offset from the quad origin, so such a point came back out of range and
  `cell_to_index()` stopped on it.

* The default ISEA orientation carries `atan(phi)` at full double precision,
  with `phi = (1 + sqrt(5)) / 2`, placing vertex 0 at latitude
  58.282525588538995 degrees. The constant was stored rounded at the eighth
  decimal, 1.5e-09 degrees from the analytic value, which is 0.16 mm on Earth
  and reaches the assignment of a point that close to a cell boundary.

* `cell_to_sf()` on a mixed-aperture grid reads its cells with the aperture
  sequence the grid was built from. The polygon path took a single aperture,
  which a mixed sequence does not have, and decoded every mixed grid as
  aperture 3: cell IDs past that grid's range stopped the call and the rest
  named different cells. `is_pentagon()` read the same collapsed aperture and
  is on the sequence path as well.

* Aperture-7 hexagons carry the orientation of their own lattice, unrotated at
  even resolutions and turned 19.1 degrees at odd ones. The corners were placed
  30 degrees off the lattice at every resolution, which aimed them at the
  neighbouring cell centres rather than the points between them. Corner
  positions now come from the grid form, so every aperture reads its rotation
  and scale from one place.

* A cell at an icosahedral vertex drops the corner that falls in the
  icosahedron's angular deficit on a rotated lattice too. The dropped corner
  was fixed to the unrotated layout, so aperture-7 pentagons and aperture-3
  pentagons at odd resolutions kept a corner no neighbour reaches and dropped
  one they share.

* `hexify_cell_to_sf()` returns the same boundary whichever way `return_sf` is
  set. The data frame path was a second implementation, without the pentagon
  handling, the fallback for a corner outside its quad, or the extension of
  polar cells to the pole, so it stopped on every pentagon and drew edge cells
  differently. Both paths read one boundary routine now, and the data frame
  gives a pentagon six rows rather than seven.

* `get_parent()` and `get_children()` honour `levels` on pure ISEA grids, where
  they moved one level whatever was asked for.

* `get_children(get_parent(x))` contains `x`. `get_parent_index()` stripped two
  levels off a Z3 or aperture-3 z-order index whenever its digit count was
  even, and one instead of two off an aperture-7 z-order index;
  `get_children_indices()` rebuilt aperture-3 and aperture-4 children from a
  coordinate box rather than appending the child digit, which is not the
  inverse of stripping one. Over every cell ID the round trip goes from 4.3 per
  cent to 100 per cent for aperture 3 at resolution 2, and from an error to
  98.2 per cent for aperture 7, where the cells below a pentagon remain: a
  pentagon has six aperture-7 children, so a Z7 index does not name seven
  distinct descendants there.

* `quad_ij_to_cell()` re-expresses a coordinate that has stepped outside its
  quad in the quad that owns it, and returns `NA` where the icosahedron folds
  at a vertex and no quad does. Such a pair was packed as it stood, giving cell
  IDs past the end of the grid.

* The compiled entry points all carry the `cpp_` prefix. Five unprefixed
  duplicates at the foot of `rcpp_index.cpp` had no callers, and
  `compileAttributes()` wrote each one into `R/RcppExports.R`, where
  `cell_to_index()` collided with the exported `cell_to_index(cell_id, grid)`,
  and which definition survived came down to the order R sources `R/`.

## Documentation

* The README and the theory vignette say ISEA takes apertures 3, 4 and 7 in any
  sequence, which has held since 0.8.0. They named 4/3 as the one mixed sequence
  ISEA reads.

# hexify 0.8.1

## New features

* `crs` takes any coordinate reference system sf reads: an EPSG code as before,
  or a 'PROJ' or 'WKT' string. A grid left to its default is read on its own
  body, so a grid built with `radius_km` carries a longlat CRS on the sphere of
  that radius and an Earth grid keeps EPSG:4326. EPSG codes name Earth
  reference systems, so a grid on another body had no code to name it and was
  reported in 'WGS84'.

* Every function that hands out coordinates reads the grid's CRS through
  `grid_crs()`, the way kilometres go through `grid_radius_km()`. `cell_to_sf()`,
  `as_sf()`, `hex_summarize(geometry = TRUE)` and `hex_extract()` all return a
  grid's own coordinates rather than assuming Earth.

## Bug fixes

* `plot_globe()` returns in about a second on the grids its own examples use,
  where it took eleven minutes. Cells on the hemisphere edge come out of the
  orthographic transform as rings of two points, which GEOS rejects for the
  whole set, so every call fell through to a per-cell repair loop that rewrote
  the whole table once per cell. Those cells are now dropped before the repair,
  which lets the batch call through, and the loop that remains as a fallback
  assembles its geometries once. The cells the function returns are unchanged.

## Documentation

* `hexify_heatmap()` documents the basemap formats it takes, including the
  rasters it already accepted, and names `sf::st_read()` and `terra::rast()` as
  the way a file becomes one.

* The quickstart says what an sf object is where it first uses one, and the
  bodies section of `vignette("workflows")` covers coordinates on another body.

* Reported by Christian Carey.

# hexify 0.8.0

## Breaking changes

* Aperture-7 cell IDs and Z7 index strings change. The IDs now span `[0, 7^res)`
  within a quad and the index carries the hierarchy seed in its leading field,
  both of which the Bug fixes below explain. Aperture-7 IDs and indices stored
  from an earlier version do not name the same cells and need regenerating from
  the coordinates. Aperture-3 and aperture-4 cell IDs are unchanged.

* `hexify_assign()` assigns about a third of points to a different cell, at every
  aperture and effective resolution, because the quantizer built its cube triple
  from the wrong pair of axes. The cell it used to name sits 0.68 to 0.87 cell
  spacings from the point, further than the circumradius, so the point fell
  outside it. Its `id` is now the Z3 index string and its `face` the quad, and
  the `match_dggrid_parity` argument is gone; it was never wired to an effect.

* `cell_to_lonlat()` returns the centres of the two vertex-quad pentagons at
  (11.25, 58.28) and (-168.75, -58.28), where those cells are under the ISEA
  default orientation, rather than at the geographic poles.

## New features

* Grids cover any body, not just Earth (#56, requested by Christian Carey).
  `hex_grid(radius_km = )` takes a radius in kilometres or a body name --
  `"mars"`, `"moon"`, `"titan"`, `"europa"` and thirteen others, at the IAU mean
  radii (Archinal et al. 2018) tabulated by JPL Solar System Dynamics. The grid
  object carries the radius, so everything reporting kilometres follows it:
  `dgearthstat()`, `cell_area()`, `hexify_compare_resolutions()`, the sf
  exports, and the resolution a target `area_km2` picks. Cell geometry is
  angular and unchanged -- a coordinate lands in the same cell on every body --
  and Earth remains the default, sized against the WGS84 ellipsoid area as
  before.

* Both backends take a radius. H3 reports a cell's area as its solid angle times
  Earth's radius squared, so another radius scales those areas by the square of
  the radius ratio, exactly; `cell_area()`, `dgearthstat()`,
  `hexify_compare_resolutions(type = "h3")`, the resolution an `area_km2` picks
  and `h3_crosswalk()` all follow. An H3 cell ID names a position in H3's
  topology, which Uber's H3 reads on Earth, so the IDs of a grid on another body
  are that topology on that body and are not interchangeable with Earth H3 data;
  `hex_grid()` says so once per session and `h3_crosswalk()` needs both grids on
  the same body.

* `hex_grid()` takes any aperture sequence, so mixed grids beyond ISEA43H are
  reachable from R (#57). A family name splits the levels in two, as `"4/3"`
  already did -- `"4/7"`, `"7/4"`, `"3/7"` -- and a vector names one aperture
  per resolution level, `aperture = c(4, 4, 7, 3)`. Cell IDs, centres,
  neighbours, the geometric hierarchy and the resolution-for-area inversion all
  follow the sequence. The cell count of a sequence is
  `10 * prod(apertures) + 2`.

* Cell IDs pack onto the substrate sublattice of any grid form. A sequence with
  an odd number of aperture-7 levels leaves the cells on a lattice of norm 7 or
  21 rather than the norm-1 or norm-3 lattices pure apertures give; all four are
  now the single congruence `j = c * i (mod N)`, which the aligned and Class II
  packings turn out to be the `N = 1` and `N = 3` cases of.

* Mixed aperture sequences accept aperture 7 alongside 3 and 4, in any order
  (#55, requested by Christian Carey). Scale and lattice orientation now come
  from one model shared by the pure-aperture and mixed-sequence code: refining
  by aperture `a` multiplies the lattice generator by an Eisenstein integer of
  norm `a` (`1 + w` for 3, `2` for 4, `2 + w` for 7), so orientation is the
  product of the steps taken rather than a per-level Class I/II flag. The
  C++ entry points are now `hex_quantize_mixed()`, `hex_center_mixed()` and
  `hex_corners_mixed()`, replacing the `_ap34` names, which no longer describe
  what they accept.

## Bug fixes

* Aperture-7 cells now round-trip through their Z7 index (#53). A cell was
  encoded by walking it up the aperture-7 hierarchy to resolution 0, and decoded
  by walking the digits back down from the origin, but the walk does not arrive
  at the origin for every cell: a quad is a rhombus while the aperture-7 parents
  are hexagons, so the quad boundary cuts through the parents of the cells along
  it, and those arrive at one of the six neighbours of the origin instead. Two
  thirds of all cells do, and the encoder discarded where it landed, so
  `cell_to_index()` followed by `cell_to_quad_ij()` returned the starting cell
  for only 33 to 54 per cent of them. The index now carries the arrival point in
  its leading field as `quad + 12 * seed`, with `seed` one of the seven unit
  digits, and decoding seeds the walk from it. Every cell round-trips at
  resolutions 1 to 9, distinct cells keep distinct indices, and the index length
  is unchanged at `2 + resolution`. A cell whose ancestry stays inside its quad
  has `seed = 0` and keeps the plain two-digit quad, so the two cells DGGRID's
  own encoder merges onto `0045310` now read `2752310` and `2545310`. Read the
  quad from an index with `BB %% 12` or `hexify_index_to_cell()`. A coordinate
  that lies outside the quad it is given is now rejected rather than encoded to
  a string that names a different cell.

* Aperture-7 cell IDs ran outside the documented `[1, 10 * 7^res + 2]` range and
  `validate_cell_id()` rejected them: a uniform sample at resolution 2 reached
  2647 against a nominal count of 492, and at resolutions 1 and 2 the IDs also
  overran the padded space the encoder itself allotted, so `cell_to_quad_ij()`
  errored on ids the assignment had just produced. Aperture 7 stored each cell
  in a square bounding box sized by a fitted constant, roughly three times the
  cell count, because the surrogate lattice sits at 19.1 degrees to the quad
  frame. Cells now index by walking the quad's own substrate box, where the
  centres are the sublattice `2u + v = 0 (mod 7)` and a row holds one seventh of
  the box, so the IDs of a quad span `[0, 7^res)` with no gaps and the cell count
  is `10 * 7^res + 2` for every aperture.

* Aperture 7 gave the same cell two addresses at odd resolutions, once from each
  quad an edge cell straddles: a uniform sample at resolution 1 produced 122
  distinct `(quad, i, j)` for a grid of 72 cells. One cell there covers seven
  substrate points, and the quad was decided from the sampled point rather than
  the cell, so a cell whose points fall on both sides of an edge belonged to
  both. The cell centre now picks the quad, which puts each cell in exactly one:
  every ID in `[1, 10 * 7^res + 2]` names one distinct cell and round-trips
  through `cell_to_quad_ij()`.

* ISEA neighbour offsets stepped along the wrong lattice directions. The aligned
  lattice writes `(i, j)` as `x = i - j/2`, so `(-1, 1)` and `(1, -1)` stand
  `sqrt(3)` units away rather than one, and aperture 3's odd resolutions carry
  cells on the 30 degree lattice, where neighbours are `sqrt(3)` substrate units
  apart. Aperture 3 at resolution 5 returned as few as two neighbours for a cell
  and one ID past the end of the grid; apertures 3 and 4 reported second-ring
  cells, at 1.6 to 2.4 times the cell spacing, as neighbours of the 12 pentagons.
  A neighbour that under-runs its quad's frame has no image in it and used to be
  dropped; it now steps to the owning quad through DGGRID's edge table. Both
  poles are read from an adjacent quad's corner, where the offsets reach the
  five surrounding cells, instead of from their own single-cell frame, which the
  south pole overran into IDs past the end of the grid and the north pole left
  empty. Resolution 0, where all 12 cells are pentagons and each quad holds one
  cell, reads the icosahedron's vertex graph. Every aperture at resolutions 0-4
  now reports exactly 12 pentagons, hexagons elsewhere, a symmetric adjacency,
  and no neighbour centre further than 1.12 cell spacings.

* `cell_to_lonlat()` reported the two vertex-quad pentagons (quads 0 and 11) at
  the geographic poles (#58). Those cells sit at icosahedron vertex 0 and its
  antipode, which are the poles only under a pole-aligned orientation; under the
  ISEA default they are at (11.25, 58.28) and (-168.75, -58.28), about 3500 km
  away. Both now come from folding the cell through the quad frame, as the
  polygon and mixed-sequence paths already did, so a point assigned to one gets
  its own cell's centre back. `index_to_cell()`'s inverse carried the same
  hardcoded pair and is fixed with it.

* `hexify_assign()` returned centres thousands of kilometres outside the
  assigned cell for a share of points -- 113 of 400 uniformly sampled points
  further than two centre spacings at effective resolution 2 (#58). It was the
  one path that quantized against raw triangle-local coordinates instead of
  folding into the non-negative quad frame, so it needed its own Z3 digit scheme
  (magnitude digits plus two trailing sign digits) whose round trip through
  `z3::encode()`/`decode()` did not preserve `(i, j)`. It now runs on the same
  pipeline as `hexify()`: `lonlat_to_cell()` for the cell, `cell_to_lonlat()`
  for the centre, `cell_to_index()` for the ID. Every sampled point is now
  within one circumradius of its cell centre at effective resolutions 1-8. The
  reported `id` is the Z3 index string, `face` is the quad (0-11), and the
  `match_dggrid_parity` argument is gone -- it was never wired to an effect, and
  the shared pipeline is the DGGRID-verified one. `cpp_hex_index_z3_*()` are
  removed with the scheme they served.

* Quantization built the cube triple `(i, j, -i-j)` from the grid coordinates.
  The `i` and `j` axes point at 0 and 120 degrees, so they are two of the three
  cube axes rather than the adjacent pair round-and-fix expects, and the step
  that corrects the coordinate with the largest rounding error corrected the
  wrong one: against brute-force nearest-centre over 3000 points, the quantizer
  returned a cell other than the nearest for 1026 of them. The triple is
  `(i - j, j, -i)`, which now matches brute-force nearest-centre on every
  sampled point. This reached `hex_quantize_ap3/4/7()`, the mixed-sequence
  quantizer and `hexify_assign()`, whose cell assignment changes for about a
  third of points; the quantization is scale-free, so the share is the same at
  every effective resolution. Cell IDs and hierarchical indices are
  unaffected -- they quantize through `quad_xy_to_ij()`, which uses the exact
  region-classifying quantizer.

* The aperture-7 grids returned by the `hex_*_ap7()` helpers were not
  centre-nested: each level applied a fixed `kAp7RotDeg` offset plus a 30-degree
  alternation, which is not an aperture-7 refinement (the norm-7 Eisenstein
  generators sit at +/- `kAp7RotDeg`, not 30 degrees), so a cell centre was
  about 0.85 of a child cell spacing away from the nearest centre one level
  down. Levels now alternate the two norm-7 generators, putting even
  resolutions at 0 degrees and odd resolutions at `kAp7RotDeg` on a `sqrt(7)`
  substrate -- the same convention as the exact-integer aperture-7 route used
  for cell IDs, whose divisor is `7^ceil(res/2)`.

* Mixed sequences computed orientation from the count of aperture-3 levels
  including the base level, and treated any sequence ending in aperture 4 as
  unrotated. A sequence such as `c(4, 4, 3)` was therefore 30 degrees away from
  the grid it describes, and `c(4, 3, 4)` was not nested in `c(4, 3)`.

* The aperture-7 inverse in `quad_ij_to_xy()` scaled by a `sqrt(7)`/`sqrt(21)`
  substrate while the forward quantization used the exact `7^ceil(res/2)` one.
  Callers that reached it -- `hexify_index_to_lonlat()`,
  `hexify_quad_ij_to_xy()` and `hexify_quad_ij_to_icosa_tri()` -- placed
  aperture-7 cells in the wrong location, or outside the quad entirely, where
  the conversion then failed. DGGRID's `DgHexGrid2DS` toggles Class III on
  every aperture-7 level, so even resolutions are an unrotated Class I grid and
  odd ones carry one aperture-7 level; the inverse now takes the same
  exact-integer route as the forward. The float-rotation Class III helpers that
  encoded the old substrate (`quantize_class3i()`, `quantize_class3ii()`,
  `substrate_to_surrogate_ap7()`, `surrogate_to_substrate_ap7()`) were unused
  and have been removed.

* `hexify_index_to_lonlat()` returned the icosahedron vertex instead of the
  pole for cells in the polar pentagon quads, at every aperture. It now answers
  those directly, as `cell_to_lonlat()` does.

# hexify 0.7.5

## Bug fixes

* `kAp7RotDeg` was mislabeled as `atan(sqrt(3/7))` (~33.2 deg) and its numeric
  literal drifted from the true value starting at the 12th significant digit,
  present since the initial release. The correct closed form is
  `atan(sqrt(3)/5)`, now cross-checked against DGGRID's `M_AP7_ROT_DEGS` and
  updated to full double precision. A duplicate literal in
  `coordinate_transforms.cpp` was removed in favor of deriving from the same
  constant (reported by Christian Carey).

# hexify 0.7.4

## Documentation

* Documented the bijective aperture-7 Z7 format, exact round-trip guarantee,
  index structure, and the limited pentagon-region difference from DGGRID's
  non-injective raw Z7 encoding.

## New features

* Hierarchical navigation now works for mixed aperture `"4/3"` (ISEA43H) grids:
  `get_parent()`, `get_children()`, and `cell_to_index()` no longer error on
  these grids (#31). Mixed 4/3 grids have no DGGRID-standard hierarchical index,
  so hexify defines the hierarchy geometrically -- a cell's parent is the
  coarser cell whose lattice point contains the cell's centre. This is the
  correct relationship for how these grids quantise (a single scaled
  quantisation per resolution, not a step-composed subdivision), and it
  round-trips `cell -> cell_to_index() -> cell` exactly, including the seam and
  icosahedron-vertex cells that a purely local walk would miss.

# hexify 0.7.3

**Bug fixes, package hygiene, and test coverage**

## Fixes

* Fixed `hexify_lonlat_to_index()`/`hexify_index_to_lonlat()` for aperture 3:
  quantization skipped the quad-frame fold that aperture 4/7 already had,
  producing indices inconsistent with `lonlat_to_cell()`/`cell_to_index()`
  for essentially all points.
* Fixed `is_pentagon()` for ISEA aperture 7, which undercounted or threw at
  resolution >= 1; it now decodes each cell's own (i, j) instead of relying
  on an unreliable forward computation.
* Fixed `hex_compact()` emitting duplicate cell IDs when the input already
  contained a parent cell alongside all 7 of its children.
* Fixed `hexify_assign()`'s Z3 backend discarding the sign of quantized
  (i, j), causing sign-flipped points to collide on the same cell ID.
* Fixed NaN/Inf inputs reaching undefined behavior in `snyder_forward()`'s
  sort comparator and silently returning garbage from hex quantization
  instead of erroring.
* `cell_to_index()`/`get_parent()`/`get_children()` no longer produce a
  confusing internal bound error (or, with a naive fix, silently wrong cell
  IDs) for mixed aperture `"4/3"` grids. In 0.7.3 these raised a clear
  "not implemented" error; 0.7.4 implements the navigation geometrically
  (see above).
* `hex_extract()` now respects `cells=`/`boundary=` when `grid` is a
  HexData object (previously silently ignored).
* `hex_browse()`'s data.frame input mode no longer crashes on duplicate
  `cell_id` rows.
* `plot_globe()`'s `resolve_center()` now validates a named `center`
  vector's names instead of silently building an invalid PROJ string.
* Hierarchical-index functions (`hexify_cell_to_index()` and siblings) now
  validate `resolution`/`aperture` like `hexify_lonlat_to_cell()` already
  did; the C++ layer now enforces the max resolution for aperture 3/4 the
  same way it already did for aperture 7.
* `as_dggrid()` now accepts a modern `HexGridInfo` object, matching
  `dgverify()`.
* `HexData`'s validity check no longer has a blind spot for empty
  `cell_id`/`cell_center` paired with non-empty `data`.
* `hex_grid()` and legacy `hexify_grid()` now give clear errors for
  non-numeric/NA/non-positive `resolution`, `crs`, or `area`, instead of
  base-R errors or a silent `resolution = NaN`.
* `plot_globe(exclude_antarctica = TRUE)` now warns instead of silently
  no-op'ing for custom `land_data` without a recognized country-name
  column.
* `import_h3()` now drops NA cell IDs (with a warning) when `data` is
  attached, matching `hexify()`'s existing NA-coordinate handling.
* `prepare_fill_column()` now warns when `breaks=` is supplied for a
  discrete value column instead of silently discarding it.

## Package hygiene

* Removed dead code (`index_to_cell_internal()`), an orphaned test
  fixture, and a committed rendered vignette; fixed a native-pipe usage
  that required a newer R than the package declares; fixed stale roxygen
  text; documented `HexData[`'s `drop = FALSE` default explicitly.
* C++: replaced a raw `malloc`/`free` digit buffer with `std::vector`;
  an out-of-range digit now throws instead of silently clamping; added a
  resolution/string-length consistency check to `z3::decode()`;
  `get_children_indices()` no longer swallows unrelated exceptions behind
  a catch-all around the max-resolution boundary check.

## Tests

* Strengthened ~20 `hexify_heatmap()` tests that previously only checked
  the return type is a ggplot object to assert on the actual plot
  data/mapping (fill column, scale type, colors, bins, limits, labels).
* Added coverage for the `Raster*`/`SpatRaster` basemap path and for
  `hex_browse()`'s value-to-fill-color mapping.

# hexify 0.7.2

* Fixed aperture 4/7 lon/lat-to-index and index-to-lon/lat conversion:
  quantization was operating on raw icosahedron triangle coordinates
  instead of folding into the quad frame first, producing wrong indices
  for any point on triangles numbered 12-19.
* Implemented `cpp_hex_index_z3_quantize_digits()`, `cpp_hex_index_z3_center()`,
  and the Z3 corners helper, which previously returned placeholder zeros
  instead of real quantized digits/coordinates.
* Fixed `hex_zonal()` row misalignment: results are now keyed off the
  deduplicated `hex_sf$cell_id` order instead of the pre-dedup input,
  and `cells` input now drops `NA`/duplicate values before lookup.
* Fixed `grid_clip()` argument order in `hex_zonal()`'s boundary path.

# hexify 0.7.1

**Code quality and documentation**

## Improvements

* Refactored `rcpp_aperture.cpp`: extracted parameterized helpers for
  quantize/center/corners, eliminating copy-paste across 3 apertures.
* Consolidated input validation: conversion functions now use shared
  `validate_aperture()` / `validate_resolution()` from `constants.R`.
* Optimized `is_pentagon()` for ISEA grids: computes pentagon cell IDs
  directly from quad coordinates instead of looping through the lon/lat
  pipeline.
* Stricter aperture parsing in `hexify()`: rejects unexpected values with
  an informative error instead of silently coercing via `as.integer()`.
* Deduplicated `theme_minimal()` calls in plot methods with `.theme_clean()`
  helper.

## Documentation

* Added "Why Hexagonal Grids?" section to README (equal area, uniform
  adjacency, low shape distortion).
* Added "Known Limitations" section to README (H3 resolution cap, pentagons,
  projection precision).

## Tests

* New `test-edge-cases.R` with 32 tests covering poles, antimeridian,
  dateline, equator, roundtrip stability, pentagon invariants, neighbor
  symmetry, and resolution-0 cell counts.

# hexify 0.7.0

**Spatial analysis primitives**

## Hotfix

* Fixed aperture 7 cell encoding to use Class III substrate quantization.
* Updated `max_cell_id()` for aperture 7 bounding box.

## New features

* New `get_neighbors()`: returns k-ring (disk) of neighboring cells for
  both ISEA and H3 grids. Supports `k > 1` for multi-ring expansion,
  `distances = TRUE` for ring distance output, and vectorized cell input.
  ISEA backend uses axial coordinate offsets with lon/lat fallback for
  cross-quad boundaries. H3 backend uses vendored `gridDisk` /
  `gridDiskDistances` / `gridRingUnsafe`.
* New `hex_summarize()`: cell-level data aggregation with tidyeval support.
  Groups by cell_id, applies user-defined summary expressions, returns a
  data.frame with cell centers, areas, and point counts. Supports
  `geometry = TRUE` for sf output.
* New C++ bindings: `cpp_h3_gridDisk`, `cpp_h3_gridDiskDistances`,
  `cpp_h3_gridRingUnsafe`, `cpp_get_neighbors_isea`, `cpp_get_neighbors_z7`

# hexify 0.6.5

* Fixed empty translation unit warning in vendored H3 `h3Assert.c`
  (clang 21 `-Wempty-translation-unit`)

# hexify 0.6.4

* Removed compiled object files (`.o`) from source tarball that caused
  installation failure on Linux (Debian) and NOTE on all platforms
* Wrapped `plot_globe()` examples in `\donttest{}` to reduce check time
  (was 608s on win-builder)
* Reworded DESCRIPTION to avoid "vendored" spelling flag

# hexify 0.6.3

* `cell_to_sf()` now applies `sf::st_wrap_dateline()` automatically, fixing
  horizontal streaks on flat map projections (Plate Carrée, Robinson, etc.) for
  hexagons crossing the ±180° antimeridian
* `as_sf(x, geometry = "polygon")` now routes all grids (ISEA and H3) through
  `cell_to_sf()`, ensuring consistent antimeridian handling
* `hexify_cell_to_sf()` gains antimeridian normalization matching `cell_to_sf()`
* Vignettes no longer require manual `st_wrap_dateline()` calls

# hexify 0.6.2

**Native H3 backend — zero external dependencies**

* Vendored H3 v4.4.1 C source, replacing the `h3o` R package dependency
* H3 is now always available — no optional install, no Suggests
* Uses H3 experimental polygon fill for full spatial coverage in `hexify()`
* New native C++ bindings: `h3_lat_lng_to_cell`, `h3_cell_to_boundary`,
  `h3_cell_to_parent`, `h3_cell_to_children`, `h3_polygon_to_cells`,
  `h3_cell_area_km2`, `h3_cell_to_lat_lng`

# hexify 0.6.1

* Warn when `aperture` is passed with `type = "h3"` (ignored parameter)
* Expanded H3 resolution guidance in `hex_grid()` documentation
* Extended H3 vignette resolution table to full range (0-15)
* Added `h3_crosswalk()` example to H3 vignette

# hexify 0.6.0

**ISEA–H3 crosswalk and per-cell area**

* New `h3_crosswalk()`: bidirectional mapping between ISEA and H3 cell IDs,
  with automatic resolution matching and per-cell area comparison
* New `cell_area()`: returns geodesic area (km²) for each cell — constant for
  ISEA (equal-area), location-dependent for H3, with session-scoped caching
* HexData `$cell_area_km2`, `[["cell_area_km2"]]`, and `as.data.frame()` now
  return per-cell areas for H3 grids instead of the grid-wide average
* Internal: extracted `closest_h3_resolution()` helper shared by `hex_grid()`
  and `h3_crosswalk()`

# hexify 0.5.0

**H3 grid support**

* Added H3 (Uber) as a first-class grid type: `hex_grid(resolution = 8, type = "h3")`
* All core functions work with H3 grids: `hexify()`, `cell_to_sf()`, `grid_rect()`,
  `grid_global()`, `grid_clip()`, `get_parent()`, `get_children()`
* H3 support requires the `h3o` package (Suggests, not required for ISEA workflows)
* New `hexify_compare_resolutions(type = "h3")` for H3 resolution table
* `dgearthstat()` now accepts HexGridInfo objects directly
* New `grid_type` slot on HexGridInfo: `"isea"` (default) or `"h3"`
* HexData `cell_id` slot supports character (H3) and numeric (ISEA) cell IDs
* Backward compatible: all existing ISEA workflows unchanged

# hexify 0.3.10

**Hotfix for geometry issues**

* Fixed invalid pentagon geometries that caused gaps in global grids
* Fixed antimeridian-crossing polygons using `st_wrap_dateline()`
* Added polar cap sampling to `grid_global()` to include cells above ±85° latitude

# hexify 0.3.6

* Reduced test suite runtime for CRAN by skipping detailed consistency tests
  (full tests still run locally via NOT_CRAN=true)
* Fixed CRAN incoming check NOTE: "Overall checktime 15 min > 10 min"

# hexify 0.3.5

* Simplified plot examples to reduce runtime
* Added non-standard files to .Rbuildignore
* Fixed slow example NOTE

# hexify 0.3.4

**Hotfix for CRAN UBSAN check failure**

* Fixed undefined behavior in Snyder ISEA projection causing NaN values
  (UBSAN error on M1 Mac CRAN check: "nan is outside the range of representable
  values of type 'long long'" at coordinate_transforms.cpp:243-244)
* Root cause: floating-point precision in face assignment could project points
  onto geometrically invalid triangle faces near icosahedron edges
* Solution: added projection validation with face-retry logic - if validation
  fails (z > DH), automatically tries adjacent faces until valid

# hexify 0.3.3

* CRAN submission fixes

# hexify 0.3.2

* Documentation improvements
* Fixed hierarchical index functions
* Increased test coverage to 90%

# hexify 0.3.1

* Minor bug fixes

# hexify 0.3.0

* Major refactoring of coordinate transformation system
* Improved performance for large grids

# hexify 0.2.0

* Added dggridR compatibility functions (`as_dggrid()`, `from_dggrid()`)
* New hierarchical indexing support with H-index functions
* Improved coordinate conversion pipeline
* Added grid statistics functions

# hexify 0.1.0

* Initial release
* ISEA discrete global grid implementation
* Support for apertures 3, 4, 7, and mixed 4/3
* Compatible with dggridR output
