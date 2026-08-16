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
