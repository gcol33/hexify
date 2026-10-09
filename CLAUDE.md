# hexify Development Notes

## Grid Coverage

hexify supports **all** major hexagonal DGGS through two backends:

- **ISEA** (built-in C++): apertures 3, 4, 7, and any mixed sequence of them — resolutions 0-30
  where the cell count fits a signed 64-bit integer (icosahedron: ap3 30, ap4 29, ap7 21),
  on Snyder's equal-area projection or Fuller's (`hex_grid(projection = "fuller")`),
  on the icosahedron or, with Snyder's projection, the octahedron
  (`hex_grid(polyhedron = "octahedron")`).
  `hex_grid(aperture = "4/3")` / `"4/7"` / `"7/4"` name a family (first `floor(res/2)` levels
  take the first aperture), and `aperture = c(4, 4, 7, 3)` names one aperture per level,
  stored as `"4,4,7,3"`: "," marks a per-level spelling and "/" a family, so `"4,7"` and
  `"4/7"` stay apart (same grid at resolution 2, different parents at resolution 1).
- **H3** (vendored H3 v4.4.1 C source in `src/h3`): fixed aperture 7 — resolutions 0-15

This covers every hexagonal grid system that matters:
- ISEA3H, ISEA4H, ISEA7H, ISEA43H = ISEA backend with different aperture settings
- FULLER3H, FULLER4H, FULLER7H, FULLER43H = the same with `projection = "fuller"`
- IVEA3H, IVEA7H (DGGAL, van Leeuwen & Strebe 2006) = the same with `projection = "ivea"`;
  DGGAL's grids sit at `orientation = c(11.20, 58.282525588538995, 0)` on the authalic sphere
- H3 = H3 backend
- IGEO7/Z7 (Kmoch, Sahr, Chan, Uuemaa 2025) = equal-area aperture-7 hex grid with Z7 indexing = already covered by `hex_grid(aperture = 7)` (ISEA7H with Z7 index). hexify uses the same Z7 hierarchical indexing (7-digit encoding, 0-6 per level) in `src/index_z7.cpp`. On the icosahedron its strings equal DGGRID's for every cell (`paper/bench/bench_dggrid_z7.R`, fixture `tests/testthat/data/dggrid_z7.csv`); DGGRID's encoder is bijective, so do not diverge from it. Octahedral grids, for which IGEO7 defines nothing, write quad + 6 * seed.
- OpenEAGGR ISEA3H = already covered by `hex_grid(aperture = 3)`
- rHEALPix = diamond-based, not hexagonal — out of scope

No additional grid backends needed.

## Cell IDs

ISEA cell IDs are `bit64::integer64` in R and `int64_t` in C++ (`src/cell_id.h`):
a double is exact only below 2^53, which aperture 4 passes at resolution 25 and
aperture 7 at 18 (#99). Every Rcpp entry point taking IDs calls
`require_cell_ids()`, so a plain double never reaches the decoder as raw bits;
every R function taking IDs passes them through `as_cell_id()` (accepts
integer64, whole doubles below 2^53, digit strings). `quad_frame()` refuses a
grid whose cell count passes 2^63 - 1 (`grid_cell_count()`), the one source of
the resolution cap (`isea_cell_count()`, `isea_max_resolution()` in R).

integer64 loses its class silently in `unlist()`, `ifelse()`,
`vapply(numeric(1))`, `rep_len()`, `for (x in ids)`, `paste()`/`sprintf()` and
`c(<double>, <integer64>)`: use `cell_id_unlist()`, `seq_along()` loops,
`rep(length.out =)` and `as.character()`. `match`, `%in%`, `order`, `table`,
`factor` and the set functions are imported from bit64, so they are safe inside
the package.

## Z7 integer forms and neighbours

`src/index_z7.cpp` works on a `z7::Label` (leading field, resolution, digits);
the string, IGEO7's packed 64-bit index (base cell in bits 63-60, three bits
per digit, 7 below the cell's resolution, resolutions 0-20) and the monotonic
ID (base * 7^r + digits in base 7, resolutions 0-21) are forms of it.
`cell_to_index(form = "int" | "hex" | "monotonic")` and `index_to_cell()`
expose them on aperture-7 icosahedral grids. The packed index crosses into R as
integer64 with its 64 bits unchanged: base cells 08-11 read negative, and base
cell 08's resolution-20 pentagon (bits 2^63) is bit64's NA; `"hex"` is exact
for every index. The slides' Table 3 and DGGRID's INT64 output
(`tests/testthat/data/dggrid_z7_int.csv`) are the test vectors.

`z7::neighbor()` is generalized balanced ternary addition on a label: digit
sum and carry tables built from `downAp7`/`downAp7r` (they equal H3's
NEW_DIGIT_II/III and NEW_ADJUSTMENT_II/III), a 12 x 6 base-cell table
(`kBaseNeighbor`, `kBaseTurns`) for a carry out of resolution 1, and the 60-degree
turn across a pentagon's seam. `neighbors_in_frame()` takes an aperture-7
icosahedral cell's neighbours from it whenever a lattice step leaves the quad;
every other grid sends such a step through lon/lat. The two agree on every cell
at resolutions 0-6 (`tests/testthat/test-z7-neighbors.R`, the lon/lat path
reached through `cpp_get_neighbors_isea_lonlat()`). Apertures 3 and 4 keep the
lattice: Z3 and Z-order digits are the base-3/base-2 digits of the quad's
(i, j), so digit arithmetic there is the integer step done in O(r).

DGGRID built on Windows reads an input SEQNUM above 2^32 modulo 2^32
(`sscanf("%lu")` into a 32-bit `unsigned long`, `SubOpBasicMulti.cpp`,
`SubOpGenHelper.cpp`); its SEQNUM output is right. Benches feeding DGGRID
SEQNUMs stay below 2^32 or go through GEO input. Repro (git-ignored):
`dev_notes/dggrid_seqnum_input_2pow32.R`.

## Bodies

A grid is sized on any sphere: `hex_grid(radius_km = )` takes a radius in km or a
body name (`"mars"`, `"moon"`, `"titan"`, ...). The grid object carries the
radius in its `radius_km` slot, and every function reporting kilometres reads it
through `grid_radius_km()` rather than the Earth constant. Cell geometry is
angular and radius-free, so only areas, diagonals, spacings and the
resolution-for-area inversion change.

Both backends take a radius. The vendored H3 C library computes a cell's area as
its solid angle times `EARTH_RADIUS_KM` (`src/h3/area.c`), so hexify scales H3
areas at the R boundary with `scale_area_to_body()` -- by the square of the
radius ratio -- rather than touching vendored code. H3 cell IDs remain H3's
Earth-read topology; a grid on another body is that topology on that body, and
`h3_crosswalk()` needs both grids on the same body.

## Orientation and face projection

An ISEA grid carries its solid in the `polyhedron` slot (`"icosahedron"` or
`"octahedron"`), the solid's orientation in the `orientation` slot,
`c(vert0_lon, vert0_lat, azimuth)` (DGGRID's `dggs_vert0_*`), and its face
projection in the `projection` slot, `"isea"` (Snyder), `"ivea"` (van Leeuwen
and Strebe's vertex-oriented equal-area) or `"fuller"`. `FACE_PROJECTIONS`,
`EQUAL_AREA_PROJECTIONS` and `ICOSAHEDRON_PROJECTIONS` (`R/orientation.R`)
say which is which; R code asks `is_equal_area_projection()`, never a name.

`src/polyhedron.cpp` describes each solid (icosahedron, octahedron,
tetrahedron) by its vertices, faces and Snyder constants, and derives the
quad scheme from them (`SolidTopology`, read through `topo()`): quad q is
vertex q, quads 1..V-2 are diamonds of two faces, quads 0 and V-1 are
single-cell vertex quads, and every face pairing, region table and edge map
is derived, not hand-written. Code below the projection reads counts and
tables from `topo()`, never 20/12/10. The tetrahedron has no quad scheme
(`has_quads` false: its four valence-3 vertices would all have to be vertex
quads), so it carries the projection only.

The C++ layer keeps one face table per (solid, orientation) and reads the
active one through `poly()`; `project_core()` and `face_xy_to_ll()` read the
active projection (`active_projection()`). Every Rcpp entry point that reaches
either takes `icosa` as its FIRST argument and calls `activate_icosa(icosa)`,
or `activate_grid(icosa)` where it needs cells (it stops on a solid without a
grid) (`src/rcpp_icosa.h`). R passes `icosa_arg(g)` =
`c(orientation, projection code, solid code)` (plus the flattening for a
grid on an ellipsoid, see Earth model), `numeric(0)` for the
icosahedron's default orientation that `hexify_build_icosa()` sets on ISEA,
`projection_icosa(projection, polyhedron)` = `c(projection, solid)` for a
solid's default orientation, or `standard_icosa(polyhedron)`. Test-only entry
points call `activate_default_icosa()` instead. A new entry point reaching
`poly()`, `topo()` or a projection must do one of these, or it reads
whatever the previous call left.

Fuller (`src/projection_fuller.cpp`) is written from Gray (1995) and Crider
(2008), DGGRID's `DgProjFuller` used as a reference only. Both projections
share the face frame: azimuth measured from the face's first vertex, plane
triangle of unit edge, so everything downstream of the face coordinates is
projection-free. Fuller is not equal-area: `cell_area()` sums each cell's
solid angle (`cpp_cell_solid_angle`).

IVEA (`src/projection_ivea.cpp`) is written from van Leeuwen and Strebe
(2006), DGGAL's `icoVertexGreatCircle.ec` used as a reference only. It works
on any triangular face from `SolidTopology::vgc`, so it runs on every solid;
its inverse is closed form, with a 2D Newton on the forward kept as the
checker. It creases on the lines from a face's centre to its corners and to
its edge midpoints, so face pieces (`radius_crossings()`) and the globe mesh
(`n_face_sectors()`) are cut at six lines, against Snyder's three.
The cell-ID hierarchy is orientation-free, so the mixed-aperture hierarchy
(`R/aperture_mixed_hierarchy.R`) runs in the standard orientation.

## Globe widget: cells found per pixel, values in a table

`hex_globe()` (and `hex_globe_png()` through the companion package hexglobe,
`C:/Users/GillesC/Documents/dev/hexglobe`, local git) draw ISEA cells in
`inst/wgsl/globe.wgsl`: the fragment shader projects the pixel to its face,
places it in the face's quad and takes the nearest multiple of the
generator (`locate()`). hexglobe draws with the shader the widget carries,
so the two stay in step through `globe_scene()`; a change to the bindings
changes `src/rust/lib.rs` and `src/init.c` there too.

Values are not read by cell ID. `cpp_globe_table()` (`src/globe_table.cpp`)
lays every resolution from the grid's down to 0 out as one block per
diamond quad: row u holds substrate column u, padded past the box with the
cells that own those points across the quad's edges (the inverse of
`canonicalize_q2d()`'s edge maps). The shader reads a cell's word where it
finds the cell, with no move into the owning quad; `spot_id()` still makes
that move for the readout's ID. The table is slot by slot when the cells
fill it, else packed by a perfect hash (key mod M + offset[bucket]) with the
key stored; `mix32()` must match in C++ and WGSL. A coarser level is gathered
from the finer block through a fixed stencil when slot by slot, pushed from
the given cells when hashed; a test asserts the two agree cell by cell, and
the solid's vertex cells gather from every quad around them. Level 0 words are f32
values (NA one NaN, absent 0xFFFFFFFF); coarser levels pack ramp position
and cover in 16 bits each, and the shader fills from the finest level whose
cells are at least a pixel wide. R builds the Grid uniform once
(`globe_grid_uniform()`), the widget and hexglobe only upload it.

## Earth model: sphere or ellipsoid

By default latitude is read on the sphere (as DGGRID and H3 do); on WGS84
the cells' areas then run 0.9955 to 1.0090 of the reported area
(`paper/bench/bench_ellipsoid_area.R`). `hex_grid(ellipsoid = )` ("WGS84",
"GRS80", "earth", the IAU spheroids "mars", "jupiter", "saturn", "uranus",
"neptune", or `c(a_km, f)`) stores `c(a_km, f)` in the `ellipsoid` slot,
sizes the grid on the authalic radius (`radius_km`), and converts geodetic
to authalic latitude where lon/lat enter the projection and back where they
leave it, so ISEA and IVEA cells are equal-area on the ellipsoid.

- C++: `src/authalic.cpp`, one implementation. `use_ellipsoid(f)` is set by
  `activate_icosa()` from a sixth icosa element (the flattening;
  `icosa_arg()` appends it only when f > 0, so `ellipsoid = NULL` keeps the
  five-element argument and is bit-identical to a sphere-only build).
  `snyder_forward()`, `snyder_forward_to_face()`, `which_face()` and
  `face_xy_to_ll()` take and give geodetic latitude; code that stays on the
  sphere uses `snyder_forward_sphere()` and `face_xy_to_sphere_ll()`. The
  conversion is a sine series fitted at activation to the closed form in
  q(phi) (`authalic_lat_exact()`, Newton inverse `geodetic_lat_exact()`,
  kept as the checker, `cpp_sphere_latitude(exact = TRUE)`).
- R: code doing sphere geometry on lon/lat (sampling, drawing, nets, globe)
  converts through `unit_vec(lon, lat, icosa)`, `vec_lonlat(P, icosa)`,
  `sphere_lat()` and `geodetic_lat()`. The orientation is on the sphere
  (`vert0_lat` authalic); a `region` is given geodetic.
- `orientation = "ogc"` is OGC's and DGGAL's (11.20E, authalic arctan(phi)).
  OGC's ISEA3H/ISEA7H/IVEA3H/IVEA7H = `hex_grid(aperture = 3 or 7,
  projection = "isea" or "ivea", ellipsoid = "WGS84", orientation = "ogc")`.

## OGC zone identifiers (textZIRS, uint64ZIRS)

`cell_to_index(form = "textZIRS" | "uint64ZIRS")` and `index_to_cell()`
(`R/ogc_zirs.R`) write and read the identifiers DGGAL gives OGC's ISEA3H,
ISEA7H, IVEA3H and IVEA7H, for any aperture-3/7 icosahedron grid. Root
rhombus r is hexify quad r/2 + 1 (even r) or (r - 1)/2 + 6 (odd), 10 and 11
the vertex quads 0 and 11; sub-rhombi come from the quad-plane centre. Odd
aperture-7 letters follow DGGAL's own naming of pentagon, polar and
first-row/first-column parents (see the file's header). They equal DGGAL's
for every zone up to ISEA3H level 10 and ISEA7H level 5
(`paper/bench/bench_dggal_zirs.R`); do not diverge from DGGAL there.

## DGGRID corner bug: aperture 7, odd resolutions

Found 2026-10-05; not filed upstream yet. At aperture 7, odd resolutions
(Class III), DGGRID's corners for a cell straddling a face edge are wrong.
Its output then has gaps and overlaps; hexify's does not. hexify is correct
here; do not "fix" hexify towards DGGRID.

- Size: 140 of 3,432 cells at res 3 (up to 5.5 km), 980 of 168,072 at res 5
  (up to 1.4 km), two corners per cell. Cell assignment and centres agree.
- Proof: every hexify corner is shared by exactly three cells. Each affected
  DGGRID cell has two corners no other DGGRID cell carries, and its two
  DGGRID neighbours put that corner where hexify does (139/140 cells at
  res 3, 980/980 at res 5).
- Cause: `DgIDGGBase::setAddVertices` → `DgQ2DDtoVertex2DDConverter`
  (`DgIDGGutil.cpp:632`) chooses the face from the cell's own quad, so a
  corner past the quad's far edge is computed on the near face, with the
  inverse projection extended past its edge. hexify's `quad_point_lonlat()`
  (`src/rcpp_cell.cpp`) first moves the point into the quad that owns it.
- Effect on benches: `bench_dggrid_agreement.R` tests aperture 7 only at even
  resolutions, so it never shows this. The orientation bench tests res 5, so
  its corner-gap column picks up DGGRID's error, not hexify's.
- Local repro (dev_notes/ is git-ignored): `dev_notes/ap7_odd_corner_repro.R`,
  `dev_notes/ap7_odd_corner_share.R`,
  `dev_notes/dggrid_ap7_odd_corners_2026-10-05.md`.

## Vendored H3 C Library

The H3 backend uses vendored C source from Uber's H3 library in `src/h3/`. This is a direct copy of upstream code — **do not rewrite or restyle it**. Benefits:
- Zero external dependencies at install time
- Easy to update: drop in new upstream C files when Uber releases a new version
- Apache 2.0 license (compatible with MIT), kept at `src/h3/LICENSE`

Only make minimal targeted fixes when R CMD check flags specific symbols (e.g., `sprintf` → `snprintf`). Never refactor vendored code for style.

## Build & Check

Use Windows R (not WSL R):
```bash
"/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe" -e 'devtools::document()'
"/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe" -e 'devtools::check(args = "--no-manual")'
```

pkgdown site:
```bash
"/mnt/c/Program Files/R/R-4.6.0/bin/Rscript.exe" -e 'source("~/.R/build_pkgdown.R"); build_pkgdown_site()'
```

## Git Push

Use Windows git/gh for remote operations (WSL SSH agent doesn't persist):
```bash
cmd.exe /C "cd /d C:\Users\Gilles Colling\documents\dev\hexify && git push origin main"
cmd.exe /C "cd /d C:\Users\Gilles Colling\documents\dev\hexify && gh release create v0.x.x --title \"title\" --notes \"notes\""
```

## SoftwareX paper

The manuscript lives in `C:/GillesC/Documents/writing/papers/paper_hexify_2026` (own git repo; see its `CLAUDE.md`). Scripts behind any number in the paper go in this repository under `paper/bench/`, so the tagged release carries them.

## University of Vienna and open-source software

Recorded 2026-10-05 from university documents and Austrian law read that day.

- **Rights:** UrhG §40b gives the employer an unrestricted exploitation right to
  computer programs an employee writes "in Erfüllung seiner dienstlichen
  Obliegenheiten", unless agreed otherwise; the right to be named as author stays
  with the author. The university's RDM Policy (Rectorate, 8 Sept 2021) §3:
  usage rights to research data lie "im Regelfall" with the university, citing
  §40b for programs. Whether hexify counts as a job duty depends on the predoc
  contract (not checked).
- **Open release is allowed:** RDM Policy §3 entitles staff to publish research
  data under open licences in repositories unless legal, contractual, ethical or
  other reasons (such as the university's commercial interests) stand against it.
  The RDM FAQ (v1, 5 Oct 2021) counts "software and code" as research data and
  names GPL as an example licence; MIT is not mentioned. The Technology Transfer
  Office's Guide to Patents (printed June 2024, p. 18): "Can I publish the
  software open source? Yes, the University supports an open-source policy",
  recommending a check of third-party requirements, all authors' contributions
  and incoming licences (for hexify: the Apache-2.0 H3 code).
- **No approval step** is stated in any public document. The "open-source
  policy" the guide refers to was not found publicly (may be intranet-only).
- **Reporting:** UG 2002 §106(3) requires service inventions to be reported to
  the Rectorate; software as such is not named. The invention-disclosure form has
  a software section; a separate Software Disclosure Form existed but its URL
  now returns 404.
- **Archiving:** the RDM Policy names PHAIDRA as an example repository and asks
  for persistent identifiers. Zenodo is not mentioned.
- **LICENSE copyright line:** no university rule found. "Copyright (c) Gilles
  Colling" is consistent with the authorship right; whether the university wants
  to be named is open.
- **Contacts (from the documents):** Technology Transfer Office
  techtransfer@univie.ac.at (Guide to Patents p. 23), now under
  https://forschungsservice.univie.ac.at/transfer; research data management
  rdm@univie.ac.at (RDM FAQ).
- Sources:
  https://www.ris.bka.gv.at/NormDokument.wxe?Abfrage=Bundesnormen&Gesetzesnummer=10001848&Paragraf=40b ;
  https://www.ris.bka.gv.at/NormDokument.wxe?Abfrage=Bundesnormen&Gesetzesnummer=20002128&Paragraf=106 ;
  https://datenmanagement.univie.ac.at/fileadmin/user_upload/p_forschungsdatenmanagement/Dokumente/RDM_Policy_UNIVIE_v1_de.pdf ;
  https://rdm.univie.ac.at/fileadmin/user_upload/p_forschungsdatenmanagement/Dokumente/RDM_FAQ__UNIVIE_v1_en.pdf ;
  https://forschungsservice.univie.ac.at/fileadmin/user_upload/forschungsservice/Dokumente/Guide_to_Patents.pdf
- Open for the user: MIT acceptable? notice to the TTO needed? copyright line?
  any IP clause in the predoc contract?
