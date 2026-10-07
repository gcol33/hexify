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
- H3 = H3 backend
- IGEO7/Z7 (2025, Sahr) = equal-area aperture-7 hex grid with Z7 indexing = already covered by `hex_grid(aperture = 7)` (ISEA7H with Z7 index). hexify uses the same Z7 hierarchical indexing (7-digit encoding, 0-6 per level) in `src/index_z7.cpp`. On the icosahedron its strings equal DGGRID's for every cell (`paper/bench/bench_dggrid_z7.R`, fixture `tests/testthat/data/dggrid_z7.csv`); DGGRID's encoder is bijective, so do not diverge from it. Octahedral grids, for which IGEO7 defines nothing, write quad + 6 * seed.
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
projection in the `projection` slot, `"isea"` (Snyder) or `"fuller"`.

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
`c(orientation, projection code, solid code)`, `numeric(0)` for the
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
The cell-ID hierarchy is orientation-free, so the mixed-aperture hierarchy
(`R/aperture_mixed_hierarchy.R`) runs in the standard orientation.

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
