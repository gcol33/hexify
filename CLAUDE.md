# hexify Development Notes

## Grid Coverage

hexify supports **all** major hexagonal DGGS through two backends:

- **ISEA** (built-in C++): apertures 3, 4, 7, and any mixed sequence of them — resolutions 0-30.
  `hex_grid(aperture = "4/3")` / `"4/7"` / `"7/4"` name a family (first `floor(res/2)` levels
  take the first aperture), and `aperture = c(4, 4, 7, 3)` names one aperture per level.
- **H3** (vendored H3 v4.4.1 C source in `src/h3`): fixed aperture 7 — resolutions 0-15

This covers every hexagonal grid system that matters:
- ISEA3H, ISEA4H, ISEA7H, ISEA43H = ISEA backend with different aperture settings
- H3 = H3 backend
- IGEO7/Z7 (2025, Sahr) = equal-area aperture-7 hex grid with Z7 indexing = already covered by `hex_grid(aperture = 7)` (ISEA7H with Z7 index). hexify uses the same Z7 hierarchical indexing (7-digit encoding, 0-6 per level) in `src/index_z7.cpp`.
- OpenEAGGR ISEA3H = already covered by `hex_grid(aperture = 3)`
- rHEALPix = diamond-based, not hexagonal — out of scope

No additional grid backends needed.

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

## Orientation

An ISEA grid carries its icosahedron's orientation in the `orientation` slot,
`c(vert0_lon, vert0_lat, azimuth)` (DGGRID's `dggs_vert0_*`). The C++ layer
keeps one face table per orientation (`src/icosahedron.cpp`) and reads the
active one through `ico()`. Every Rcpp entry point that reaches `ico()` takes
`orient` as its FIRST argument and calls `activate_orientation(orient)`
(`src/rcpp_orientation.h`) on entry; R passes `orient_arg(g)`, or
`numeric(0)` for the default orientation that `hexify_build_icosa()` sets.
Projection-level and test-only entry points call
`hexify::use_default_orientation()` instead. A new entry point reaching
`ico()` must do one of the two, or it reads whatever the previous call left.
The cell-ID hierarchy is orientation-free, so the mixed-aperture hierarchy
(`R/aperture_mixed_hierarchy.R`) runs in the standard orientation.

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
