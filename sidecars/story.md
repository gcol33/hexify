# Story: hexify

## Claim
Assigns points on the globe to equal-area hexagonal grid cells (ISEA), with H3 grids for interoperability.

## Audience
Ecologists, epidemiologists and geographers who bin point data into spatial cells in R, mostly with sf; some know dggridR or H3.

## Relation
- to: dggridR (ISEA grids, DGGRID); H3 (bundled C library)
- phrase: the discrete global grid system of Sahr, White and Kimerling (2003) that dggridR also provides
- same: ISEA projection and cell numbering; checked against dggridR reference tables in tests/testthat/test-dggrid-compat.R
- differs: computed in C++ in the package with no external libraries; H3 bundled; other planetary radii
- home: vignettes/theory.Rmd (method); reference index, dggridR Compatibility group (conversion)

## Arcs

### A1 Equal area at any cell size
- question: will my cell counts be comparable from equator to pole?
- answer: ISEA cells have one area per resolution; hex_grid() picks the resolution closest to a requested area
- evidence: vignettes/theory.Rmd; tests
- limits: twelve pentagons per grid
- home: vignettes/theory.Rmd
- also: README, quickstart

### A2 Points to cells to maps
- question: how do I get from a table of coordinates to a map per cell?
- answer: hexify(), hex_grid(), cell_to_sf(), grid_rect()/grid_clip(), plot()/hexify_heatmap()
- home: vignettes/quickstart.Rmd
- also: README Usage, workflows.Rmd, visualization.Rmd

### A3 H3 grids alongside ISEA
- question: can I use H3 cells and compare them to ISEA?
- answer: bundled H3 v4.4.1, h3_crosswalk()
- limits: H3 cells are not equal-area; aperture 7 only, resolution up to 15
- home: vignettes/h3.Rmd
- also: README

### A4 Same grid as dggridR
- question: will my cell IDs match dggridR?
- answer: test suite compares against dggridR reference tables; as_dggrid()/from_dggrid()
- evidence: tests/testthat/test-dggrid-compat.R
- home: reference index (dggridR Compatibility)
- also: README trust paragraph

## Non-goals
None stated in the docs.

## Documents
- DESCRIPTION: claim, ISEA, H3, Hex9, dggridR compatibility
- README / landing: claim, relation, first example, A1-A4 in one paragraph or section each, routing
- Get Started (quickstart): one path from points to a map
- articles: workflows, visualization, globe (how-to); h3 (how-to); theory, lambert, ogc-conformance (explanation)

## Open
- Whether the README should carry a Hex9 section once its article exists.
