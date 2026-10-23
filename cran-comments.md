## R CMD check results

0 errors | 0 warnings | 1 note

The note comes from the incoming feasibility check:

* It lists "octahedral" in DESCRIPTION as possibly misspelled. It names the
  octahedron the Hex9 grid is built on and is spelled as intended.

* It lists 'hexglobe' in Suggests as not in a mainstream repository.
  'hexglobe' draws hexify's globes to a PNG on the graphics card through
  'wgpu'; only `hex_globe_png(renderer = "wgpu")` uses it, after
  `requireNamespace()`, and it is served from https://gcol33.r-universe.dev,
  declared in Additional_repositories.

## Submission (v0.9.0)

The version on CRAN is 0.8.2. This release carries the work published on
GitHub since then (0.8.3, 0.8.4 and 0.9.0), which NEWS.md lists in full. The
parts that reach a user's code or stored results:

* ISEA cell IDs are 64-bit integers (`bit64::integer64`; bit64 is a new
  import), exact at every resolution. As doubles they were exact only below
  2^53, so aperture 4 from resolution 25 and aperture 7 from resolution 18
  returned IDs that decoded to other cells. Functions taking cell IDs also
  accept whole numbers below 2^53 and digit strings. A grid's resolution is
  capped where its cell count passes 2^63 - 1.

* The aperture-7 Z7 index strings on the icosahedron now equal the ones
  DGGRID writes, and a grid with one aperture per level stores its sequence
  comma-separated ("4,4,7,3"). Cell IDs are unchanged by both.

* ISEA cell polygons follow the true, densified cell edges, and at odd
  aperture-7 resolutions a point goes to the cell whose centre is nearest,
  which changes the cell ID of points near cell edges there.

* `as_sf()` is replaced by methods on `sf::st_as_sf()`, `plot_globe()` by
  `plot(<grid>)`, and a few low-level helpers are removed or internal; NEWS.md
  names each.

* New: Fuller's and van Leeuwen and Strebe's face projections, grids on the
  octahedron, the Hex9 octahedral grid (Griffin 2026, address tables from
  'libhex9' under Apache-2.0, recorded in Authors@R and inst/COPYRIGHTS),
  orientations, grids on an ellipsoid, OGC zone identifiers, cell metrics and
  an interactive WebGPU globe (`hex_globe()`).

`hex9_warp_download()` saves a 19 MB data file from the 'libhex9' repository
to `tools::R_user_dir("hexify", "data")` when the user calls it; its example
runs only in an interactive session, and no other function writes outside
`tempdir()`.

## Test environments

* local: Windows 11, R 4.6.1, R CMD check --as-cran
* win-builder: R-devel (2026-10-22 r90748), R-release (4.6.1)

## Downstream dependencies

There are currently no downstream dependencies for this package.
