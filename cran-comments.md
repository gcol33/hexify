## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission (v0.8.1)

The version on CRAN is 0.6.5. This release carries the work published on GitHub
since then, which NEWS.md lists in full. The parts that reach a user's stored
results:

* Aperture-7 cell IDs and 'Z7' index strings changed, and `hexify_assign()` now
  names the cell a point falls in. Aperture-3 and aperture-4 cell IDs are
  unchanged. NEWS.md lists these under Breaking changes.

* A grid is sized on any body: `hex_grid(radius_km = )` takes a radius in
  kilometres or a body name, and the grid carries the coordinate reference
  system of the body it sits on.

* Apertures 3, 4 and 7 in any sequence, named as a family ("4/7") or one
  aperture per resolution level.

Uber Technologies, Inc. is now recorded in Authors@R as contributor and
copyright holder for the H3 C library included in src/h3, which inst/COPYRIGHTS
documents with its file paths and license.

## Test environments

* local: Windows 11, R 4.6.0, R CMD check --as-cran
* GitHub Actions: ubuntu-latest (release, devel), macos-latest (release)

## Downstream dependencies

There are currently no downstream dependencies for this package.
