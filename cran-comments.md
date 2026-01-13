## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission addressing CRAN reviewer feedback.

The NOTE is "unable to verify current time" which is a network check issue, not a package problem.

## Changes since last submission

Addressed all reviewer feedback:

1. Fixed DESCRIPTION quoting (proper names without quotes, software with quotes)
2. Added \value tags to HexData-methods.Rd and HexGridInfo-methods.Rd
3. Removed examples from unexported functions (hexify_compare_resolutions, hexify_get_resolution)
4. Replaced all \dontrun{} with executable examples or removed where not applicable
5. Fixed par() reset in visualization vignette (now uses oldpar <- par(...); par(oldpar))

## Test environments

* local Windows 11, R 4.5.2
* win-builder R-devel (2026-01-08 r89292 ucrt)
* win-builder R-release

## Downstream dependencies

There are currently no downstream dependencies for this package.
