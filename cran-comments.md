## R CMD check results

0 errors | 0 warnings | 1 note

* This is a hotfix release addressing geometry issues in global grids.

The NOTE is "Days since last update" which is expected for a patch release.

## Changes in v0.3.10

1. Fixed invalid pentagon geometries that caused visual gaps in global grids
2. Fixed antimeridian-crossing polygons using `st_wrap_dateline()`
3. Added polar cap sampling to `grid_global()` to include cells above ±85° latitude

## Test environments

* local Windows 11, R 4.5.2 (tests: ~13s)
* GitHub Actions (ubuntu, macOS, windows)

## Downstream dependencies

There are currently no downstream dependencies for this package.
