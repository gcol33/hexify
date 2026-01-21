## R CMD check results

0 errors | 0 warnings | 1 note

* This is a resubmission after test runtime issues.

The NOTE is "Days since last update" which is expected for a resubmission.

## Changes in v0.3.7

1. Removed slow spatial operation tests (grid_global, grid_rect, grid_clip)
2. Reduced comprehensive roundtrip test combinations
3. Test runtime now ~13 seconds (previously exceeded 10 min limit)

## Test environments

* local Windows 11, R 4.5.2 (tests: ~13s)
* GitHub Actions (ubuntu, macOS, windows)

## Downstream dependencies

There are currently no downstream dependencies for this package.
