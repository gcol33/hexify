## R CMD check results

0 errors | 0 warnings | 1 note

* This is a hotfix release addressing long test runtime (>10 min).

The NOTE is "Days since last update: 1" which is expected for a hotfix resubmission.

## Changes in v0.3.6 (hotfix)

1. Reduced test runtime from 13 min to ~3 min by adding skip_on_cran() to
   detailed consistency tests (full tests still run locally via NOT_CRAN=true)
2. Fixed CRAN incoming check NOTE: "Overall checktime 15 min > 10 min"

## Test environments

* local Windows 11, R 4.5.2 (tests: 199s)
* mac-builder R-release (tests: 182s)
* R-hub ubuntu-release (passed)

## Downstream dependencies

There are currently no downstream dependencies for this package.
