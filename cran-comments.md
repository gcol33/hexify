## R CMD check results

0 errors | 0 warnings | 1 note

* This is a hotfix release addressing the clang-UBSAN issue identified in v0.3.4.

The NOTE is "Days since last update: 0" which is expected for a hotfix resubmission.

## Changes in v0.3.5 (hotfix)

1. Fixed clang-UBSAN NaN casting error in Snyder projection C++ code
2. Simplified plot examples to reduce runtime (fixes slow example NOTE)

## Test environments

* local Windows 11, R 4.5.2
* win-builder R-devel (2026-01-20)
* R-hub clang-asan (passed)
* R-hub gcc-asan (passed)

## Downstream dependencies

There are currently no downstream dependencies for this package.
