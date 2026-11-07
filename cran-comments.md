## Submission (v0.9.1)

This release fixes the failure of 0.9.0 on r-release-macos-arm64 shown on the
CRAN check page: the examples took 18 minutes there and the tests ran past the
30-minute limit. It comes 15 days after 0.9.0 for that reason.

The `?cell_metrics` example hung on one aperture-3 cell. Cell walls are halved
until two estimates of their length and of the area beside them agree, and the
area's rounding floor shrank with the piece. On arm64 a wall ending on a face
edge rounded onto the next face and left a piece of zero length, which never
passed and was halved 2^20 times. The floor is now absolute, as the length's
already was, and a wall ending on a face edge ends there. On an Apple M4 Pro
(macOS, R 4.5.3) all examples now run in 9 s and the tests in 58 s. The same
floor made cells finer than about resolution 18 slow on every platform; they
now take milliseconds. Perimeters change by at most 2e-13 relative.

The linux-arm64 additional issue is likely the same hang; its log was not
reachable to confirm.

Where the compiler fuses multiply-adds (arm64, x86-64 with FMA), two inlined
copies of one projection step could round their last digit differently, so a
test asserting that `projection_stages()` ends on the forward projection bit
for bit failed. Both now call one compiled function. Projected coordinates and
cell IDs are unchanged where the compiler does not fuse.

## R CMD check results

0 errors | 0 warnings | 1 note

The note comes from the incoming feasibility check:

* Days since last update: 15 (see above).

* It lists 'hexglobe' in Suggests as not in a mainstream repository.
  'hexglobe' draws hexify's globes to a PNG on the graphics card through
  'wgpu'; only `hex_globe_png(renderer = "wgpu")` uses it, after
  `requireNamespace()`, and it is served from https://gcol33.r-universe.dev,
  declared in Additional_repositories.

`hex9_warp_download()` saves a 19 MB data file from the 'libhex9' repository
to `tools::R_user_dir("hexify", "data")` when the user calls it; its example
runs only in an interactive session, and no other function writes outside
`tempdir()`.

## Test environments

* local: Windows 11, R 4.6.1, R CMD check --as-cran
* local: macOS arm64 (Apple M4 Pro), R 4.5.3, R CMD check --as-cran
  (examples 14 s, tests 56 s)
* win-builder: R-devel (2026-11-06 r90845), R-release (4.6.2)
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)

The mac builder returned "currently unavailable" on 2026-11-07, so the macOS
arm64 check ran on my own machine.

## Downstream dependencies

There are currently no downstream dependencies for this package.
