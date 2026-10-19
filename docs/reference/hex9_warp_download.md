# Download Hex9's Warp Field

Hex9's own projection (`hex_grid(aperture = 9, projection = "akw")`) is
Kaseorg's octahedral projection with an area-correcting warp that
Griffin trained by optimal transport and publishes with 'libhex9', his
reference implementation (Griffin 2026, Sec. 11b). The field is 19 MB of
stored displacements and gradients, which hexify reads but does not
ship; this function downloads it once from the 'libhex9' repository,
checks it, and saves it where hexify looks for it.

## Usage

``` r
hex9_warp_download(path = hex9_warp_path(), quiet = FALSE)
```

## Arguments

- path:

  Where to save the field. The default is hexify's user data directory
  (`tools::R_user_dir("hexify", "data")`), or the option
  `hexify.hex9_warp` when set; hexify reads it from the same place.

- quiet:

  Passed to
  [`utils::download.file()`](https://rdrr.io/r/utils/download.file.html).

## Value

The path of the saved file, invisibly.

## Details

The file is 'libhex9”s `core/Sphere_l6_fund.f64g.h9warp` at the commit
of version 2.3.1, which the Hex9 agreement checks in hexify were run
against. 'libhex9' is Apache-2.0 licensed. Its field is defined on the
unit sphere: 'libhex9”s WGS84 functions first turn geodetic latitude
into authalic latitude, while hexify reads latitude as spherical
latitude, as 'libhex9”s `_sphere` functions do.

## References

Griffin, B. (2026). Hex9: A Quasi-Authalic, Quasi-Continuous Hexagonal
DGGS on the Reference Ellipsoid. arXiv:2608.00022.

## Examples

``` r
if (interactive()) {
  hex9_warp_download()
  g <- hex_grid(resolution = 5, aperture = 9, projection = "akw")
  lonlat_to_cell(16.37, 48.21, g)
}
```
