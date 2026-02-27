# Convert effective resolution to area

Calculates approximate cell area for ISEA3H (aperture 3) grids. Based on
calibration: eff_res 10 = 863.8006 km^2.

## Usage

``` r
hexify_eff_res_to_area(eff_res)
```

## Arguments

- eff_res:

  Effective resolution (0.5 \* resolution for aperture 3)

## Value

Area in km^2

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)
