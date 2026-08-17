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
[`dg_closest_res_to_area()`](https://gillescolling.com/hexify/reference/dg_closest_res_to_area.md),
[`dgearthstat()`](https://gillescolling.com/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gillescolling.com/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_resolution()`](https://gillescolling.com/hexify/reference/hexify_eff_res_to_resolution.md),
[`hexify_resolution_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_resolution_to_eff_res.md)
