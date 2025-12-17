# Convert effective resolution to index resolution

For aperture 3: res_index = 2 \* eff_res - 1 for odd index resolutions.

## Usage

``` r
hexify_eff_res_to_resolution(eff_res)
```

## Arguments

- eff_res:

  Effective resolution

## Value

Index resolution (integer)

## See also

Other grid statistics:
[`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md),
[`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md),
[`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md),
[`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md),
[`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md),
[`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md),
[`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)
