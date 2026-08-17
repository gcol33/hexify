# Summarize Data by Hex Cell

Aggregates data within each hexagonal cell, similar to
`dplyr::group_by(cell_id) |> summarize(...)`. Returns a data.frame with
one row per unique cell, including cell center coordinates and area.

## Usage

``` r
hex_summarize(hex_data, ..., .fns = NULL, geometry = FALSE)
```

## Arguments

- hex_data:

  A HexData object (from
  [`hexify()`](https://gillescolling.com/hexify/reference/hexify.md)).

- ...:

  Named summary expressions (tidyeval). Each expression is evaluated per
  cell group with the group's columns available by name. Examples:
  `mean_temp = mean(temperature)`,
  `n_species = length(unique(species))`,
  `max_elev = max(elevation, na.rm = TRUE)`. Provide either `...` or
  `.fns`, not both.

- .fns:

  Optional named list of functions for formula-style aggregation.
  Example: `.fns = list(mean_temp = ~mean(temperature))`. Provide either
  `...` or `.fns`, not both.

- geometry:

  Logical. If `TRUE`, attach cell center points as an sf geometry column
  (requires sf). Default `FALSE`.

## Value

A data.frame with columns:

- cell_id:

  Unique cell identifier

- cell_cen_lon, cell_cen_lat:

  Cell center coordinates

- cell_area_km2:

  Cell area in km^2

- n_points:

  Number of data points in this cell

- ...:

  User-defined summary columns

If `geometry = TRUE`, returns an sf object with POINT geometry.

## Details

The function works entirely in R (no C++ needed). It groups by `cell_id`
and evaluates the summary expressions within each group.

If no summary expressions are provided, returns cell counts only.

## See also

[`hexify()`](https://gillescolling.com/hexify/reference/hexify.md) for
creating HexData objects,
[`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md)
for finding neighboring cells

## Examples

``` r
# \donttest{
df <- data.frame(
  lon = runif(100, -10, 10),
  lat = runif(100, 40, 55),
  temperature = rnorm(100, 15, 5),
  species = sample(letters[1:5], 100, replace = TRUE)
)
hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 500)

# Count points per cell
hex_summarize(hd)

# Custom summaries
hex_summarize(hd, mean_temp = mean(temperature),
                  n_species = length(unique(species)))
# }
```
