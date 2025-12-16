# hexify

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**Fast Hexagonal Grid Assignment for Geographic Data**

`hexify` assigns geographic points to equal-area hexagonal cells using
the ISEA discrete global grid system. It produces output identical to
`dggridR` but with a simpler, modern interface. Supports apertures 3, 4,
7, and mixed 4/3.

## Quick Start

``` r

library(hexify)

# Any data.frame with coordinates
df <- data.frame(
  site = c("Vienna", "Paris", "Madrid"),
  lon = c(16.37, 2.35, -3.70),
  lat = c(48.21, 48.86, 40.42)
)

# Assign to ~1000 km² hexagonal cells
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)
result
#>     site   lon   lat cell_id cell_cen_lon cell_cen_lat cell_area cell_diag
#> 1 Vienna 16.37 48.21   12847     16.42035     48.26151    863.94     31.58
#> 2  Paris  2.35 48.86   12532      2.31894     48.89826    863.94     31.58
#> 3 Madrid -3.70 40.42   22178     -3.71892     40.38721    863.94     31.58
```

## Why hexify?

- **Simple API**: One function
  ([`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md))
  for the common use case
- **dggridR compatible**: Produces identical `cell_id` (SEQNUM) values
- **Modern R**: Works with data.frames and sf objects
- **Fast**: C++ implementation of ISEA Snyder projection

## Features

- **Point-to-cell assignment**: Given lon/lat, return cell ID and center
- **Polygon generation**: Native C++ polygon generation (faster than
  dggridR)
- **Area or diagonal control**: Specify target cell area (km²) or
  diagonal (km)
- **sf integration**: Handles any CRS, auto-transforms to WGS84
- **Multiple apertures**: ISEA3H, ISEA4H, ISEA7H, and ISEA43H grids
  (apertures 3, 4, 7, and mixed 4/3)

## Installation

``` r

# install.packages("pak")
pak::pak("gcol33/hexify")
```

## Usage

### Basic: Data Frame

``` r

library(hexify)

df <- data.frame(
  lon = c(0, 10, 20),
  lat = c(0, 45, -30)
)

# By area (km²)
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# By diagonal (long diagonal in km)
result <- hexify(df, lon = "lon", lat = "lat", diagonal = 50)

# Different apertures
result_ap4 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 4)
result_ap7 <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = 7)

# Mixed aperture (ISEA43H)
result_mixed <- hexify(df, lon = "lon", lat = "lat", area = 1000, aperture = "4/3")
```

### With sf Objects

``` r

library(sf)

# Any CRS works - auto-transforms to WGS84
pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
result <- hexify(pts, area = 1000)

# sf class preserved
class(result)
#> [1] "sf"         "data.frame"
```

### Output Columns

[`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md) adds
five columns to your data:

| Column         | Type    | Description                     |
|----------------|---------|---------------------------------|
| `cell_id`      | integer | Unique cell identifier (SEQNUM) |
| `cell_cen_lon` | numeric | Cell center longitude           |
| `cell_cen_lat` | numeric | Cell center latitude            |
| `cell_area`    | numeric | Actual cell area in km²         |
| `cell_diag`    | numeric | Actual cell diagonal in km      |

## dggridR Compatibility

hexify produces identical cell assignments to dggridR:

``` r

library(dggridR)
library(hexify)

# Same result from both packages
dggs <- dggridR::dgconstruct(res = 10, aperture = 3)
ref <- dggridR::dgGEO_to_SEQNUM(dggs, lon, lat)

result <- hexify(df, lon = "lon", lat = "lat", area = 863)  # res 10
all(result$cell_id == ref$seqnum)
#> TRUE
```

## Aperture Options

| Aperture | Grid Type | Cell Count Formula | Use Case |
|----|----|----|----|
| 3 | ISEA3H | 10 × 3^res + 2 | Default, compatible with dggridR |
| 4 | ISEA4H | 10 × 4^res + 2 | Faster cell count growth |
| 7 | ISEA7H | 10 × 7^res + 2 | Densest grid per resolution |
| “4/3” | ISEA43H | 10 × 4^m × 3^(res-m) + 2 | Mixed aperture (m = mixed_aperture_level) |

## Resolution Reference (Aperture 3)

| Resolution | Cells     | Area (km²) | Diagonal (km) |
|------------|-----------|------------|---------------|
| 0          | 12        | 42,506,000 | 9,908         |
| 1          | 32        | 15,939,750 | 6,067         |
| 2          | 92        | 5,544,261  | 3,578         |
| 3          | 272       | 1,875,265  | 2,081         |
| 4          | 812       | 628,167    | 1,204         |
| 5          | 2,432     | 209,734    | 696           |
| 6          | 7,292     | 69,950     | 402           |
| 7          | 21,872    | 23,321     | 232           |
| 8          | 65,612    | 7,774      | 134           |
| 9          | 196,832   | 2,591      | 77            |
| 10         | 590,492   | 864        | 45            |
| 11         | 1,771,472 | 288        | 26            |
| 12         | 5,314,412 | 96         | 15            |

### Polygon Generation

Generate hexagon polygons for visualization:

``` r

library(sf)

# After running hexify(), get polygons for unique cells
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# Generate sf polygons
polys <- hexify_cell_to_sf(result$cell_id, resolution = 10, aperture = 3)
plot(st_geometry(polys), col = "lightblue", border = "blue")

# Or generate a grid over a region
grid <- hex_grid_rect(
  minlon = -10, maxlon = 20,
  minlat = 35, maxlat = 60,
  area = 5000
)
```

## See Also

- [dggridR](https://github.com/r-barnes/dggridR) - Full DGGRID
  functionality
- [DGGRID](https://github.com/sahrk/DGGRID) - Original C++
  implementation

## License

MIT
