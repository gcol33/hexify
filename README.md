# hexify

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**Fast Hexagonal Grid Assignment for Geographic Data**

`hexify` assigns geographic points to equal-area hexagonal cells using the ISEA3H discrete global grid system. It produces output identical to `dggridR` but with a simpler, modern interface.

## Quick Start

```r
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
#>     site   lon   lat hex_id hex_cen_lon hex_cen_lat
#> 1 Vienna 16.37 48.21  12847    16.42035    48.26151
#> 2  Paris  2.35 48.86  12532     2.31894    48.89826
#> 3 Madrid -3.70 40.42  22178    -3.71892    40.38721
```

## Why hexify?

- **Simple API**: One function (`hexify()`) for the common use case
- **dggridR compatible**: Produces identical `hex_id` (SEQNUM) values
- **Modern R**: Works with data.frames and sf objects
- **Fast**: C++ implementation of ISEA Snyder projection

## Features

- **Point-to-cell assignment**: Given lon/lat, return cell ID and center
- **Area or spacing control**: Specify target cell area (km²) or spacing (km)
- **sf integration**: Handles any CRS, auto-transforms to WGS84
- **Aperture 3 grids**: ISEA3H compatible with dggridR

## Installation

```r
# Install from GitHub
remotes::install_github("gcol33/hexify")
```

## Usage

### Basic: Data Frame

```r
library(hexify)

df <- data.frame(
  lon = c(0, 10, 20),
  lat = c(0, 45, -30)
)

# By area (km²)
result <- hexify(df, lon = "lon", lat = "lat", area = 1000)

# By spacing (long diagonal in km)
result <- hexify(df, lon = "lon", lat = "lat", spacing = 50)
```

### With sf Objects

```r
library(sf)

# Any CRS works - auto-transforms to WGS84
pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
result <- hexify(pts, area = 1000)

# sf class preserved
class(result)
#> [1] "sf"         "data.frame"
```

### Output Columns

`hexify()` adds three columns to your data:

| Column | Type | Description |
|--------|------|-------------|
| `hex_id` | integer | Unique cell identifier (SEQNUM) |
| `hex_cen_lon` | numeric | Cell center longitude |
| `hex_cen_lat` | numeric | Cell center latitude |

## dggridR Compatibility

hexify produces identical cell assignments to dggridR:

```r
library(dggridR)
library(hexify)

# Same result from both packages
dggs <- dggridR::dgconstruct(res = 10, aperture = 3)
ref <- dggridR::dgGEO_to_SEQNUM(dggs, lon, lat)

result <- hexify(df, lon = "lon", lat = "lat", area = 863)  # res 10
all(result$hex_id == ref$seqnum)
#> TRUE
```

## Resolution Reference

| Resolution | Cells | Area (km²) | Spacing (km) |
|------------|-------|------------|--------------|
| 5 | 2,432 | 209,903 | 695 |
| 6 | 7,292 | 69,968 | 401 |
| 7 | 21,872 | 23,323 | 232 |
| 8 | 65,612 | 7,774 | 134 |
| 9 | 196,832 | 2,591 | 77 |
| 10 | 590,492 | 864 | 45 |
| 11 | 1,771,472 | 288 | 26 |
| 12 | 5,314,412 | 96 | 15 |

## Limitations

- Currently supports aperture 3 only (ISEA3H)
- No polygon/boundary generation (use dggridR for that)
- No neighbor/parent/child operations in main API

## See Also

- [dggridR](https://github.com/r-barnes/dggridR) - Full DGGRID functionality
- [DGGRID](https://github.com/sahrk/DGGRID) - Original C++ implementation

## License

MIT
