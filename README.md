# hexify

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/gcol33/hexify/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/hexify/actions/workflows/R-CMD-check.yaml)

**Equal-area hexagonal grids for spatial analysis**

Assign geographic points to hexagonal cells using the ISEA discrete global grid system. All cells have the same area regardless of latitude.

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/hexify")
```

## Quick Start

```r
library(hexify)

# Your data
cities <- data.frame(
  name = c("Vienna", "Paris", "Madrid"),
  lon = c(16.37, 2.35, -3.70),
  lat = c(48.21, 48.86, 40.42)
)

# Create a grid and assign points
grid <- hex_grid(area_km2 = 10000)
result <- hexify(cities, lon = "lon", lat = "lat", grid = grid)

# Visualize
plot(result)
```

## Key Features

- **Equal-area cells**: No latitude distortion
- **Simple workflow**: Define grid once, reuse everywhere
- **Fast C++ core**: Handles millions of points
- **sf integration**: Works with any CRS
- **dggridR compatible**: Same cell IDs for interoperability

## Learn More

- [Quick Start](https://gcol33.github.io/hexify/articles/quickstart.html) - Basic usage and concepts
- [Visualization](https://gcol33.github.io/hexify/articles/visualization.html) - Plotting options
- [Workflows](https://gcol33.github.io/hexify/articles/workflows.html) - Grid generation, spatial joins, multi-resolution analysis

## License

MIT
