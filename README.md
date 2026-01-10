# hexify

[![CRAN status](https://www.r-pkg.org/badges/version/hexify)](https://CRAN.R-project.org/package=hexify)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/hexify)](https://cran.r-project.org/package=hexify)
[![R-CMD-check](https://github.com/gcol33/hexify/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gcol33/hexify/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/hexify/graph/badge.svg)](https://app.codecov.io/gh/gcol33/hexify)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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

- [Quick Start](https://gillescolling.com/hexify/articles/quickstart.html) - Basic usage and concepts
- [Visualization](https://gillescolling.com/hexify/articles/visualization.html) - Plotting options
- [Workflows](https://gillescolling.com/hexify/articles/workflows.html) - Grid generation, spatial joins, multi-resolution analysis

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools should be free and open. I started these projects for my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT
