# Create a Hexagonal Grid Specification

Creates a HexGridInfo object that stores all parameters needed for
hexagonal grid operations. Use this to define the grid once and pass it
to all downstream functions.

## Usage

``` r
hex_grid(
  area_km2 = NULL,
  resolution = NULL,
  aperture = 3,
  type = c("isea", "h3"),
  resround = "nearest",
  crs = 4326L
)
```

## Arguments

- area_km2:

  Target cell area in square kilometers. Mutually exclusive with
  `resolution`.

- resolution:

  Grid resolution level (0-30 for ISEA, 0-15 for H3). Mutually exclusive
  with `area_km2`.

- aperture:

  Grid aperture: 3 (default), 4, 7, or "4/3" for mixed. Ignored for H3
  grids (fixed at 7).

- type:

  Grid type: "isea" (default) or "h3". H3 grids require the h3o package.

- resround:

  Resolution rounding when using `area_km2`: "nearest" (default), "up",
  or "down".

- crs:

  Coordinate reference system EPSG code (default 4326 = 'WGS84').

## Value

A HexGridInfo object containing the grid specification.

## Details

Exactly one of `area_km2` or `resolution` must be provided.

When `area_km2` is provided, the resolution is calculated automatically
using the cell count formula: N = 10 \* aperture^res + 2 (ISEA) or by
matching the closest H3 resolution.

H3 grids use the Uber H3 hierarchical hexagonal system. Unlike ISEA
grids, H3 cells are NOT exactly equal-area (area varies by ~3-5\\
location). H3 support requires the h3o package.

## One Grid, Many Datasets

A HexGridInfo acts as a shared spatial reference system - like a CRS,
but discrete and equal-area. Define the grid once, then attach multiple
datasets without repeating parameters:


    # Step 1: Define the grid once
    grid <- hex_grid(area_km2 = 1000)

    # Step 2: Attach multiple datasets to the same grid
    birds <- hexify(bird_obs, lon = "longitude", lat = "latitude", grid = grid)
    mammals <- hexify(mammal_obs, lon = "lon", lat = "lat", grid = grid)
    climate <- hexify(weather_stations, lon = "x", lat = "y", grid = grid)

    # No aperture, resolution, or area needed after step 1 - the grid
    # travels with the data.

    # Step 3: Work at the cell level
    # Once hexified, lon/lat no longer matter - cell_id is the shared key
    bird_counts <- aggregate(species ~ cell_id, data = as.data.frame(birds), length)
    mammal_richness <- aggregate(species ~ cell_id, data = as.data.frame(mammals),
                                 function(x) length(unique(x)))

    # Join datasets by cell_id - guaranteed to align because same grid
    combined <- merge(bird_counts, mammal_richness, by = "cell_id")

    # Step 4: Visual confirmation
    # All datasets produce identical grid overlays
    plot(birds)   # See the grid
    plot(mammals) # Same grid, different data

## See also

[`hexify`](https://gcol33.github.io/hexify/reference/hexify.md) for
assigning points to cells,
[`HexGridInfo-class`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
for class documentation

## Examples

``` r
# Create grid by target area
grid <- hex_grid(area_km2 = 1000)
print(grid)

# Create grid by resolution
grid <- hex_grid(resolution = 8, aperture = 3)

# Create grid with different aperture
grid4 <- hex_grid(area_km2 = 500, aperture = 4)

# Create mixed aperture grid
grid43 <- hex_grid(area_km2 = 1000, aperture = "4/3")

# Use grid in hexify
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
```
