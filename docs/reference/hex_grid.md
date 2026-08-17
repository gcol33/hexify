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
  crs = 4326L,
  radius_km = EARTH_RADIUS_KM
)
```

## Arguments

- area_km2:

  Target cell area in square kilometers. Mutually exclusive with
  `resolution`.

- resolution:

  Grid resolution level (0-30 for ISEA, 0-15 for H3). Mutually exclusive
  with `area_km2`. For H3, typical use cases by resolution:

  - 0-3: continental/country scale

  - 4-7: regional/city scale

  - 8-10: neighborhood/block scale (FCC uses 8-9)

  - 11-15: building/sub-meter scale

- aperture:

  Grid aperture: 3 (default), 4, 7, a mixed family such as "4/3", "4/7"
  or "7/4", or one aperture per resolution level as a vector, e.g.
  `c(4, 4, 7, 3)`. A family name refines by the first aperture for the
  first `floor(resolution / 2)` levels and by the second for the rest,
  which is how DGGRID arranges ISEA43H. A per-level vector needs
  `resolution` rather than `area_km2`. Ignored for H3 grids (fixed at
  7).

- type:

  Grid type: "isea" (default) or "h3".

- resround:

  Resolution rounding when using `area_km2`: "nearest" (default), "up",
  or "down".

- crs:

  Coordinate reference system EPSG code (default 4326 = 'WGS84').

- radius_km:

  Radius of the body the grid covers, in kilometers, or the name of a
  body: "mercury", "venus", "earth" (default), "moon", "mars", "ceres",
  "jupiter", "io", "europa", "ganymede", "callisto", "saturn",
  "enceladus", "titan", "uranus", "neptune", "pluto".

## Value

A HexGridInfo object containing the grid specification.

## Details

Exactly one of `area_km2` or `resolution` must be provided.

When `area_km2` is provided, the resolution is calculated automatically
using the cell count formula: N = 10 \* aperture^res + 2 (ISEA) or by
matching the closest H3 resolution.

H3 grids use the Uber H3 hierarchical hexagonal system. Unlike ISEA
grids, H3 cells are NOT exactly equal-area (area varies by ~3-5\\
location).

## Other Bodies

A grid is a partition of the sphere, and `radius_km` sets the sphere it
is measured on. Cell geometry – which cell a coordinate lands in, where
cell centres and corners sit, the hierarchy, the neighbours – is angular
and identical on every body; the radius sets the kilometer figures: cell
area, diagonal, spacing, and the resolution that `area_km2` picks.
Earth's area comes from the 'WGS84' ellipsoid, every other radius gives
the sphere area 4*pi*r^2.


    mars <- hex_grid(area_km2 = 1000, radius_km = "mars")
    hex_grid(resolution = 8, radius_km = 3389.5)   # the same grid

Both backends take a radius. 'H3' reports a cell's area as its solid
angle times Earth's radius squared, so another radius scales those areas
by the square of the radius ratio, exactly. One caveat carries: an 'H3'
cell ID names a position in 'H3”s topology, which 'Uber”s 'H3' reads on
Earth, so the IDs of a grid on another body are that topology on that
body and are not interchangeable with Earth 'H3' data.

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

[`hexify`](https://gillescolling.com/hexify/reference/hexify.md) for
assigning points to cells,
[`HexGridInfo-class`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
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

# Mix in aperture 7, either as a family or level by level
grid47 <- hex_grid(area_km2 = 1000, aperture = "4/7")
grid_seq <- hex_grid(resolution = 4, aperture = c(4, 4, 7, 3))

# Grid on another body, by name or by radius
mars <- hex_grid(area_km2 = 1000, radius_km = "mars")
titan <- hex_grid(resolution = 6, radius_km = 2574.76)

# Use grid in hexify
df <- data.frame(lon = c(0, 10, 20), lat = c(45, 50, 55))
result <- hexify(df, lon = "lon", lat = "lat", grid = grid)
```
