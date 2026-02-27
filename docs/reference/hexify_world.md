# Simplified World Map

A lightweight sf object containing simplified world country borders,
suitable for use as a basemap when visualizing hexagonal grids.

## Usage

``` r
hexify_world
```

## Format

An sf object with 177 features and 15 fields:

- name:

  Country short name

- name_long:

  Country full name

- admin:

  Administrative name

- sovereignt:

  Sovereignty

- iso_a2:

  ISO 3166-1 alpha-2 country code

- iso_a3:

  ISO 3166-1 alpha-3 country code

- iso_n3:

  ISO 3166-1 numeric code

- continent:

  Continent name

- region_un:

  UN region

- subregion:

  UN subregion

- region_wb:

  World Bank region

- pop_est:

  Population estimate

- gdp_md:

  GDP in millions USD

- income_grp:

  Income group classification

- economy:

  Economy type

- geometry:

  MULTIPOLYGON geometry in 'WGS84' (EPSG:4326)

## Source

Simplified from Natural Earth 1:110m Cultural Vectors
(<https://www.naturalearthdata.com/>)

## Examples

``` r
library(sf)

# Plot the built-in world map
plot(st_geometry(hexify_world), col = "lightgray", border = "white")

# Filter by continent
europe <- hexify_world[hexify_world$continent == "Europe", ]
plot(st_geometry(europe))
```
