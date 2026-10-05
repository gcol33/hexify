# Validation and benchmarks for the SoftwareX paper

Every number and figure in the hexify SoftwareX article comes from a script in
this directory. Each script writes a CSV to `results/` and a `.meta.txt` beside
it with the date, the hexify commit and installed version, the R and package
versions and the processor.

| Script | What it measures | Output |
|---|---|---|
| `bench_dggrid_agreement.R` | Cell assignment, centres and corners against DGGRID for apertures 3, 4, 7, ISEA43H and a mixed sequence, 200,000 random points per grid; a cell whose corners differ by more than 1 m is checked against its neighbours in both programs | `dggrid_agreement.csv`, `dggrid_disagreements.csv` |
| `bench_dggrid_orientation.R` | The same comparison under three given icosahedron orientations and three DGGRID `REGION_CENTER` placements, 50,000 random points per grid | `dggrid_orientation.csv` |
| `make_dggrid_orientation_fixture.R` | DGGRID cells and centres under non-standard orientations, for the test suite | `tests/testthat/data/dggrid_orientation.csv` |
| `bench_h3_agreement.R` | H3 cells and centres against the H3 library (h3r), resolutions 0 to 15 | `h3_agreement.csv` |
| `bench_cell_area.R` | Cell area against latitude, measured with s2, for ISEA, H3 and a 1-degree grid of about 12,400 km2 per cell | `cell_area_by_latitude.csv`, `cell_area_summary.csv` |
| `bench_speed.R` | Run time of point assignment (one million points) and polygon generation (10,000 cells) for hexify, dggridR, h3r and h3o | `speed.csv` |
| `example_world_cities.R` | The article's example: population of `maps::world.cities` per equal-area cell | `example_world_cities.csv`, `.gpkg` |

## Requirements

- hexify installed from this repository (`R CMD INSTALL --preclean .`).
- dggridR, h3r, h3o, maps and sf from CRAN.
- DGGRID built from source (<https://github.com/sahrk/DGGRID>), because
  dggridR supports apertures 3 and 4 only. Configure with
  `cmake -S . -B build -DWITH_GDAL=OFF` and build; point `DGGRID_EXE` at
  `build/src/apps/dggrid/dggrid`.

Run each script from the repository root with `Rscript paper/bench/<script>.R`.
Run `bench_speed.R` alone on an otherwise idle machine.
