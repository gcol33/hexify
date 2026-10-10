# Validation and benchmarks for the SoftwareX paper

Every number and figure in the hexify SoftwareX article comes from a script in
this directory. Each script writes a CSV to `results/` and a `.meta.txt` beside
it with the date, the hexify commit and installed version, the R and package
versions and the processor.

| Script | What it measures | Output |
|---|---|---|
| `bench_dggrid_agreement.R` | Cell assignment, centres and corners against DGGRID for apertures 3, 4, 7, ISEA43H and a mixed sequence, 200,000 random points per grid; a cell whose corners differ by more than 1 m is checked against its neighbours in both programs. The argument `fuller` runs the same comparison on DGGRID's FULLER grids | `dggrid_agreement.csv`, `dggrid_disagreements.csv`; with `fuller`, `dggrid_agreement_fuller.csv`, `dggrid_disagreements_fuller.csv` |
| `bench_dggal_agreement.R` | Every zone of DGGAL's IVEA3H (levels 1 to 8) and IVEA7H (1 to 5), with ISEA3H and ISEA7H as the control, against hexify's IVEA and ISEA grids on DGGAL's orientation: whether zones map one-to-one onto cells, centre distance, and distance from each DGGAL vertex to the nearest corner of the same hexify cell. `dggal_zones.py` lists the zones through DGGAL's Python package | `dggal_agreement.csv` |
| `make_dggal_fixture.R` | DGGAL's IVEA3H zones at level 3 and a sample of IVEA7H at level 3, centroids and vertices, for the test suite | `tests/testthat/data/dggal_ivea.csv` |
| `bench_dggrid_orientation.R` | The same comparison under three given icosahedron orientations and three DGGRID `REGION_CENTER` placements, 50,000 random points per grid | `dggrid_orientation.csv` |
| `bench_dggrid_z7.R` | Aperture-7 Z7 strings against DGGRID's (IGEO7) for every cell at resolutions 0 to 7, and 0 to 5 under a given orientation and on FULLER; also each program's round trip and the parent named by DGGRID's string | `dggrid_z7.csv` |
| `make_dggrid_z7_int_fixture.R` | DGGRID's IGEO7 Z7 indices as digit strings and as packed 64-bit integers (INT64, hexadecimal), every cell at resolutions 0 to 3 and the cells of 200 random points at 5, 10, 15 and 20, for the test suite | `tests/testthat/data/dggrid_z7_int.csv` |
| `bench_z7_neighbors.R` | Aperture-7 neighbours by IGEO7 digit arithmetic against the quad lattice with quad-leaving steps sent through lon/lat: time per cell for 20,000 cells along quad edges and 20,000 cells of random points at even resolutions 2 to 20, and every cell at resolutions 5 and 6, with the count of cells on whose neighbours the two agree | `z7_neighbors.csv` |
| `make_dggrid_fixture.R` | DGGRID cells and centres under non-standard orientations and on the FULLER projection, and DGGRID's Z7 strings, for the test suite | `tests/testthat/data/dggrid_reference.csv`, `tests/testthat/data/dggrid_z7.csv` |
| `bench_h3_agreement.R` | H3 cells and centres against the H3 library (h3r), resolutions 0 to 15 | `h3_agreement.csv` |
| `bench_cell_area.R` | Cell area against latitude, measured with s2, for ISEA, H3 and a 1-degree grid of about 12,400 km2 per cell | `cell_area_by_latitude.csv`, `cell_area_summary.csv` |
| `bench_grid_metrics.R` | Compactness, intercell distance and cell wall midpoint ratio of every cell and wall of ISEA and FULLER grids of apertures 3, 4 and 7 and of H3, with the spread of cell areas | `grid_metrics.csv` |
| `bench_dggrid_walls.R` | Cell wall midpoint ratio of hexify's and DGGRID's aperture-7 cells at resolutions 3 and 5, both read as corners joined by great-circle arcs; counts the walls whose ratio differs between the two programs and how many of them border a cell whose DGGRID corners differ from hexify's | `dggrid_wall_ratio.csv` |
| `bench_speed.R` | Run time of point assignment (one million points) and polygon generation (10,000 cells) for hexify, dggridR, h3r and h3o | `speed.csv` |
| `example_world_cities.R` | The article's example: population of `maps::world.cities` per equal-area cell | `example_world_cities.csv`, `.gpkg` |
| `bench_kmoch_distortion.R` | Normalized area and isoperimetric quotient of every cell of ISEA and FULLER grids of apertures 3, 4, 7 and 4/3 and of H3, through the pipeline of Kmoch et al. (2022): corner polygons, antimeridian-crossing cells removed, each cell in a Lambert azimuthal equal-area plane on WGS84; beside their published tables, and cell by cell against their per-cell files (downloaded from Zenodo record 6634479 to `data/kmoch2022/`) | `kmoch_distortion.csv`, `kmoch_percell.csv` |
| `bench_hierarchy_commutation.R` | How often a point's cell taken k levels up differs from the cell of the point at that resolution (Griffin 2026, Sec. 12), 10,000 uniform points, k = 1 to 6, against the lattice's one-level rates 4/9, 3/8 and 1/14 | `hierarchy_commutation.csv` |
| `bench_parent_shares.R` | How each cell divides among the coarser cells it overlaps, one level up: every share of every cell measured on the sphere as `hex_aggregate()` measures Fuller and H3 shares (cells clipped where straight, pieces measured on the sphere), against the lattice's 1/3, 1/2, 11/12 and 1/12 and against an independent s2 intersection, on ISEA grids of apertures 3, 4, 7, 4/3 and a per-level sequence, IVEA7H, an octahedral grid, FULLER3H/4H/7H and H3; and a uniform density aggregated by `hex_aggregate()` under both rules against each parent's area | `parent_shares.csv` |
| `bench_ogc_conformance.R` | Distance between cell centre and spherical area centroid (OGC Topic 21 Req 27), and each cell's area on WGS84 relative to the sphere of equal area (Req 28) | `centroid_offset.csv`, `wgs84_area_scale.csv` |
| `bench_trunc_icosa.R` | The truncated icosahedron carried into the icosahedron's face frame by slicing Snyder's eighteen spherical right triangles of a face onto plane triangles of proportional area: area shares, where its vertices land on the frame, the offset along the hexagon-pentagon edges, and the angular deformation over a face (100,000 points) and within 0.5 degrees of a vertex, beside ISEA, IVEA and FULLER (2,000,000 points on the sphere) and the bound 2 asin(1/11) for equal-area projections near a vertex | `trunc_icosa.csv` |
| `bench_globe_table.R` | The globe widget's cells' table against the per-cell upload of a base commit (default `7996043`, installed beside the working tree): R build time (median of five), page size, bytes of the cells on the graphics card, and the graphics card's time per frame (40 frames per submission, mean of the middle ten of twenty) at 800 and 1600 pixels in headless Chrome, for values on every cell, on a sample of cells of a fine grid, and none | `globe_table.csv` |
| `bench_dggal_zirs.R` | Every zone of DGGAL's ISEA3H and IVEA3H (levels 0 to 10) and ISEA7H and IVEA7H (0 to 5) against hexify's OGC grids (`ellipsoid = "WGS84"`, `orientation = "ogc"`): one-to-one mapping onto cells, distance between DGGAL's WGS84 centroid and the cell centre, and how many textZIRS and uint64ZIRS identifiers `cell_to_index()` writes as DGGAL does and `index_to_cell()` reads back. `dggal_zone_ids.py` lists the zones | `dggal_zirs.csv` |
| `make_dggal_zirs_fixture.R` | A sample of DGGAL's zones with their identifiers and WGS84 centroids, for the test suite | `tests/testthat/data/dggal_zirs.csv` |
| `bench_ellipsoid_area.R` | Geodesic area on WGS84 (GeographicLib through geosphere) of every cell of ISEA3H, ISEA4H and ISEA7H grids, over the area `cell_area()` reports, for the default grid and for `ellipsoid = "WGS84"`, overall and by 10-degree latitude band | `ellipsoid_area.csv`, `ellipsoid_area_by_latitude.csv` |
| `bench_lambert_nudge.R` | Lambert's construction of Snyder's projection (`projection_stages()`) over one face of the icosahedron and of the octahedron, 20,301 points of a triangular lattice on the plane triangle: the swing (the distance from the sphere point to the Lambert point) over the face edge's chord, the nudge (Lambert point to nudged point) over the nudged triangle's edge, the distance from the Lambert point to Snyder's in the plane and in space over the plane triangle's edge, and the turn of the azimuth, largest and mean, with where the largest lies | `lambert_nudge.csv` |
| `fig_lambert_stages.R` | Still frames of one face of the aperture-3 resolution-3 grid at the Lambert step (curved, equal-area) and at Snyder's (straight), with cell walls and Tissot's ellipses, at one scale, text-free | `fig-lambert-curved.png`, `fig-lambert-straight.png` |

## Requirements

- hexify installed from this repository (`R CMD INSTALL --preclean .`).
- dggridR, h3r, h3o, maps and sf from CRAN; geosphere for `bench_ellipsoid_area.R`.
- For the DGGAL scripts, a Python (3.6 to 3.13) with `pip install dggal`;
  point `DGGAL_PYTHON` at it.
- DGGRID built from source (<https://github.com/sahrk/DGGRID>), because
  dggridR supports apertures 3 and 4 only. Configure with
  `cmake -S . -B build -DWITH_GDAL=OFF` and build; point `DGGRID_EXE` at
  `build/src/apps/dggrid/dggrid`.

Run each script from the repository root with `Rscript paper/bench/<script>.R`.
Run `bench_speed.R` alone on an otherwise idle machine.
