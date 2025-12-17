# Package index

## Core Classes

S4 classes for grid specification and hexified data

- [`` `$`( ``*`<HexGridInfo>`*`)`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
  [`names(`*`<HexGridInfo>`*`)`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
  [`show(`*`<HexGridInfo>`*`)`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
  [`as.list(`*`<HexGridInfo>`*`)`](https://gcol33.github.io/hexify/reference/HexGridInfo-class.md)
  : HexGridInfo Class
- [`grid_info(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`cells(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`n_cells(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`nrow(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`ncol(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`dim(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`names(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`` `$`( ``*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`` `$<-`( ``*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`` `[`( ``*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`` `[[`( ``*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`` `[[<-`( ``*`<HexData>`*`,`*`<ANY>`*`,`*`<missing>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`show(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`as.data.frame(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  [`as.list(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
  : HexData Class
- [`hex_grid()`](https://gcol33.github.io/hexify/reference/hex_grid.md)
  : Create a Hexagonal Grid Specification
- [`hexify()`](https://gcol33.github.io/hexify/reference/hexify.md) :
  Assign hexagonal DGGS cell IDs to geographic points
- [`hexify_assign()`](https://gcol33.github.io/hexify/reference/hexify_assign.md)
  : Assign hex cells (ISEA3H, aperture 3) for lon/lat

## HexData Methods

Methods for working with HexData objects

- [`grid_info()`](https://gcol33.github.io/hexify/reference/grid_info.md)
  : Get Grid Specification
- [`cells()`](https://gcol33.github.io/hexify/reference/cells.md) : Get
  Cell IDs
- [`n_cells()`](https://gcol33.github.io/hexify/reference/n_cells.md) :
  Get Number of Cells
- [`as_sf()`](https://gcol33.github.io/hexify/reference/as_sf.md) :
  Convert HexData to sf Object
- [`as_tibble.HexData()`](https://gcol33.github.io/hexify/reference/as_tibble.HexData.md)
  : Convert HexData to tibble
- [`is_hex_grid()`](https://gcol33.github.io/hexify/reference/is_hex_grid.md)
  : Check if object is HexGridInfo
- [`is_hex_data()`](https://gcol33.github.io/hexify/reference/is_hex_data.md)
  : Check if object is HexData

## Grid Generation

Functions for creating hexagonal grids over regions

- [`grid_rect()`](https://gcol33.github.io/hexify/reference/grid_rect.md)
  : Generate a rectangular grid of hexagons
- [`grid_global()`](https://gcol33.github.io/hexify/reference/grid_global.md)
  : Generate a global hexagon grid
- [`grid_clip()`](https://gcol33.github.io/hexify/reference/grid_clip.md)
  : Clip hexagon grid to polygon boundary
- [`cell_to_sf()`](https://gcol33.github.io/hexify/reference/cell_to_sf.md)
  : Convert cell IDs to sf polygons
- [`lonlat_to_cell()`](https://gcol33.github.io/hexify/reference/lonlat_to_cell.md)
  : Convert longitude/latitude to cell ID
- [`cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude

## Visualization

Plotting and mapping functions

- [`hexify_ggplot()`](https://gcol33.github.io/hexify/reference/hexify_ggplot.md)
  : Create a ggplot2 visualization of HexData
- [`hexify_heatmap()`](https://gcol33.github.io/hexify/reference/hexify_heatmap.md)
  : Create a heatmap visualization of hexagonal grid cells
- [`plot_world()`](https://gcol33.github.io/hexify/reference/plot_world.md)
  : Quick world map plot
- [`hexify_world`](https://gcol33.github.io/hexify/reference/hexify_world.md)
  : Simplified World Map

## Grid Statistics

Resolution comparison and grid statistics

- [`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md)
  : Get grid statistics for Earth coverage
- [`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md)
  : Compare grid resolutions
- [`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md)
  **\[deprecated\]** : Print resolution comparison table
- [`dgverify()`](https://gcol33.github.io/hexify/reference/dgverify.md)
  : Verify grid object

## dggridR Compatibility

Functions for interoperability with dggridR package

- [`as_dggrid()`](https://gcol33.github.io/hexify/reference/as_dggrid.md)
  : Convert hexify grid to dggridR-compatible grid object
- [`from_dggrid()`](https://gcol33.github.io/hexify/reference/from_dggrid.md)
  : Convert dggridR grid object to hexify_grid
- [`dggrid_is_compatible()`](https://gcol33.github.io/hexify/reference/dggrid_is_compatible.md)
  : Validate dggridR grid compatibility with hexify
- [`dggrid_43h_sequence()`](https://gcol33.github.io/hexify/reference/dggrid_43h_sequence.md)
  : Create DGGRID 43H aperture sequence
- [`hexify_grid()`](https://gcol33.github.io/hexify/reference/hexify_grid.md)
  : Create a hexagonal grid specification

## Low-Level Coordinate Conversions

Internal coordinate transformation functions

- [`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md)
  : Forward Snyder projection
- [`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md)
  : Forward projection to specific face
- [`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md)
  : Inverse Snyder projection
- [`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)
  : Determine which face contains a point

## Low-Level Grid Functions

Internal grid construction and projection functions

- [`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md)
  : Initialize icosahedron geometry
- [`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md)
  : Get icosahedron face centers
- [`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md)
  : Generate a rectangular grid of hexagon polygons
- [`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md)
  : Generate a global grid of hexagon polygons
- [`hexify_grid_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md)
  : Convert longitude/latitude to cell ID using a grid object
- [`hexify_grid_cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude using a grid object
- [`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md)
  : Build an sf POLYGON from six (lon, lat) corner pairs
- [`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md)
  : Get inverse projection statistics

## Resolution & Index Utilities

Internal resolution and index functions

- [`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md)
  : Find closest resolution for target cell area
- [`hexify_area_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_area_to_eff_res.md)
  : Convert area to effective resolution
- [`hexify_eff_res_to_area()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_area.md)
  : Convert effective resolution to area
- [`hexify_eff_res_to_resolution()`](https://gcol33.github.io/hexify/reference/hexify_eff_res_to_resolution.md)
  : Convert effective resolution to index resolution
- [`hexify_resolution_to_eff_res()`](https://gcol33.github.io/hexify/reference/hexify_resolution_to_eff_res.md)
  : Convert index resolution to effective resolution
- [`hexify_compare_indices()`](https://gcol33.github.io/hexify/reference/hexify_compare_indices.md)
  : Compare two indices
- [`hexify_default_index_type()`](https://gcol33.github.io/hexify/reference/hexify_default_index_type.md)
  : Get default index type for aperture
- [`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md)
  : Check if index type is valid for aperture

## Hierarchy Functions

Parent-child cell navigation

- [`get_parent()`](https://gcol33.github.io/hexify/reference/get_parent.md)
  : Get parent cell
- [`get_children()`](https://gcol33.github.io/hexify/reference/get_children.md)
  : Get children cells
- [`cell_to_index()`](https://gcol33.github.io/hexify/reference/cell_to_index.md)
  : Convert cell ID to hierarchical index string
- [`hexify_get_children()`](https://gcol33.github.io/hexify/reference/hexify_get_children.md)
  : Get children indices
- [`hexify_get_parent()`](https://gcol33.github.io/hexify/reference/hexify_get_parent.md)
  : Get parent index
- [`hexify_get_precision()`](https://gcol33.github.io/hexify/reference/hexify_get_precision.md)
  : Get current precision settings
- [`hexify_get_resolution()`](https://gcol33.github.io/hexify/reference/hexify_get_resolution.md)
  : Get index resolution
- [`hexify_set_precision()`](https://gcol33.github.io/hexify/reference/hexify_set_precision.md)
  : Set inverse projection precision
- [`hexify_z7_canonical()`](https://gcol33.github.io/hexify/reference/hexify_z7_canonical.md)
  : Get canonical form of Z7 index

## Utilities

Testing and debugging functions

- [`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)
  : Round-trip accuracy test
- [`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md)
  : Set verbose mode for inverse projection

## Internal

Internal documentation topics

- [`hexify-package`](https://gcol33.github.io/hexify/reference/hexify-package.md)
  : hexify
- [`hexify-conversions`](https://gcol33.github.io/hexify/reference/hexify-conversions.md)
  : Coordinate Conversions
- [`hexify-grid`](https://gcol33.github.io/hexify/reference/hexify-grid.md)
  : Core Grid Construction
- [`hexify-stats`](https://gcol33.github.io/hexify/reference/hexify-stats.md)
  : Grid Statistics
- [`plot(`*`<HexData>`*`,`*`<missing>`*`)`](https://gcol33.github.io/hexify/reference/plot-HexData-missing-method.md)
  : Plot HexData objects
