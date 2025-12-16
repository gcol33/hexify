# Package index

## Core Classes

- [`` `$`( ``*`<HexGrid>`*`)`](https://gcol33.github.io/hexify/reference/HexGrid-class.md)
  [`names(`*`<HexGrid>`*`)`](https://gcol33.github.io/hexify/reference/HexGrid-class.md)
  [`show(`*`<HexGrid>`*`)`](https://gcol33.github.io/hexify/reference/HexGrid-class.md)
  [`as.list(`*`<HexGrid>`*`)`](https://gcol33.github.io/hexify/reference/HexGrid-class.md)
  : HexGrid Class
- [`grid(`*`<HexData>`*`)`](https://gcol33.github.io/hexify/reference/HexData-class.md)
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

## Class Methods

- [`grid()`](https://gcol33.github.io/hexify/reference/grid.md) : Get
  Grid Specification
- [`cells()`](https://gcol33.github.io/hexify/reference/cells.md) : Get
  Cell IDs
- [`n_cells()`](https://gcol33.github.io/hexify/reference/n_cells.md) :
  Get Number of Cells
- [`as_sf()`](https://gcol33.github.io/hexify/reference/as_sf.md) :
  Convert HexData to sf Object
- [`as_tibble.HexData()`](https://gcol33.github.io/hexify/reference/as_tibble.HexData.md)
  : Convert HexData to tibble
- [`autoplot.HexData()`](https://gcol33.github.io/hexify/reference/autoplot.HexData.md)
  : Create a ggplot2 visualization of HexData
- [`plot(`*`<HexData>`*`,`*`<missing>`*`)`](https://gcol33.github.io/hexify/reference/plot-HexData-missing-method.md)
  : Plot HexData objects

## Grid Helpers

- [`lonlat_to_cell()`](https://gcol33.github.io/hexify/reference/lonlat_to_cell.md)
  : Convert longitude/latitude to cell ID
- [`cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude
- [`cell_to_sf()`](https://gcol33.github.io/hexify/reference/cell_to_sf.md)
  : Convert cell IDs to sf polygons
- [`cell_to_index()`](https://gcol33.github.io/hexify/reference/cell_to_index.md)
  : Convert cell ID to hierarchical index string
- [`grid_rect()`](https://gcol33.github.io/hexify/reference/grid_rect.md)
  : Generate a rectangular grid of hexagons
- [`grid_global()`](https://gcol33.github.io/hexify/reference/grid_global.md)
  : Generate a global hexagon grid
- [`get_parent()`](https://gcol33.github.io/hexify/reference/get_parent.md)
  : Get parent cell
- [`get_children()`](https://gcol33.github.io/hexify/reference/get_children.md)
  : Get children cells
- [`extract_grid()`](https://gcol33.github.io/hexify/reference/extract_grid.md)
  : Extract grid from various objects
- [`is_hex_grid()`](https://gcol33.github.io/hexify/reference/is_hex_grid.md)
  : Check if object is HexGrid
- [`is_hex_data()`](https://gcol33.github.io/hexify/reference/is_hex_data.md)
  : Check if object is HexData
- [`hexify_df()`](https://gcol33.github.io/hexify/reference/hexify_df.md)
  : Extract plain data frame from hexify result
- [`new_hex_data()`](https://gcol33.github.io/hexify/reference/new_hex_data.md)
  : Create a HexData Object (Internal)

## Coordinate Conversions

- [`hexify_lonlat_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_cell.md)
  : Convert longitude/latitude to cell ID
- [`hexify_lonlat_to_index()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_index.md)
  : Convert longitude/latitude to index string
- [`hexify_lonlat_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_plane.md)
  : Convert longitude/latitude to PLANE coordinates
- [`hexify_lonlat_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_lonlat_to_quad_ij.md)
  : Convert longitude/latitude to Quad IJ coordinates
- [`hexify_cell_id_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_cell_id_to_quad_ij.md)
  : Get cell info from cell ID
- [`hexify_cell_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_icosa_tri.md)
  : Convert Cell ID to Icosa Triangle coordinates
- [`hexify_cell_to_index()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_index.md)
  : Convert cell coordinates to index string
- [`hexify_cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude
- [`hexify_cell_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_plane.md)
  : Convert Cell ID to PLANE coordinates
- [`hexify_cell_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_ij.md)
  : Convert Cell ID to Quad IJ coordinates
- [`hexify_cell_to_quad_xy()`](https://gcol33.github.io/hexify/reference/hexify_cell_to_quad_xy.md)
  : Convert Cell ID to Quad XY coordinates
- [`hexify_quad_ij_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_cell.md)
  : Convert Quad IJ coordinates to cell ID
- [`hexify_quad_ij_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_icosa_tri.md)
  : Convert Quad IJ to Icosa Triangle coordinates
- [`hexify_quad_ij_to_xy()`](https://gcol33.github.io/hexify/reference/hexify_quad_ij_to_xy.md)
  : Convert Quad IJ to Quad XY (continuous coordinates)
- [`hexify_quad_xy_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_cell.md)
  : Convert Quad XY coordinates to Cell ID
- [`hexify_quad_xy_to_icosa_tri()`](https://gcol33.github.io/hexify/reference/hexify_quad_xy_to_icosa_tri.md)
  : Convert Quad XY to Icosa Triangle coordinates
- [`hexify_grid_cell_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_grid_cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude using a grid object
- [`hexify_grid_global()`](https://gcol33.github.io/hexify/reference/hexify_grid_global.md)
  : Generate a global grid of hexagon polygons
- [`hexify_grid_rect()`](https://gcol33.github.io/hexify/reference/hexify_grid_rect.md)
  : Generate a rectangular grid of hexagon polygons
- [`hexify_grid_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_grid_to_cell.md)
  : Convert longitude/latitude to cell ID using a grid object
- [`hexify_icosa_tri_to_plane()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_plane.md)
  : Convert Icosa Triangle coordinates to PLANE coordinates
- [`hexify_icosa_tri_to_quad_ij()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_quad_ij.md)
  : Convert Icosa Triangle to Quad IJ coordinates
- [`hexify_icosa_tri_to_quad_xy()`](https://gcol33.github.io/hexify/reference/hexify_icosa_tri_to_quad_xy.md)
  : Convert Icosa Triangle to Quad XY coordinates
- [`hexify_index_to_cell()`](https://gcol33.github.io/hexify/reference/hexify_index_to_cell.md)
  : Convert index string to cell coordinates
- [`hexify_index_to_lonlat()`](https://gcol33.github.io/hexify/reference/hexify_index_to_lonlat.md)
  : Convert index string to longitude/latitude
- [`hexify_forward()`](https://gcol33.github.io/hexify/reference/hexify_forward.md)
  : Forward Snyder projection
- [`hexify_forward_to_face()`](https://gcol33.github.io/hexify/reference/hexify_forward_to_face.md)
  : Forward projection to specific face
- [`hexify_inverse()`](https://gcol33.github.io/hexify/reference/hexify_inverse.md)
  : Inverse Snyder projection
- [`hexify_which_face()`](https://gcol33.github.io/hexify/reference/hexify_which_face.md)
  : Determine which face contains a point

## Grid Construction

- [`hexify_build_icosa()`](https://gcol33.github.io/hexify/reference/hexify_build_icosa.md)
  : Initialize icosahedron geometry
- [`hexify_face_centers()`](https://gcol33.github.io/hexify/reference/hexify_face_centers.md)
  : Get icosahedron face centers
- [`hexify_to_polygons()`](https://gcol33.github.io/hexify/reference/hexify_to_polygons.md)
  : Generate polygons directly from hexify result
- [`hexify_to_sf()`](https://gcol33.github.io/hexify/reference/hexify_to_sf.md)
  : Convert hexify result to sf object
- [`hex_corners_to_sf()`](https://gcol33.github.io/hexify/reference/hex_corners_to_sf.md)
  : Build an sf POLYGON from six (lon, lat) corner pairs

## Resolution & Area

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
- [`hexify_compare_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_compare_resolutions.md)
  : Compare grid resolutions
- [`hexify_default_index_type()`](https://gcol33.github.io/hexify/reference/hexify_default_index_type.md)
  : Get default index type for aperture
- [`hexify_is_valid_index_type()`](https://gcol33.github.io/hexify/reference/hexify_is_valid_index_type.md)
  : Check if index type is valid for aperture
- [`hexify_print_resolutions()`](https://gcol33.github.io/hexify/reference/hexify_print_resolutions.md)
  : Print resolution comparison table
- [`hexify_projection_stats()`](https://gcol33.github.io/hexify/reference/hexify_projection_stats.md)
  : Get inverse projection statistics
- [`calculate_resolution_for_area()`](https://gcol33.github.io/hexify/reference/calculate_resolution_for_area.md)
  : Calculate resolution for target area

## Hierarchy Functions

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

## Grid Statistics

- [`dgearthstat()`](https://gcol33.github.io/hexify/reference/dgearthstat.md)
  : Get grid statistics for Earth coverage
- [`dg_closest_res_to_area()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_area.md)
  : Find closest resolution for target cell area
- [`dg_closest_res_to_cls()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_cls.md)
  : Find closest resolution for target CLS
- [`dg_closest_res_to_spacing()`](https://gcol33.github.io/hexify/reference/dg_closest_res_to_spacing.md)
  : Find closest resolution for target cell spacing
- [`dgverify()`](https://gcol33.github.io/hexify/reference/dgverify.md)
  : Verify grid object

## dggridR Compatibility

- [`dggrid_is_compatible()`](https://gcol33.github.io/hexify/reference/dggrid_is_compatible.md)
  : Validate dggridR grid compatibility with hexify
- [`as_dggrid()`](https://gcol33.github.io/hexify/reference/as_dggrid.md)
  : Convert hexify grid to dggridR-compatible grid object
- [`from_dggrid()`](https://gcol33.github.io/hexify/reference/from_dggrid.md)
  : Convert dggridR grid object to hexify_grid

## Visualization

- [`plot_world()`](https://gcol33.github.io/hexify/reference/plot_world.md)
  : Quick world map plot
- [`hexify_map()`](https://gcol33.github.io/hexify/reference/hexify_map.md)
  : Plot hexagonal grid cells with optional basemap
- [`hexify_plot()`](https://gcol33.github.io/hexify/reference/hexify_plot.md)
  : Quick plot of hexify results
- [`hexify_heatmap()`](https://gcol33.github.io/hexify/reference/hexify_heatmap.md)
  : Create a heatmap visualization of hexagonal grid cells
- [`hexify_world`](https://gcol33.github.io/hexify/reference/hexify_world.md)
  : Simplified World Map

## Utilities

- [`hexify_roundtrip_test()`](https://gcol33.github.io/hexify/reference/hexify_roundtrip_test.md)
  : Round-trip accuracy test
- [`hexify_set_verbose()`](https://gcol33.github.io/hexify/reference/hexify_set_verbose.md)
  : Set verbose mode for inverse projection
- [`index_to_cell_internal()`](https://gcol33.github.io/hexify/reference/index_to_cell_internal.md)
  : Decode a cell index to face, i, j, and resolution

## Internal

- [`hexify-package`](https://gcol33.github.io/hexify/reference/hexify-package.md)
  : hexify
- [`hexify-conversions`](https://gcol33.github.io/hexify/reference/hexify-conversions.md)
  : Coordinate Conversions
- [`hexify-grid`](https://gcol33.github.io/hexify/reference/hexify-grid.md)
  : Core Grid Construction
- [`hexify-stats`](https://gcol33.github.io/hexify/reference/hexify-stats.md)
  : Grid Statistics
- [`print(`*`<hexify_grid>`*`)`](https://gcol33.github.io/hexify/reference/print.hexify_grid.md)
  : Print method for hexify_grid objects
