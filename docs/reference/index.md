# Package index

## Core

Create grids and assign points to cells

- [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md)
  : Create a Hexagonal Grid Specification
- [`hexify()`](https://gillescolling.com/hexify/reference/hexify.md) :
  Assign hexagonal DGGS cell IDs to geographic points
- [`hexify_assign()`](https://gillescolling.com/hexify/reference/hexify_assign.md)
  : Assign hex cells ('ISEA3H', aperture 3) for lon/lat
- [`HexGridInfo-class`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  : HexGridInfo Class
- [`HexData-class`](https://gillescolling.com/hexify/reference/HexData-class.md)
  : HexData Class

## Working with Results

Extract and convert HexData

- [`grid_info()`](https://gillescolling.com/hexify/reference/grid_info.md)
  : Get Grid Specification
- [`st_crs(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/st_crs.HexGridInfo.md)
  [`st_crs(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/st_crs.HexGridInfo.md)
  : Coordinate reference system of a grid
- [`cells()`](https://gillescolling.com/hexify/reference/cells.md) : Get
  Cell IDs
- [`n_cells()`](https://gillescolling.com/hexify/reference/n_cells.md) :
  Get Number of Cells
- [`summary(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/hexify-summary.md)
  [`print(`*`<hexify_grid_summary>`*`)`](https://gillescolling.com/hexify/reference/hexify-summary.md)
  [`summary(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/hexify-summary.md)
  [`print(`*`<hexify_data_summary>`*`)`](https://gillescolling.com/hexify/reference/hexify-summary.md)
  : Summary of a grid or of gridded data
- [`st_as_sf(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/st_as_sf.HexData.md)
  [`st_as_sf(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/st_as_sf.HexData.md)
  : Convert gridded data or a grid to sf
- [`as_tibble(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/as_tibble.HexData.md)
  : Convert HexData to tibble
- [`is_hex_grid()`](https://gillescolling.com/hexify/reference/is_hex_grid.md)
  : Check if object is HexGridInfo
- [`is_hex_data()`](https://gillescolling.com/hexify/reference/is_hex_data.md)
  : Check if object is HexData
- [`hex_summarize()`](https://gillescolling.com/hexify/reference/hex_summarize.md)
  : Summarize Data by Hex Cell
- [`hex_smooth()`](https://gillescolling.com/hexify/reference/hex_smooth.md)
  : Smooth cell values over neighbouring cells
- [`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md)
  : Zonal Statistics for Hex Cells
- [`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md)
  : Extract Raster Values at Hex Cell Centers

## Grid Generation

Create grids over regions

- [`grid_rect()`](https://gillescolling.com/hexify/reference/grid_rect.md)
  : Generate a rectangular grid of hexagons
- [`grid_global()`](https://gillescolling.com/hexify/reference/grid_global.md)
  : Generate a global hexagon grid
- [`grid_clip()`](https://gillescolling.com/hexify/reference/grid_clip.md)
  : Clip hexagon grid to polygon boundary
- [`cell_to_sf()`](https://gillescolling.com/hexify/reference/cell_to_sf.md)
  : Convert cell IDs to sf polygons
- [`lonlat_to_cell()`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md)
  : Convert longitude/latitude to cell ID
- [`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude
- [`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md)
  : Get parent cell
- [`get_children()`](https://gillescolling.com/hexify/reference/get_children.md)
  : Get children cells
- [`get_siblings()`](https://gillescolling.com/hexify/reference/get_siblings.md)
  : Sibling cells

## Cell Operations

Relations between cells and cell sets

- [`get_neighbors()`](https://gillescolling.com/hexify/reference/get_neighbors.md)
  : Get Neighboring Cells
- [`hex_distance()`](https://gillescolling.com/hexify/reference/hex_distance.md)
  : Grid Distance Between Cells
- [`is_pentagon()`](https://gillescolling.com/hexify/reference/is_pentagon.md)
  : Detect Pentagon Cells
- [`hex_compact()`](https://gillescolling.com/hexify/reference/hex_compact.md)
  : Compact Hex Cells
- [`hex_uncompact()`](https://gillescolling.com/hexify/reference/hex_uncompact.md)
  : Uncompact Hex Cells
- [`hex_sample()`](https://gillescolling.com/hexify/reference/hex_sample.md)
  : Random points inside cells

## H3 Interoperability

Cross-grid mapping and per-cell area

- [`h3_crosswalk()`](https://gillescolling.com/hexify/reference/h3_crosswalk.md)
  : Crosswalk Between ISEA and H3 Cell IDs
- [`cell_area()`](https://gillescolling.com/hexify/reference/cell_area.md)
  : Compute per-cell area in km²
- [`import_h3()`](https://gillescolling.com/hexify/reference/import_h3.md)
  : Import External H3 Cell IDs into hexify
- [`dggrs_definition()`](https://gillescolling.com/hexify/reference/dggrs_definition.md)
  : Definition of a grid as a Discrete Global Grid Reference System

## Cell Metrics

Shape and spacing of a grid’s cells

- [`cell_metrics()`](https://gillescolling.com/hexify/reference/cell_metrics.md)
  : Shape metrics of cells: area, perimeter and compactness
- [`wall_metrics()`](https://gillescolling.com/hexify/reference/wall_metrics.md)
  : Spacing metrics of neighbouring cells: wall length, centre distance
  and cell wall midpoint ratio
- [`projection_distortion()`](https://gillescolling.com/hexify/reference/projection_distortion.md)
  : Distortion of a grid's face projection

## Visualization

Plotting functions

- [`hexify_heatmap()`](https://gillescolling.com/hexify/reference/hexify_heatmap.md)
  : Create a ggplot2 visualization of hexagonal grid cells
- [`plot_grid()`](https://gillescolling.com/hexify/reference/plot_grid.md)
  : Plot hexagonal grid clipped to a polygon boundary
- [`plot_world()`](https://gillescolling.com/hexify/reference/plot_world.md)
  : Quick world map plot
- [`plot(`*`<HexGridInfo>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/plot-HexGridInfo-missing-method.md)
  : Plot a grid on the sphere, on its solid, or on the unfolded net
- [`net_layout()`](https://gillescolling.com/hexify/reference/net_layout.md)
  : Lay out the faces of a grid's solid in the plane
- [`net_project()`](https://gillescolling.com/hexify/reference/net_project.md)
  : Place points on a net layout
- [`net_cells()`](https://gillescolling.com/hexify/reference/net_cells.md)
  : The parts of cells on a net layout
- [`hex_rays()`](https://gillescolling.com/hexify/reference/hex_rays.md)
  : Ray glyphs for cell values
- [`hex_triangles()`](https://gillescolling.com/hexify/reference/hex_triangles.md)
  : The triangles of cells towards each neighbour
- [`hex_globe()`](https://gillescolling.com/hexify/reference/hex_globe.md)
  : Interactive globe of a grid
- [`hex_globe_png()`](https://gillescolling.com/hexify/reference/hex_globe_png.md)
  : Save a globe as a PNG image
- [`hexify_world`](https://gillescolling.com/hexify/reference/hexify_world.md)
  : Simplified World Map
- [`globe_centers`](https://gillescolling.com/hexify/reference/globe_centers.md)
  : Globe center presets
- [`hex_browse()`](https://gillescolling.com/hexify/reference/hex_browse.md)
  : Interactive Hex Map
- [`plot(`*`<HexData>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/plot-HexData-missing-method.md)
  : Plot HexData objects

## dggridR Compatibility

Interoperability with dggridR

- [`as_dggrid()`](https://gillescolling.com/hexify/reference/as_dggrid.md)
  : Convert hexify grid to 'dggridR'-compatible grid object
- [`from_dggrid()`](https://gillescolling.com/hexify/reference/from_dggrid.md)
  : Convert 'dggridR' grid object to hexify_grid
- [`dggrid_is_compatible()`](https://gillescolling.com/hexify/reference/dggrid_is_compatible.md)
  : Validate 'dggridR' grid compatibility with hexify
- [`dggrid_43h_sequence()`](https://gillescolling.com/hexify/reference/dggrid_43h_sequence.md)
  : Create DGGRID 43H aperture sequence
- [`hexify_grid()`](https://gillescolling.com/hexify/reference/hexify_grid.md)
  : Create a hexagonal grid specification
- [`dgearthstat()`](https://gillescolling.com/hexify/reference/dgearthstat.md)
  : Get grid statistics for whole-body coverage
- [`dgverify()`](https://gillescolling.com/hexify/reference/dgverify.md)
  : Verify grid object
- [`hexify_compare_resolutions()`](https://gillescolling.com/hexify/reference/hexify_compare_resolutions.md)
  : Compare grid resolutions

## Low-Level API

Direct access to projection and coordinate transforms

- [`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md)
  : Forward face projection
- [`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md)
  : Forward projection to specific face
- [`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md)
  : Inverse face projection
- [`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)
  : Determine which face contains a point
- [`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md)
  : Set the default icosahedron orientation
- [`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md)
  : Get face centers of a solid
- [`hexify_grid_rect()`](https://gillescolling.com/hexify/reference/hexify_grid_rect.md)
  : Generate a rectangular grid of hexagon polygons
- [`hexify_grid_global()`](https://gillescolling.com/hexify/reference/hexify_grid_global.md)
  : Generate a global grid of hexagon polygons
- [`hexify_grid_to_cell()`](https://gillescolling.com/hexify/reference/hexify_grid_to_cell.md)
  : Convert longitude/latitude to cell ID using a grid object
- [`hexify_grid_cell_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_grid_cell_to_lonlat.md)
  : Convert cell ID to longitude/latitude using a grid object
- [`hex_corners_to_sf()`](https://gillescolling.com/hexify/reference/hex_corners_to_sf.md)
  : Build an sf POLYGON from six (lon, lat) corner pairs
- [`dg_closest_res_to_area()`](https://gillescolling.com/hexify/reference/dg_closest_res_to_area.md)
  : Find closest resolution for target cell area
- [`hexify_area_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_area_to_eff_res.md)
  : Convert area to effective resolution
- [`hexify_eff_res_to_area()`](https://gillescolling.com/hexify/reference/hexify_eff_res_to_area.md)
  : Convert effective resolution to area
- [`hexify_eff_res_to_resolution()`](https://gillescolling.com/hexify/reference/hexify_eff_res_to_resolution.md)
  : Convert effective resolution to index resolution
- [`hexify_resolution_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_resolution_to_eff_res.md)
  : Convert index resolution to effective resolution
- [`hexify_compare_indices()`](https://gillescolling.com/hexify/reference/hexify_compare_indices.md)
  : Compare two indices
- [`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md)
  : Get default index type for aperture
- [`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md)
  : Check if index type is valid for aperture
- [`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)
  : Convert cell ID to hierarchical index string
- [`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md)
  : Get children indices
- [`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md)
  : Get parent index
- [`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md)
  : Get index resolution
- [`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)
  : Get canonical form of Z7 index
- [`hexify_roundtrip_test()`](https://gillescolling.com/hexify/reference/hexify_roundtrip_test.md)
  : Round-trip accuracy test

## Internal

Package internals

- [`hexify-package`](https://gillescolling.com/hexify/reference/hexify-package.md)
  : hexify
- [`hexify-conversions`](https://gillescolling.com/hexify/reference/hexify-conversions.md)
  : Coordinate Conversions
- [`hexify-grid`](https://gillescolling.com/hexify/reference/hexify-grid.md)
  : Core Grid Construction
- [`hexify-stats`](https://gillescolling.com/hexify/reference/hexify-stats.md)
  : Grid Statistics
- [`` `$`( ``*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-methods.md)
  [`n_cells(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-methods.md)
  [`names(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-methods.md)
  [`show(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-methods.md)
  [`as.list(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-methods.md)
  : HexGridInfo S4 Methods
- [`grid_info(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`cells(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`n_cells(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`nrow(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`ncol(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`dim(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`names(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`` `$`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`` `$<-`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`` `[`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`` `[[`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`` `[[<-`( ``*`<HexData>`*`,`*`<ANY>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`show(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`as.data.frame(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  [`as.list(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-methods.md)
  : HexData S4 Methods
- [`new_hex_data()`](https://gillescolling.com/hexify/reference/new_hex_data.md)
  : Create a HexData Object (Internal)
- [`extract_grid()`](https://gillescolling.com/hexify/reference/extract_grid.md)
  : Extract grid from various objects
- [`calculate_resolution_for_area()`](https://gillescolling.com/hexify/reference/calculate_resolution_for_area.md)
  : Calculate resolution for target area
- [`hexify_grid_to_HexGridInfo()`](https://gillescolling.com/hexify/reference/hexify_grid_to_HexGridInfo.md)
  : Convert legacy hexify_grid to HexGridInfo
- [`HexGridInfo_to_hexify_grid()`](https://gillescolling.com/hexify/reference/HexGridInfo_to_hexify_grid.md)
  : Convert HexGridInfo to legacy hexify_grid
