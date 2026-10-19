# Package index

## Core

Create grids and assign points to cells

- [`hex_grid()`](https://gillescolling.com/hexify/reference/hex_grid.md)
  : Create a Hexagonal Grid Specification
- [`hex9_warp_download()`](https://gillescolling.com/hexify/reference/hex9_warp_download.md)
  : Download Hex9's Warp Field
- [`hexify()`](https://gillescolling.com/hexify/reference/hexify.md) :
  Assign Hexagonal DGGS Cell IDs to Geographic Points
- [`hexify_assign()`](https://gillescolling.com/hexify/reference/hexify_assign.md)
  : Assign Hex Cells ('ISEA3H', Aperture 3) for lon/lat
- [`` `$`( ``*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`n_cells(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`names(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`show(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`summary(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`print(`*`<hexify_grid_summary>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`as.list(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`st_as_sf(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`st_crs(`*`<HexGridInfo>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  [`plot(`*`<HexGridInfo>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/HexGridInfo-class.md)
  : HexGridInfo Class
- [`grid_info(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`cells(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`n_cells(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`nrow(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`ncol(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`dim(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`names(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`` `$`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`` `$<-`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`` `[`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`` `[[`( ``*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`` `[[<-`( ``*`<HexData>`*`,`*`<ANY>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`show(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`summary(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`print(`*`<hexify_data_summary>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`as.data.frame(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`as.list(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`as_tibble(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`st_as_sf(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`st_crs(`*`<HexData>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  [`plot(`*`<HexData>`*`,`*`<missing>`*`)`](https://gillescolling.com/hexify/reference/HexData-class.md)
  : HexData Class

## Working with Results

Extract and convert HexData

- [`grid_info()`](https://gillescolling.com/hexify/reference/grid_info.md)
  : Get Grid Specification
- [`cells()`](https://gillescolling.com/hexify/reference/cells.md) : Get
  Cell IDs
- [`n_cells()`](https://gillescolling.com/hexify/reference/n_cells.md) :
  Get Number of Cells
- [`is_hex_grid()`](https://gillescolling.com/hexify/reference/is_hex_grid.md)
  : Check If an Object Is a HexGridInfo
- [`is_hex_data()`](https://gillescolling.com/hexify/reference/is_hex_data.md)
  : Check If an Object Is a HexData
- [`hex_summarize()`](https://gillescolling.com/hexify/reference/hex_summarize.md)
  : Summarize Data by Hex Cell
- [`hex_smooth()`](https://gillescolling.com/hexify/reference/hex_smooth.md)
  : Smooth Cell Values Over Neighbouring Cells
- [`hex_aggregate()`](https://gillescolling.com/hexify/reference/hex_aggregate.md)
  : Aggregate Cell Values to a Coarser Resolution
- [`hex_zonal()`](https://gillescolling.com/hexify/reference/hex_zonal.md)
  : Zonal Statistics for Hex Cells
- [`hex_extract()`](https://gillescolling.com/hexify/reference/hex_extract.md)
  : Extract Raster Values at Hex Cell Centers

## Grid Generation

Create grids over regions

- [`grid_rect()`](https://gillescolling.com/hexify/reference/grid_rect.md)
  : Generate a Rectangular Grid of Hexagons
- [`grid_global()`](https://gillescolling.com/hexify/reference/grid_global.md)
  : Generate a Global Hexagon Grid
- [`grid_clip()`](https://gillescolling.com/hexify/reference/grid_clip.md)
  : Clip Hexagon Grid to Polygon Boundary
- [`cell_to_sf()`](https://gillescolling.com/hexify/reference/cell_to_sf.md)
  : Convert Cell IDs to sf Polygons
- [`lonlat_to_cell()`](https://gillescolling.com/hexify/reference/lonlat_to_cell.md)
  : Convert Longitude/Latitude to Cell ID
- [`cell_to_lonlat()`](https://gillescolling.com/hexify/reference/cell_to_lonlat.md)
  : Convert Cell ID to Longitude/Latitude
- [`get_parent()`](https://gillescolling.com/hexify/reference/get_parent.md)
  : Get Parent Cell
- [`get_children()`](https://gillescolling.com/hexify/reference/get_children.md)
  : Get Children Cells
- [`get_siblings()`](https://gillescolling.com/hexify/reference/get_siblings.md)
  : Sibling Cells

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
  : Random Points Inside Cells

## H3 Interoperability

Cross-grid mapping and per-cell area

- [`h3_crosswalk()`](https://gillescolling.com/hexify/reference/h3_crosswalk.md)
  : Crosswalk Between ISEA and H3 Cell IDs
- [`cell_area()`](https://gillescolling.com/hexify/reference/cell_area.md)
  : Compute Per-Cell Area in km²
- [`import_h3()`](https://gillescolling.com/hexify/reference/import_h3.md)
  : Import External H3 Cell IDs into hexify
- [`dggrs_definition()`](https://gillescolling.com/hexify/reference/dggrs_definition.md)
  : Definition of a Grid as a Discrete Global Grid Reference System

## Cell Metrics

Shape and spacing of a grid’s cells

- [`cell_metrics()`](https://gillescolling.com/hexify/reference/cell_metrics.md)
  : Shape Metrics of Cells: Area, Perimeter and Compactness
- [`wall_metrics()`](https://gillescolling.com/hexify/reference/wall_metrics.md)
  : Spacing Metrics of Neighbouring Cells: Wall Length, Centre Distance
  and Cell Wall Midpoint Ratio
- [`projection_distortion()`](https://gillescolling.com/hexify/reference/projection_distortion.md)
  : Distortion of a Grid's Face Projection
- [`projection_stages()`](https://gillescolling.com/hexify/reference/projection_stages.md)
  : Points Through the Steps of Snyder's Projection

## Visualization

Plotting functions

- [`hexify_heatmap()`](https://gillescolling.com/hexify/reference/hexify_heatmap.md)
  : Create a ggplot2 Visualization of Hexagonal Grid Cells
- [`plot_grid()`](https://gillescolling.com/hexify/reference/plot_grid.md)
  : Plot Hexagonal Grid Clipped to a Polygon Boundary
- [`plot_world()`](https://gillescolling.com/hexify/reference/plot_world.md)
  : Quick World Map Plot
- [`net_layout()`](https://gillescolling.com/hexify/reference/net_layout.md)
  : Lay Out the Faces of a Grid's Solid in the Plane
- [`net_project()`](https://gillescolling.com/hexify/reference/net_project.md)
  : Place Points on a Net Layout
- [`net_cells()`](https://gillescolling.com/hexify/reference/net_cells.md)
  : The Parts of Cells on a Net Layout
- [`hex_rays()`](https://gillescolling.com/hexify/reference/hex_rays.md)
  : Ray Glyphs for Cell Values
- [`hex_triangles()`](https://gillescolling.com/hexify/reference/hex_triangles.md)
  : The Triangles of Cells Towards Each Neighbour
- [`hex_globe()`](https://gillescolling.com/hexify/reference/hex_globe.md)
  : Interactive Globe of a Grid
- [`hex_globe_png()`](https://gillescolling.com/hexify/reference/hex_globe_png.md)
  : Save a Globe as a PNG Image
- [`hexify_world`](https://gillescolling.com/hexify/reference/hexify_world.md)
  : Simplified World Map
- [`globe_centers`](https://gillescolling.com/hexify/reference/globe_centers.md)
  : Globe Center Presets
- [`hex_browse()`](https://gillescolling.com/hexify/reference/hex_browse.md)
  : Interactive Hex Map

## dggridR Compatibility

Interoperability with dggridR

- [`as_dggrid()`](https://gillescolling.com/hexify/reference/as_dggrid.md)
  [`format(`*`<dggs>`*`)`](https://gillescolling.com/hexify/reference/as_dggrid.md)
  [`print(`*`<dggs>`*`)`](https://gillescolling.com/hexify/reference/as_dggrid.md)
  : Convert hexify Grid to 'dggridR'-compatible Grid Object
- [`from_dggrid()`](https://gillescolling.com/hexify/reference/from_dggrid.md)
  : Convert 'dggridR' Grid Object to hexify_grid
- [`dggrid_is_compatible()`](https://gillescolling.com/hexify/reference/dggrid_is_compatible.md)
  : Validate 'dggridR' Grid Compatibility with hexify
- [`dggrid_43h_sequence()`](https://gillescolling.com/hexify/reference/dggrid_43h_sequence.md)
  : Create DGGRID 43H Aperture Sequence
- [`hexify_grid()`](https://gillescolling.com/hexify/reference/hexify_grid.md)
  [`format(`*`<hexify_grid>`*`)`](https://gillescolling.com/hexify/reference/hexify_grid.md)
  [`print(`*`<hexify_grid>`*`)`](https://gillescolling.com/hexify/reference/hexify_grid.md)
  : Create a Legacy Grid Specification
- [`dgearthstat()`](https://gillescolling.com/hexify/reference/dgearthstat.md)
  : Get Grid Statistics for Whole-Body Coverage
- [`dgverify()`](https://gillescolling.com/hexify/reference/dgverify.md)
  : Verify Grid Object
- [`hexify_compare_resolutions()`](https://gillescolling.com/hexify/reference/hexify_compare_resolutions.md)
  : Compare Grid Resolutions

## Low-Level API

Direct access to projection and coordinate transforms

- [`hexify_forward()`](https://gillescolling.com/hexify/reference/hexify_forward.md)
  : Forward Face Projection
- [`hexify_forward_to_face()`](https://gillescolling.com/hexify/reference/hexify_forward_to_face.md)
  : Forward Projection to Specific Face
- [`hexify_inverse()`](https://gillescolling.com/hexify/reference/hexify_inverse.md)
  : Inverse Face Projection
- [`hexify_which_face()`](https://gillescolling.com/hexify/reference/hexify_which_face.md)
  : Determine Which Face Contains a Point
- [`hexify_build_icosa()`](https://gillescolling.com/hexify/reference/hexify_build_icosa.md)
  : Set the Default Icosahedron Orientation
- [`hexify_face_centers()`](https://gillescolling.com/hexify/reference/hexify_face_centers.md)
  : Get Face Centers of a Solid
- [`hexify_grid_rect()`](https://gillescolling.com/hexify/reference/hexify_grid_rect.md)
  : Generate a Rectangular Grid of Hexagon Polygons
- [`hexify_grid_global()`](https://gillescolling.com/hexify/reference/hexify_grid_global.md)
  : Generate a Global Grid of Hexagon Polygons
- [`hexify_grid_to_cell()`](https://gillescolling.com/hexify/reference/hexify_grid_to_cell.md)
  : Convert Longitude/Latitude to Cell ID Using a Grid Object
- [`hexify_grid_cell_to_lonlat()`](https://gillescolling.com/hexify/reference/hexify_grid_cell_to_lonlat.md)
  : Convert Cell ID to Longitude/Latitude Using a Grid Object
- [`hex_corners_to_sf()`](https://gillescolling.com/hexify/reference/hex_corners_to_sf.md)
  : Build an sf POLYGON from Six (lon, lat) Corner Pairs
- [`dg_closest_res_to_area()`](https://gillescolling.com/hexify/reference/dg_closest_res_to_area.md)
  : Find Closest Resolution for Target Cell Area
- [`hexify_area_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_area_to_eff_res.md)
  : Convert Area to Effective Resolution
- [`hexify_eff_res_to_area()`](https://gillescolling.com/hexify/reference/hexify_eff_res_to_area.md)
  : Convert Effective Resolution to Area
- [`hexify_eff_res_to_resolution()`](https://gillescolling.com/hexify/reference/hexify_eff_res_to_resolution.md)
  : Convert Effective Resolution to Index Resolution
- [`hexify_resolution_to_eff_res()`](https://gillescolling.com/hexify/reference/hexify_resolution_to_eff_res.md)
  : Convert Index Resolution to Effective Resolution
- [`hexify_compare_indices()`](https://gillescolling.com/hexify/reference/hexify_compare_indices.md)
  : Compare Two Indices
- [`hexify_default_index_type()`](https://gillescolling.com/hexify/reference/hexify_default_index_type.md)
  : Get Default Index Type for Aperture
- [`hexify_is_valid_index_type()`](https://gillescolling.com/hexify/reference/hexify_is_valid_index_type.md)
  : Check If Index Type Is Valid for Aperture
- [`cell_to_index()`](https://gillescolling.com/hexify/reference/cell_to_index.md)
  : Convert Cell ID to Hierarchical Index String
- [`index_to_cell()`](https://gillescolling.com/hexify/reference/index_to_cell.md)
  : Convert Hierarchical Index to Cell ID
- [`hexify_get_children()`](https://gillescolling.com/hexify/reference/hexify_get_children.md)
  : Get Children Indices
- [`hexify_get_parent()`](https://gillescolling.com/hexify/reference/hexify_get_parent.md)
  : Get Parent Index
- [`hexify_get_resolution()`](https://gillescolling.com/hexify/reference/hexify_get_resolution.md)
  : Get Index Resolution
- [`hexify_z7_canonical()`](https://gillescolling.com/hexify/reference/hexify_z7_canonical.md)
  : Get Canonical Form of Z7 Index
- [`hexify_roundtrip_test()`](https://gillescolling.com/hexify/reference/hexify_roundtrip_test.md)
  : Round-Trip Accuracy Test

## Internal

Package internals

- [`hexify-package`](https://gillescolling.com/hexify/reference/hexify-package.md)
  : hexify: Equal-Area Hex Grids on the Snyder ISEA Icosahedron
- [`new_hex_data()`](https://gillescolling.com/hexify/reference/new_hex_data.md)
  : Create a HexData Object (Internal)
- [`extract_grid()`](https://gillescolling.com/hexify/reference/extract_grid.md)
  : Extract Grid from Various Objects
- [`calculate_resolution_for_area()`](https://gillescolling.com/hexify/reference/calculate_resolution_for_area.md)
  : Calculate Resolution for Target Area
- [`hexify_grid_to_HexGridInfo()`](https://gillescolling.com/hexify/reference/hexify_grid_to_HexGridInfo.md)
  : Convert Legacy hexify_grid to HexGridInfo
- [`HexGridInfo_to_hexify_grid()`](https://gillescolling.com/hexify/reference/HexGridInfo_to_hexify_grid.md)
  : Convert HexGridInfo to Legacy hexify_grid
