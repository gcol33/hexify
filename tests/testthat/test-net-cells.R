# Cells on net layouts, glue tabs, parent outlines and area legends (#90)

solid_angles <- function(g, ids) {
  lv <- isea_levels(g@aperture, g@resolution)
  cpp_cell_solid_angle(icosa_arg(g), as_cell_id(ids), lv$resolution, lv$aperture,
                       lv$ap_seq, CELL_WALL_TOLERANCE)
}

test_that("the parts of a cell add up to it on an equal-area net", {
  cases <- list(
    list(grid = hex_grid(resolution = 2, aperture = 3), layout = "plane"),
    list(grid = hex_grid(resolution = 2, aperture = 4, polyhedron = "octahedron",
                         orientation = "gosper"), layout = "gosper"),
    list(grid = hex_grid(resolution = 2, aperture = 7), layout = "land")
  )
  for (cs in cases) {
    net <- net_layout(cs$grid, cs$layout)
    parts <- net_cells(net, cs$grid)
    area <- tapply(parts$area, parts$cell_id, sum)
    ids <- as.numeric(names(area))
    edge <- face_plane_edge(icosa_arg(cs$grid))
    expect_equal(as.vector(area) * edge^2, solid_angles(cs$grid, ids), tolerance = 1e-6)
    expect_setequal(ids, seq_len(n_cells(cs$grid)))
    n_faces <- nrow(icosa_solid(icosa_arg(cs$grid))$faces)
    expect_equal(sum(parts$area), n_faces * sqrt(3) / 4, tolerance = 1e-9)
  }
})

test_that("seams split cells into groups the map shows apart", {
  g <- hex_grid(resolution = 3, aperture = 3)
  parts <- net_cells(net_layout(g, "plane"), g)
  groups <- tapply(parts$group, parts$cell_id, max)
  pieces <- tapply(parts$piece, parts$cell_id, function(p) length(unique(p)))
  # Every cut cell lies on more than one piece, and some cells cross a joined
  # face edge without being cut.
  expect_true(all(pieces[groups > 1] > 1))
  expect_gt(sum(groups > 1), 0)
  expect_gt(sum(groups == 1 & pieces > 1), 0)
  # Cell 1 is the pentagon at vertex 0, where the five top faces of the
  # plane net end in separate tips.
  expect_true(is_pentagon(1, g))
  expect_equal(unname(groups[["1"]]), 5)
})

test_that("a net laid out on a tree of faces is cut along the rest of the edges", {
  g <- hex_grid(resolution = 1, aperture = 3)
  net <- net_layout(g, "land")
  tabs <- net_tabs(net)
  # 30 edges, 19 of them kept by a spanning tree of the 20 faces
  expect_length(tabs, 11L)
  for (t in tabs) expect_equal(nrow(t), 5L)
})

test_that("the plot method draws parents, the area legend and tabs", {
  pdf(NULL)
  on.exit(dev.off())
  g <- hex_grid(resolution = 3, aperture = 3)
  expect_invisible(plot(g, surface = "net", land = FALSE, parents = 1:2,
                        parent_col = c("black", "red"), area_legend = TRUE))
  expect_invisible(plot(g, land = FALSE, parents = 1))
  expect_invisible(plot(g, surface = "solid", land = FALSE, parents = 1,
                        cells = 1:20))
  expect_invisible(plot(g, surface = "net", layout = "land", land = FALSE,
                        tabs = TRUE))
  expect_invisible(plot(g, surface = "net", layout = "rhombic", land = FALSE,
                        area_legend = TRUE))
  expect_error(plot(g, parents = 4, land = FALSE), "parents")
  expect_error(plot(g, tabs = TRUE, land = FALSE), "net")
  gf <- hex_grid(resolution = 3, aperture = 3, projection = "fuller")
  expect_error(plot(gf, surface = "net", land = FALSE, area_legend = TRUE),
               "equal-area")
})

test_that("parent outlines follow the hierarchy of the drawn cells", {
  g <- hex_grid(resolution = 3, aperture = 3)
  cells <- 1:30
  h <- parent_levels(1, g, cells, 0.05, "black", 1)
  expect_length(h, 1L)
  expect_equal(h[[1]]$grid@resolution, 2L)
  n_parents <- length(unique(get_parent(cells, g)))
  expect_equal(length(unique(h[[1]]$paths[, "cell"])), n_parents)
})
