test_that("every sampled point lies in its cell", {
  set.seed(42)
  grids <- list(hex_grid(resolution = 3, aperture = 3),
                hex_grid(resolution = 3, aperture = 4),
                hex_grid(resolution = 2, aperture = 7),
                hex_grid(resolution = 4, aperture = "4/3"),
                suppressMessages(hex_grid(resolution = 2, type = "h3")))
  for (g in grids) {
    cells <- lonlat_to_cell(c(10, -120, 0, 11.25, 179.9), c(50, -40, 90, 58.28252559, 0), g)
    pts <- hex_sample(cells, g, n = 25)
    expect_equal(pts$cell_id, rep(cells, each = 25))
    expect_equal(lonlat_to_cell(pts$lon, pts$lat, g), pts$cell_id)
  }
})

test_that("n takes one count per cell, including zero", {
  set.seed(1)
  g <- hex_grid(resolution = 3, aperture = 3)
  cells <- c(5, 17, 40)
  pts <- hex_sample(cells, g, n = c(3, 0, 2))
  expect_equal(pts$cell_id, c(5, 5, 5, 40, 40))
  expect_equal(nrow(hex_sample(cells, g, n = 0)), 0L)
  expect_error(hex_sample(cells, g, n = c(1, 2)), "one per cell")
  expect_error(hex_sample(cells, g, n = 1.5), "whole number")
})

test_that("the cap around each centre holds the whole cell", {
  for (g in list(hex_grid(resolution = 3, aperture = 3),
                 hex_grid(resolution = 2, aperture = 7),
                 suppressMessages(hex_grid(resolution = 2, type = "h3")))) {
    cells <- lonlat_to_cell(c(10, 0, 11.25, -170), c(50, 90, 58.28252559, -60), g)
    ctr <- cell_to_lonlat(cells, g)
    C <- hexify:::unit_vec(ctr$lon_deg, ctr$lat_deg)
    cap <- hexify:::cell_cap_radius(cells, g, C)
    fine <- cell_to_sf(cells, g, wrap_dateline = FALSE, densify = 1e-6)
    for (i in seq_along(cells)) {
      xy <- sf::st_coordinates(sf::st_geometry(fine)[[i]])
      B <- hexify:::unit_vec(xy[, 1], xy[, 2])
      expect_lte(max(acos(pmin(1, B %*% C[i, ]))), cap[i])
    }
  }
})

test_that("points spread evenly over the cell's area", {
  set.seed(7)
  g <- hex_grid(resolution = 2, aperture = 4)
  cell <- lonlat_to_cell(10, 50, g)
  pts <- hex_sample(cell, g, n = 20000)

  # ISEA is equal-area, so the cells of a finer grid that lie wholly inside
  # the sampled cell all have the same area and expect the same count
  fine <- hex_grid(resolution = 5, aperture = 4)
  fine_cell <- lonlat_to_cell(pts$lon, pts$lat, fine)
  ids <- unique(fine_cell)
  rings <- hexify:::isea_cell_rings(ids, 5, "4", numeric(0), 0)
  inside <- vapply(rings, function(r) {
    all(lonlat_to_cell(r[, 1], r[, 2], g) == cell) &&
      all(lonlat_to_cell(r[, 1] * 0.999 + mean(r[, 1]) * 0.001,
                         r[, 2] * 0.999 + mean(r[, 2]) * 0.001, g) == cell)
  }, logical(1))
  counts <- tabulate(match(fine_cell, ids[inside]), nbins = sum(inside))
  expect_gt(length(counts), 20)
  expect_gt(suppressWarnings(stats::chisq.test(counts)$p.value), 1e-3)
})
