# Neighbour smoothing: each step is the weighted mean of a cell and its
# listed neighbours.

test_that("one step is the mean of a cell and its neighbours", {
  g <- hex_grid(resolution = 3, aperture = 3)
  ids <- grid_global(g)$cell_id
  set.seed(2)
  v <- rnorm(length(ids))
  s <- hex_smooth(ids, v, g)
  nb <- get_neighbors(ids, g)
  expect_equal(s, vapply(seq_along(ids), function(k)
    mean(v[match(c(ids[k], nb[[k]]), ids)]), numeric(1)))

  s0 <- hex_smooth(ids, v, g, self_weight = 0)
  expect_equal(s0, vapply(nb, function(x) mean(v[match(x, ids)]), numeric(1)))

  s2 <- hex_smooth(ids, v, g, self_weight = 2)
  expect_equal(s2, vapply(seq_along(ids), function(k)
    (2 * v[k] + sum(v[match(nb[[k]], ids)])) / (2 + length(nb[[k]])), numeric(1)))
})

test_that("steps repeat the smoother, constants stay and the spread shrinks", {
  g <- hex_grid(resolution = 3, aperture = 4)
  ids <- grid_global(g)$cell_id
  set.seed(3)
  v <- rnorm(length(ids))
  expect_equal(hex_smooth(ids, v, g, steps = 3),
               hex_smooth(ids, hex_smooth(ids, v, g, steps = 2), g))
  expect_equal(hex_smooth(ids, rep(4, length(ids)), g, steps = 5),
               rep(4, length(ids)))
  expect_lt(sd(hex_smooth(ids, v, g, steps = 3)), sd(v) / 2)
  expect_identical(hex_smooth(ids, v, g, steps = 0), v)
})

test_that("unlisted cells, NA values and other groups do not contribute", {
  g <- hex_grid(resolution = 4, aperture = 7)
  centre <- lonlat_to_cell(10, 45, g)
  nb <- get_neighbors(centre, g)[[1]]
  ids <- c(centre, nb[1:3])
  expect_equal(hex_smooth(ids, c(0, 3, 6, 9), g)[1], mean(c(0, 3, 6, 9)))

  v <- c(0, 3, NA, 9)
  s <- hex_smooth(ids, v, g)
  expect_equal(s[1], mean(c(0, 3, 9)))
  expect_true(is.na(s[3]))

  s <- hex_smooth(ids, c(0, 3, 6, 9), g, group = c("a", "a", "b", NA))
  expect_equal(s[1], mean(c(0, 3)))
  expect_equal(s[3], 6)
  expect_equal(s[4], 9)

  # A lone cell with self weight 0 has nothing to average and keeps its value
  expect_equal(hex_smooth(centre, 5, g, self_weight = 0), 5)
})

test_that("pentagons average over five neighbours, H3 cells too", {
  g <- hex_grid(resolution = 2, aperture = 3)
  ids <- grid_global(g)$cell_id
  v <- seq_along(ids)
  pent <- ids[is_pentagon(ids, g)][1]
  nb <- get_neighbors(pent, g)[[1]]
  expect_length(nb, 5L)
  k <- match(pent, ids)
  expect_equal(hex_smooth(ids, v, g)[k], mean(v[match(c(pent, nb), ids)]))

  h <- hex_grid(resolution = 1, type = "h3")
  hid <- grid_global(h)$cell_id
  hv <- seq_along(hid)
  hn <- get_neighbors(hid[1], h)[[1]]
  expect_equal(hex_smooth(hid, hv, h)[1], mean(hv[match(c(hid[1], hn), hid)]))
})

test_that("bad input is refused", {
  g <- hex_grid(resolution = 2, aperture = 3)
  ids <- grid_global(g)$cell_id[1:3]
  expect_error(hex_smooth(ids, 1:2, g), "same length")
  expect_error(hex_smooth(c(ids, ids[1]), 1:4, g), "each cell once")
  expect_error(hex_smooth(ids, letters[1:3], g), "numeric")
  expect_error(hex_smooth(ids, 1:3, g, steps = 1.5), "steps")
  expect_error(hex_smooth(ids, 1:3, g, self_weight = -1), "self_weight")
  expect_error(hex_smooth(ids, 1:3, g, group = 1:2), "one entry per cell")
})
