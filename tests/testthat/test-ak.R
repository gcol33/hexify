# Kaseorg's octahedral projection ("ak")

test_that("the projection is libhex9's unwarped one", {
  d <- utils::read.csv(test_path("data", "libhex9_ak.csv"))
  ico <- projection_icosa("ak", "octahedron")
  err <- vapply(seq_len(nrow(d)), function(k) {
    f <- cpp_icosa_forward(ico, d$lon[k], d$lat[k])
    s <- cpp_face_tri_to_solid(ico, f[["face"]], f[["icosa_triangle_x"]],
                               f[["icosa_triangle_y"]])
    sqrt(sum((s[1, ] - c(d$x[k], d$y[k], d$z[k]))^2))
  }, numeric(1))
  expect_lt(max(err), 1e-12)
})

test_that("forward and inverse undo each other", {
  set.seed(12)
  lon <- runif(200, -180, 180)
  lat <- asin(runif(200, -1, 1)) * 180 / pi
  for (k in seq_along(lon)) {
    f <- hexify_forward(lon[k], lat[k], projection = "ak", polyhedron = "octahedron")
    b <- hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]], f[["face"]],
                        projection = "ak", polyhedron = "octahedron")
    expect_equal(unname(b), c(lon[k], lat[k]), tolerance = 1e-10)
  }
})

test_that("the projection takes only the octahedron", {
  expect_error(hex_grid(resolution = 3, projection = "ak"), "octahedron only")
  expect_error(hexify_forward(0, 0, projection = "ak"), "octahedron only")
  g <- hex_grid(resolution = 3, aperture = 4, polyhedron = "octahedron",
                projection = "ak")
  expect_equal(g@projection, "ak")
})

test_that("cells on it keep their IDs and vary in area", {
  g <- hex_grid(resolution = 2, aperture = 9, projection = "ak")
  ids <- bit64::as.integer64(seq_len(n_cells(g)))
  cc <- cell_to_lonlat(ids, g)
  expect_identical(lonlat_to_cell(cc$lon_deg, cc$lat_deg, g), ids)
  a <- cell_area(ids, g)
  expect_equal(sum(a), body_surface_km2(EARTH_RADIUS_KM), tolerance = 1e-6)
  expect_gt(max(a) / min(a), 1.2)
  expect_false(dggrs_definition(g)$dggh$definition$constraints$cellEqualSized)
})

# Hex9's own projection ("akw") reads libhex9's warp field, which hexify does
# not ship; these tests run where hex9_warp_download() has fetched it. The
# fixture holds, for 300 points, libhex9's place of each on the octahedron
# (hex9_project_sphere), its labels (hex9_encode_sphere) and the centre of its
# level-12 cell (paper/bench/hex9_libhex9_dump.cpp, modes bench and project).
skip_without_hex9_warp <- function() {
  skip_if_not(cpp_hex9_warp_ready() || file.exists(hex9_warp_path()),
              "Hex9's warp field is not downloaded")
}

test_that("without the warp field, projection akw says how to get it", {
  skip_if(cpp_hex9_warp_ready())
  op <- options(hexify.hex9_warp = tempfile())
  on.exit(options(op), add = TRUE)
  expect_error(hexify_forward(0, 0, projection = "akw", polyhedron = "octahedron"),
               "hex9_warp_download")
})

test_that("the warped projection is libhex9's", {
  skip_without_hex9_warp()
  d <- utils::read.csv(test_path("data", "libhex9_akw.csv"), colClasses = "character")
  ico <- projection_icosa("akw", "octahedron")
  err <- vapply(seq_len(nrow(d)), function(k) {
    f <- cpp_icosa_forward(ico, as.numeric(d$lon[k]), as.numeric(d$lat[k]))
    s <- cpp_face_tri_to_solid(ico, f[["face"]], f[["icosa_triangle_x"]],
                               f[["icosa_triangle_y"]])
    sqrt(sum((s[1, ] - as.numeric(c(d$x[k], d$y[k], d$z[k])))^2))
  }, numeric(1))
  expect_lt(max(err), 1e-12)
})

test_that("a Hex9 grid on it names points as libhex9 does", {
  skip_without_hex9_warp()
  d <- utils::read.csv(test_path("data", "libhex9_akw.csv"), colClasses = "character")
  lon <- as.numeric(d$lon)
  lat <- as.numeric(d$lat)
  for (L in c(0, 5, 12, 18)) {
    g <- hex_grid(resolution = L, aperture = 9, projection = "akw")
    expect_identical(cell_to_index(lonlat_to_cell(lon, lat, g), g), d[[paste0("L", L)]])
  }
  g <- hex_grid(resolution = 12, aperture = 9, projection = "akw")
  cc <- cell_to_lonlat(lonlat_to_cell(lon, lat, g), g)
  expect_equal(cc$lon_deg, as.numeric(d$lon12), tolerance = 1e-9)
  expect_equal(cc$lat_deg, as.numeric(d$lat12), tolerance = 1e-9)
})

test_that("the warped projection inverts and nearly equalises cell areas", {
  skip_without_hex9_warp()
  set.seed(13)
  lon <- runif(100, -180, 180)
  lat <- asin(runif(100, -1, 1)) * 180 / pi
  for (k in seq_along(lon)) {
    f <- hexify_forward(lon[k], lat[k], projection = "akw", polyhedron = "octahedron")
    b <- hexify_inverse(f[["icosa_triangle_x"]], f[["icosa_triangle_y"]], f[["face"]],
                        projection = "akw", polyhedron = "octahedron")
    expect_equal(unname(b), c(lon[k], lat[k]), tolerance = 1e-10)
  }
  g <- hex_grid(resolution = 3, aperture = 9, projection = "akw")
  ids <- bit64::as.integer64(seq_len(n_cells(g)))
  a <- cell_area(ids, g)
  expect_equal(sum(a), body_surface_km2(EARTH_RADIUS_KM), tolerance = 1e-6)
  expect_lt(max(a) / min(a), 1.001)
  expect_identical(lonlat_to_cell(cell_to_lonlat(ids, g)$lon_deg,
                                  cell_to_lonlat(ids, g)$lat_deg, g), ids)
})

test_that("on WGS84 a Hex9 grid names points as libhex9's geodetic functions do", {
  skip_without_hex9_warp()
  # The same 300 points read as geodetic degrees: hex9_encode(), and the
  # level-12 centre through hex9_unproject() (dump mode bench, SPHERE 0)
  d <- utils::read.csv(test_path("data", "libhex9_akw_wgs84.csv"),
                       colClasses = "character")
  lon <- as.numeric(d$lon)
  lat <- as.numeric(d$lat)
  for (L in c(0, 5, 12, 18)) {
    g <- hex_grid(resolution = L, aperture = 9, projection = "akw", ellipsoid = "WGS84")
    expect_identical(cell_to_index(lonlat_to_cell(lon, lat, g), g), d[[paste0("L", L)]])
  }
  g <- hex_grid(resolution = 12, aperture = 9, projection = "akw", ellipsoid = "WGS84")
  cc <- cell_to_lonlat(lonlat_to_cell(lon, lat, g), g)
  expect_equal(cc$lon_deg, as.numeric(d$lon12), tolerance = 1e-9)
  expect_equal(cc$lat_deg, as.numeric(d$lat12), tolerance = 1e-9)
})
