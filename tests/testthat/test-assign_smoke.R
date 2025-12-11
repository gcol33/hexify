test_that("hex_assign runs and returns expected columns", {
  lon <- c(16.3738, 2.3522, -3.7038)
  lat <- c(48.2082, 48.8566, 40.4168)

  # choose eff_res from area (~2591 km², similar to res=9)
  eff_res <- hexify_area_to_eff_res(2591)

  out <- hex_assign(lon, lat, eff_res, make_polygons = FALSE)

  expect_true(all(c("id","face","eff_res","center_lon","center_lat") %in% names(out)))
  expect_equal(nrow(out), length(lon))
  expect_type(out$id, "character")

  # Polygons test (skip if sf not installed)
  skip_if_not_installed("sf")
  out_sf <- hex_assign(lon, lat, eff_res, make_polygons = TRUE)
  expect_s3_class(out_sf, "sf")
  expect_true(all(sf::st_is_valid(out_sf)))  # Fixed: use all() since st_is_valid returns a vector
})
