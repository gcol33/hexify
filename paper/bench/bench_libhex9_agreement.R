# Hex9 in hexify against libhex9, Griffin's reference implementation (2.3.1).
#
# For random points uniform on the sphere, at every level 0 to 18: the cell
# of hexify's Hex9 grid on Hex9's own projection (aperture 9, projection
# "akw") holding each point, as its libhex9 label, against the label
# libhex9's sphere functions give it (hex9_encode_sphere, hex9_bin); and the
# centre hexify gives that cell against libhex9's lattice centre of it
# (hex9_cell_uv, unprojected by hex9_unproject_sphere). The same again with
# the points' degrees read as geodetic on WGS84: hexify's grid with
# ellipsoid = "WGS84" against libhex9's geodetic functions (hex9_encode,
# hex9_unproject), centres compared along the great circle of their
# geodetic coordinates on the authalic sphere. Also the projection alone:
# each point's place on the flat octahedron under hexify's "akw" against
# libhex9's (hex9_project_sphere).
#
# Needs the dump tool built from hex9_libhex9_dump.cpp (see its header),
# pointed at by LIBHEX9_DUMP, and libhex9's warp field (hex9_warp_download()).
#
# Usage: Rscript paper/bench/bench_libhex9_agreement.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

dump <- Sys.getenv("LIBHEX9_DUMP")
if (!nzchar(dump)) stop("point LIBHEX9_DUMP at the hex9_libhex9_dump tool")
N_POINTS <- 20000
LEVELS <- 0:18
RADIUS <- 6371.007180918475
pts <- sphere_points(N_POINTS, seed = 100)

input <- tempfile(fileext = ".txt")
writeLines(sprintf("%.17g %.17g", pts$lon, pts$lat), input)
libhex9 <- function(sphere) {
  read.csv(text = system2(dump, c("bench", sphere, LEVELS), stdin = input, stdout = TRUE),
           colClasses = "character")
}
ref <- list(sphere = libhex9(1), wgs84 = libhex9(0))

agreement <- function(L, earth) {
  g <- hex_grid(resolution = L, aperture = 9, projection = "akw",
                ellipsoid = if (earth == "wgs84") "WGS84" else NULL)
  r <- ref[[earth]]
  ids <- lonlat_to_cell(pts$lon, pts$lat, g)
  agree <- cell_to_index(ids, g) == r[[paste0("L", L)]]
  cc <- cell_to_lonlat(ids[agree], g)
  d <- gc_km(cc$lon_deg, cc$lat_deg, as.numeric(r[[paste0("lon", L)]][agree]),
             as.numeric(r[[paste0("lat", L)]][agree]), radius = RADIUS)
  c(labels_agree = sum(agree), centre_m_median = 1000 * median(d),
    centre_m_max = 1000 * max(d))
}

rows <- lapply(LEVELS, function(L) {
  s <- agreement(L, "sphere")
  w <- agreement(L, "wgs84")
  r <- data.frame(level = L, points = N_POINTS, labels_agree = s[["labels_agree"]],
                  share_agree = s[["labels_agree"]] / N_POINTS,
                  centre_m_median = s[["centre_m_median"]], centre_m_max = s[["centre_m_max"]],
                  wgs84_labels_agree = w[["labels_agree"]],
                  wgs84_centre_m_median = w[["centre_m_median"]],
                  wgs84_centre_m_max = w[["centre_m_max"]],
                  cell_km = sqrt(4 * pi * RADIUS^2 / (12 * 9^L)))
  print(r, digits = 4)
  r
})
out <- do.call(rbind, rows)

proj <- read.csv(text = system2(dump, c("project", 1, 1), stdin = input, stdout = TRUE))
ico <- hexify:::projection_icosa("akw", "octahedron")
gap <- vapply(seq_len(N_POINTS), function(k) {
  f <- hexify:::cpp_icosa_forward(ico, pts$lon[k], pts$lat[k])
  s <- hexify:::cpp_face_tri_to_solid(ico, f[["face"]], f[["icosa_triangle_x"]],
                                      f[["icosa_triangle_y"]])
  sqrt(sum((s[1, ] - c(proj$x[k], proj$y[k], proj$z[k]))^2))
}, numeric(1))
cat("largest gap between the two projections on the octahedron:", max(gap), "\n")
out$projection_gap_max <- max(gap)

write_result(out, "libhex9_agreement",
             extra = c(sprintf("points: %d uniform on the sphere, seed 100", N_POINTS),
                       "libhex9: 2.3.1, commit 1eb0be7, sphere and WGS84 functions",
                       paste("dump tool:", dump)))
