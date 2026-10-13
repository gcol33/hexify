# Build the README banner man/figures/hexify-hero.png.
#
# The land of hexify_world covered by aperture-3 ISEA grids whose cells shrink
# threefold from one continent group to the next, west to east, drawn on the
# Equal Earth projection so equal-area cells look equal in size.
#
# Run from the package root:  Rscript data-raw/hero.R
# Needs sf and ggplot2 only here; they are not package dependencies.

suppressMessages({
  library(sf)
  library(ggplot2)
})
pkgload::load_all(quiet = TRUE)
sf_use_s2(FALSE)

groups <- list(
  list(continents = c("North America", "South America"), resolution = 5,
       fill = "#3b1f5c", text = "#3b1f5c"),
  list(continents = c("Europe", "Africa"),               resolution = 6,
       fill = "#2c6e8f", text = "#2c6e8f"),
  list(continents = "Asia",                              resolution = 7,
       fill = "#1f9e89", text = "#16806f"),
  list(continents = "Oceania",                           resolution = 8,
       fill = "#6cc24a", text = "#3f8f22")
)

crs_ee <- "+proj=eqearth +lon_0=0 +datum=WGS84"

# One polygon per landmass piece; pieces west of 25 W (French Guiana, the
# French Caribbean) go with the Americas whatever country holds them.
land <- st_make_valid(hexify_world[hexify_world$continent != "Antarctica", "continent"]) |>
  st_cast("MULTIPOLYGON") |> st_cast("POLYGON", warn = FALSE)
west <- st_coordinates(suppressWarnings(st_centroid(land)))[, 1] < -25
land$continent[west] <- "North America"

cells <- do.call(rbind, lapply(seq_along(groups), function(i) {
  gr <- groups[[i]]
  piece <- st_union(land[land$continent %in% gr$continents, ])
  g <- hex_grid(resolution = gr$resolution, aperture = 3)
  cells <- grid_clip(piece, g, crop = FALSE)
  lon <- cell_to_lonlat(cells$cell_id, g)$lon
  cells <- cells[lon > -165 & lon < 170, ]
  cells$group <- i
  st_transform(st_wrap_dateline(cells["group"]), crs_ee)
}))

bb <- st_bbox(cells)
ylim <- st_coordinates(st_transform(st_sfc(st_point(c(0, -62)), st_point(c(0, 82)),
                                           crs = 4326), crs_ee))[, 2]

# The banner is 1628 x 720 px at 120 dpi, the size of the other package
# banners; the frame is widened about the cells so the panel inside the
# margins has the banner's shape at a 1:1 scale.
dpi <- 120
width <- 1628 / dpi
height <- 720 / dpi
margin_in <- 2 * 8 / 72
aspect <- (width - margin_in) / (height - margin_in)
xlim <- mean(c(bb[["xmin"]], bb[["xmax"]])) + c(-0.5, 0.5) * aspect * diff(ylim)

# Cell areas written in open ocean next to each group, placed as fractions of
# the cells' extent (x from the left) and of the frame (y from the bottom),
# large enough for phone width.
label_at <- rbind(c(0.085, 0.27), c(0.40, 0.08), c(0.69, 0.25), c(0.84, 0.05))
key_text <- data.frame(
  x = bb[["xmin"]] + label_at[, 1] * (bb[["xmax"]] - bb[["xmin"]]),
  y = ylim[1] + label_at[, 2] * diff(ylim),
  group = seq_along(groups),
  label = vapply(groups, function(gr) {
    a <- hex_grid(resolution = gr$resolution, aperture = 3)@area_km2
    paste0(format(signif(a, 2) / 1000, scientific = FALSE, drop0trailing = TRUE), "k km²")
  }, character(1))
)

fills <- vapply(groups, `[[`, "", "fill")
p <- ggplot() +
  geom_sf(data = cells, aes(fill = factor(group)), colour = "white", linewidth = 0.15) +
  geom_text(data = key_text, aes(x, y, label = label, colour = factor(group)),
            size = 10, fontface = "bold") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_colour_manual(values = vapply(groups, `[[`, "", "text"), guide = "none") +
  coord_sf(crs = crs_ee, xlim = xlim, ylim = ylim,
           expand = FALSE, datum = NA) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 8, 8, 8))

ggsave("man/figures/hexify-hero.png", p, width = width, height = height, dpi = dpi, bg = "white")
file.copy("man/figures/hexify-hero.png", "docs/reference/figures/hexify-hero.png",
          overwrite = TRUE)
