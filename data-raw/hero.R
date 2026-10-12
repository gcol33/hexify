# Build the README hero image man/figures/hexify-hero.png.
#
# Three ISEA3H grids of different cell areas over South America, Africa and
# Europe, and East Asia, drawn as cell outlines on the land of hexify_world in
# plate carree.
#
# Run from the package root:  Rscript data-raw/hero.R
# Needs sf and ggplot2 only here; they are not package dependencies.

suppressMessages({
  library(sf)
  library(ggplot2)
})
pkgload::load_all(quiet = TRUE)

patches <- list(
  list(bbox = c(-80, -55, -35, 12), area_km2 = 150000, colour = "#1b9e77", lw = 0.7),
  list(bbox = c(-15, -10, 38, 50),  area_km2 = 40000,  colour = "#d95f02", lw = 0.5),
  list(bbox = c(90, 15, 140, 50),   area_km2 = 60000,  colour = "#7570b3", lw = 0.5)
)

p <- ggplot() +
  geom_sf(data = hexify_world, fill = "grey92", colour = "grey80", linewidth = 0.2)
for (pt in patches) {
  cells <- grid_rect(pt$bbox, hex_grid(area_km2 = pt$area_km2))
  p <- p + geom_sf(data = cells, fill = NA, colour = pt$colour, linewidth = pt$lw)
}
p <- p +
  coord_sf(xlim = c(-170, 180), ylim = c(-64, 82), expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("man/figures/hexify-hero.png", p, width = 15, height = 6.3, dpi = 100, bg = "white")
