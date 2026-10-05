# Illustrative example of the SoftwareX paper: the population of the world's
# cities and towns (maps::world.cities) summed per equal-area cell. Prints every
# number the example section cites and writes the per-cell table and polygons
# that the paper's map figure is drawn from.
#
# Usage: Rscript paper/bench/example_world_cities.R

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "bench_common.R"))

# --- code shown in the paper -------------------------------------------------
cities <- maps::world.cities
grid <- hex_grid(area_km2 = 70000)
hd <- hexify(cities, lon = "long", lat = "lat", grid = grid)
pop <- hex_summarize(hd, people = sum(pop))
cells <- cell_to_sf(pop$cell_id, grid)
pop$parent <- get_parent(pop$cell_id, grid, levels = 2)
regions <- aggregate(people ~ parent, data = pop, FUN = sum)
h3 <- h3_crosswalk(pop$cell_id, grid)
mars <- hex_grid(area_km2 = 70000, radius_km = "mars")
# -----------------------------------------------------------------------------

top <- pop[which.max(pop$people), ]
print(grid)
cat(sprintf("places: %d; people: %s\n", nrow(cities), format(sum(cities$pop), big.mark = ",")))
cat(sprintf("grid: aperture %s, resolution %d, cell area %.0f km2, %s cells worldwide\n",
            grid@aperture, as.integer(grid@resolution), grid@area_km2, format(n_cells(grid), big.mark = ",")))
cat(sprintf("occupied cells: %d; parent regions: %d (cell area %.0f km2)\n",
            nrow(pop), nrow(regions), grid@area_km2 * 9))
in_top <- cities[hd@cell_id == top$cell_id, ]
cat(sprintf("most populous cell: %s people, centre %.2f E %.2f N, largest place %s\n",
            format(top$people, big.mark = ","), top$cell_cen_lon, top$cell_cen_lat,
            in_top$name[which.max(in_top$pop)]))
cat(sprintf("cells with fewer than 1,000 people: %d\n", sum(pop$people < 1000)))
cat(sprintf("cells holding half the people: %d\n",
            which(cumsum(sort(as.numeric(pop$people), decreasing = TRUE)) >=
                    sum(as.numeric(pop$people)) / 2)[1]))
cat(sprintf("H3 resolution chosen by crosswalk: %s; ISEA/H3 area ratio range %.2f-%.2f\n",
            paste(unique(h3r::getResolution(h3$h3_cell_id)), collapse = ","),
            min(h3$area_ratio), max(h3$area_ratio)))
cat(sprintf("Mars grid: resolution %s, cell area %.0f km2\n", mars@resolution, mars@area_km2))

write_result(sf::st_drop_geometry(pop), "example_world_cities")
sf::st_write(merge(cells, pop[, c("cell_id", "people")], by = "cell_id"),
             file.path(results_dir(), "example_world_cities.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
