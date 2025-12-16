# hexify_map.R
# Visualization functions with basemap support

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Calculate view buffer for bounding box
#' @noRd
calculate_view_buffer <- function(bbox, factor = 0.1, min_buffer = 1) {
  x_range <- bbox["xmax"] - bbox["xmin"]
  y_range <- bbox["ymax"] - bbox["ymin"]
  list(
    x = max(x_range * factor, min_buffer),
    y = max(y_range * factor, min_buffer)
  )
}

#' Check if string is a color palette name (not a hex color or named color)
#' @noRd
is_palette_name <- function(x) {
  length(x) == 1 &&
    is.character(x) &&
    !grepl("^#", x) &&
    !x %in% grDevices::colors()
}

#' Check if palette is a valid RColorBrewer palette
#' @noRd
is_brewer_palette <- function(palette_name) {
  requireNamespace("RColorBrewer", quietly = TRUE) &&
    palette_name %in% row.names(RColorBrewer::brewer.pal.info)
}

#' Apply discrete color scale to ggplot
#' @noRd
apply_discrete_scale <- function(p, colors, legend_title, na_color, n_levels) {
  if (is.null(colors)) {
    return(p + ggplot2::scale_fill_viridis_d(
      name = legend_title, na.value = na_color
    ))
  }

  if (!is_palette_name(colors)) {
    return(p + ggplot2::scale_fill_manual(
      values = colors, name = legend_title, na.value = na_color
    ))
  }

  if (is_brewer_palette(colors)) {
    max_colors <- RColorBrewer::brewer.pal.info[colors, "maxcolors"]
    pal_colors <- RColorBrewer::brewer.pal(min(n_levels, max_colors), colors)
    return(p + ggplot2::scale_fill_manual(
      values = pal_colors, name = legend_title, na.value = na_color
    ))
  }

  # Fallback: treat as viridis option name
  p + ggplot2::scale_fill_viridis_d(
    option = tolower(colors), name = legend_title, na.value = na_color
  )
}

#' Apply continuous color scale to ggplot
#' @noRd
apply_continuous_scale <- function(p, colors, legend_title, na_color) {
  if (is.null(colors)) {
    return(p + ggplot2::scale_fill_viridis_c(
      name = legend_title, na.value = na_color
    ))
  }

  if (!is_palette_name(colors)) {
    return(p + ggplot2::scale_fill_gradientn(
      colors = colors, name = legend_title, na.value = na_color
    ))
  }

  if (is_brewer_palette(colors)) {
    return(p + ggplot2::scale_fill_distiller(
      palette = colors, direction = 1,
      name = legend_title, na.value = na_color
    ))
  }

  # Fallback: treat as viridis option name
  p + ggplot2::scale_fill_viridis_c(
    option = tolower(colors), name = legend_title, na.value = na_color
  )
}

#' Convert hexify data to sf polygons
#' @noRd
prepare_hex_sf <- function(data, aperture) {
  if (inherits(data, "sf")) return(data)

  if (!is.data.frame(data) || !"cell_id" %in% names(data)) {
    stop("data must be a data.frame from hexify() or an sf object")
  }

  if (!"cell_area" %in% names(data)) {
    stop("data must contain 'cell_area' column (output from hexify()). ",
         "Use hexify_polygons() directly if you have pre-computed polygons.")
  }

  hex_sf <- hexify_to_polygons(data, aperture = aperture, return_sf = TRUE)

  # Join extra columns from original data

  extra_cols <- setdiff(names(data), c("cell_id", "geometry"))
  if (length(extra_cols) > 0) {
    cols <- c("cell_id", extra_cols)
    data_unique <- data[!duplicated(data$cell_id), cols, drop = FALSE]
    hex_sf <- merge(hex_sf, data_unique, by = "cell_id", all.x = TRUE)
  }

  hex_sf
}

#' Resolve basemap specification to sf object
#' @noRd
resolve_basemap <- function(basemap) {
  if (is.null(basemap)) return(NULL)

  if (is.character(basemap) && basemap == "world") {
    return(hexify_world)
  }

  if (is.character(basemap) && basemap == "world_hires") {
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
      stop("Package 'rnaturalearth' is required for basemap = 'world_hires'. ",
           "Install with: install.packages('rnaturalearth')")
    }
    return(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"))
  }

  if (inherits(basemap, "sf") || inherits(basemap, "sfc")) {
    return(basemap)
  }

  stop("basemap must be 'world', 'world_hires', or an sf object")
}

#' Create mask for areas outside basemap
#' @noRd
create_outside_mask <- function(hex_sf, basemap_sf, xlim, ylim) {
  if (!is.null(xlim) && !is.null(ylim)) {
    bbox_coords <- c(
      xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]
    )
  } else {
    hex_bbox <- sf::st_bbox(hex_sf)
    buffer <- calculate_view_buffer(hex_bbox)
    bbox_coords <- c(
      xmin = unname(hex_bbox["xmin"] - buffer$x),
      xmax = unname(hex_bbox["xmax"] + buffer$x),
      ymin = unname(hex_bbox["ymin"] - buffer$y),
      ymax = unname(hex_bbox["ymax"] + buffer$y)
    )
  }

  bbox_poly <- sf::st_as_sfc(sf::st_bbox(bbox_coords, crs = sf::st_crs(hex_sf)))
  basemap_union <- sf::st_union(sf::st_make_valid(basemap_sf))
  suppressWarnings(sf::st_difference(bbox_poly, basemap_union))
}

#' Build heatmap layers with masking
#' @noRd
build_masked_layers <- function(p, hex_sf, fill_col, hex_border, hex_lwd,
                                hex_alpha, basemap_sf, basemap_border,
                                basemap_lwd, mask_sf) {
  # Hexagons first

  p <- p + ggplot2::geom_sf(
    data = hex_sf,
    ggplot2::aes(fill = .data[[fill_col]]),
    color = hex_border,
    linewidth = hex_lwd,
    alpha = hex_alpha
  )

  # Mask layer to hide hexes outside land
  if (!is.null(mask_sf)) {
    p <- p + ggplot2::geom_sf(data = mask_sf, fill = "white", color = NA)
  }

  # Basemap borders on top (no fill)
  p + ggplot2::geom_sf(
    data = basemap_sf,
    fill = NA,
    color = basemap_border,
    linewidth = basemap_lwd
  )
}

#' Build standard heatmap layers (basemap under hexes)
#' @noRd
build_standard_layers <- function(p, hex_sf, fill_col, hex_border, hex_lwd,
                                  hex_alpha, basemap_sf, basemap_fill,
                                  basemap_border, basemap_lwd) {
  if (!is.null(basemap_sf)) {
    p <- p + ggplot2::geom_sf(
      data = basemap_sf,
      fill = basemap_fill,
      color = basemap_border,
      linewidth = basemap_lwd
    )
  }

  p + ggplot2::geom_sf(
    data = hex_sf,
    ggplot2::aes(fill = .data[[fill_col]]),
    color = hex_border,
    linewidth = hex_lwd,
    alpha = hex_alpha
  )
}

#' Simple sf preparation for hexify_map (no extra column merging)
#' @noRd
prepare_hex_sf_simple <- function(data, aperture) {
  if (inherits(data, "sf")) return(data)

  if (!is.data.frame(data) || !"cell_id" %in% names(data)) {
    stop("data must be a data.frame from hexify() or an sf object")
  }

  if (!"cell_area" %in% names(data)) {
    stop("data must contain 'cell_area' column (output from hexify()). ",
         "Use hexify_polygons() directly if you have pre-computed polygons.")
  }

  hexify_to_polygons(data, aperture = aperture, return_sf = TRUE)
}

#' Resolve value column name (auto-detect if NULL)
#' @noRd
resolve_value_column <- function(hex_sf, value) {
  if (!is.null(value)) {
    if (!value %in% names(hex_sf)) {
      stop("Column '", value, "' not found in data. ",
           "Available columns: ", paste(names(hex_sf), collapse = ", "))
    }
    return(value)
  }

  if ("count" %in% names(hex_sf)) return("count")
  if ("n" %in% names(hex_sf)) return("n")

  stop("No 'value' column specified and no 'count' or 'n' column found in data")
}

#' Prepare fill column with optional binning
#' @noRd
prepare_fill_column <- function(hex_sf, value, breaks, labels) {
  value_data <- hex_sf[[value]]
  is_discrete <- is.factor(value_data) || is.character(value_data)

  if (is.null(breaks) || is_discrete) {
    return(list(data = hex_sf, fill_col = value, is_discrete = is_discrete))
  }

  # Generate labels if not provided
  if (is.null(labels)) {
    labels <- generate_bin_labels(breaks)
  }

  bin_col <- paste0(value, "_bin")
  hex_sf[[bin_col]] <- cut(
    value_data, breaks = breaks, labels = labels, include.lowest = TRUE
  )

  list(data = hex_sf, fill_col = bin_col, is_discrete = TRUE)
}

#' Generate bin labels from break points
#' @noRd
generate_bin_labels <- function(breaks) {
  n_bins <- length(breaks) - 1
  labels <- character(n_bins)

  for (i in seq_len(n_bins)) {
    low <- breaks[i]
    high <- breaks[i + 1]
    labels[i] <- if (is.infinite(low) && low < 0) {
      paste0("<", high)
    } else if (is.infinite(high)) {
      paste0(">", low)
    } else {
      paste0(low, "-", high)
    }
  }

  labels
}

#' Resolve basemap specification including raster support
#' @noRd
resolve_basemap_with_raster <- function(basemap) {
  result <- list(sf = NULL, raster = NULL)
  if (is.null(basemap)) return(result)

  if (is.character(basemap) && basemap == "world") {
    result$sf <- hexify_world
    return(result)
  }

  if (inherits(basemap, "sf") || inherits(basemap, "sfc")) {
    result$sf <- basemap
    return(result)
  }

  if (inherits(basemap, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required for SpatRaster basemaps")
    }
    result$raster <- basemap
    return(result)
  }

  if (inherits(basemap, c("RasterLayer", "RasterBrick", "RasterStack"))) {
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package 'raster' is required for Raster* basemaps")
    }
    result$raster <- basemap
    return(result)
  }

  stop("basemap must be 'world', an sf object, or a raster (terra/raster)")
}

# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================

#' Plot hexagonal grid cells with optional basemap
#'
#' Creates a map visualization of hexagonal grid cells. Supports the built-in
#' world map or user-supplied basemaps (sf vectors or raster images).
#'
#' @param data Data frame from hexify() containing cell_id and cell_area columns,
#'   or an sf object with hexagon polygons
#' @param basemap Optional basemap. Can be:
#'   \itemize{
#'     \item \code{NULL}: No basemap (default)
#'     \item \code{"world"}: Use built-in \code{hexify_world} map
#'     \item An sf object: User-supplied vector map
#'     \item A SpatRaster (terra) or RasterLayer (raster): User-supplied raster
#'   }
#' @param fill Fill color for hexagons (single color or column name for mapping)
#' @param border Border color for hexagons
#' @param lwd Line width for hexagon borders
#' @param alpha Transparency for hexagon fill (0-1)
#' @param basemap_fill Fill color for basemap polygons (if vector)
#' @param basemap_border Border color for basemap polygons (if vector)
#' @param basemap_lwd Line width for basemap borders
#' @param aperture Grid aperture (default 3), used if data is from hexify()
#' @param xlim Optional x-axis (longitude) limits as c(min, max)
#' @param ylim Optional y-axis (latitude) limits as c(min, max)
#' @param main Plot title
#' @param ... Additional arguments passed to plot()
#'
#' @return NULL invisibly. Creates a plot as side effect.
#'
#' @details
#' This function provides a simple way to visualize hexagonal grids with
#' geographic context. For more sophisticated visualizations, use
#' \code{\link{hexify_to_polygons}} to get an sf object and plot with ggplot2,
#' tmap, or other mapping packages.
#'
#' The function automatically:
#' \itemize{
#'   \item Converts hexify() output to polygons if needed
#'   \item Adjusts aspect ratio for latitude
#'   \item Clips basemap to data extent (with buffer)
#' }
#'
#' @section Basemap options:
#' \describe{
#'   \item{Built-in world map}{Use \code{basemap = "world"} for the included
#'     simplified Natural Earth map. No additional packages required.}
#'   \item{Custom sf vector}{Pass any sf object as basemap for custom boundaries,
#'     regions, or detailed coastlines.
#'     install.packages("rnaturalearth")
#'   }
#'   \item{Raster basemap}{Pass a SpatRaster (terra) or RasterLayer (raster)
#'     for satellite imagery or other raster backgrounds. Requires terra or
#'     raster package.}
#' }
#'
#' @family visualization
#' @seealso \code{\link{hexify_plot}} for simple base R plots,
#'   \code{\link{hexify_heatmap}} for ggplot2-based heatmaps
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#'
#' # Sample data
#' cities <- data.frame(
#'   name = c("Vienna", "Paris", "Madrid"),
#'   lon = c(16.37, 2.35, -3.70),
#'   lat = c(48.21, 48.86, 40.42)
#' )
#' result <- hexify(cities, lon = "lon", lat = "lat", area = 5000)
#'
#' # Simple plot without basemap
#' hexify_map(result)
#'
#' # With built-in world map
#' hexify_map(result, basemap = "world")
#'
#' # Custom colors
#' hexify_map(result, basemap = "world",
#'            fill = "steelblue", border = "darkblue",
#'            basemap_fill = "ivory", basemap_border = "gray50")
#'
#' # With user-supplied sf basemap
#' library(rnaturalearth)
#' europe <- ne_countries(continent = "Europe", returnclass = "sf")
#' hexify_map(result, basemap = europe)
#'
#' # Zoom to specific region
#' hexify_map(result, basemap = "world",
#'            xlim = c(-10, 25), ylim = c(35, 55))
#' }
hexify_map <- function(data,
                    basemap = NULL,
                    fill = "steelblue",
                    border = "gray30",
                    lwd = 1,
                    alpha = 0.7,
                    basemap_fill = "gray95",
                    basemap_border = "gray70",
                    basemap_lwd = 1,
                    aperture = 3L,
                    xlim = NULL,
                    ylim = NULL,
                    main = NULL,
                    ...) {

  # Convert hexify output to polygons if needed
  hex_sf <- prepare_hex_sf_simple(data, aperture)

  # Get bounding box and calculate buffer
  hex_bbox <- sf::st_bbox(hex_sf)
  buffer <- calculate_view_buffer(hex_bbox)

  # Determine plot limits
  if (is.null(xlim)) {
    xlim <- c(hex_bbox["xmin"] - buffer$x, hex_bbox["xmax"] + buffer$x)
  }
  if (is.null(ylim)) {
    ylim <- c(hex_bbox["ymin"] - buffer$y, hex_bbox["ymax"] + buffer$y)
  }

  # Handle basemap (supports raster and vector)
  basemap_resolved <- resolve_basemap_with_raster(basemap)
  basemap_sf <- basemap_resolved$sf
  basemap_raster <- basemap_resolved$raster

  # Calculate aspect ratio for latitude
  mean_lat <- mean(ylim)
  asp <- 1 / cos(mean_lat * pi / 180)

  # Set up plot
  plot(xlim, ylim, type = "n",
       xlim = xlim, ylim = ylim,
       xlab = "Longitude", ylab = "Latitude",
       asp = asp, main = main, ...)

  # Draw basemap first (underneath hexagons)
  if (!is.null(basemap_raster)) {
    # Raster basemap
    if (inherits(basemap_raster, "SpatRaster")) {
      terra::plot(basemap_raster, add = TRUE)
    } else {
      raster::plot(basemap_raster, add = TRUE)
    }
  }

  if (!is.null(basemap_sf)) {
    # Clip basemap to plot extent for efficiency
    # Strip names from xlim/ylim to avoid NA issues in st_bbox
    clip_box <- sf::st_bbox(c(xmin = unname(xlim[1]),
                               xmax = unname(xlim[2]),
                               ymin = unname(ylim[1]),
                               ymax = unname(ylim[2])),
                             crs = sf::st_crs(4326))
    clip_poly <- sf::st_as_sfc(clip_box)

    # Use planar geometry for clipping (avoids S2 edge crossing errors)
    s2_was_used <- sf::sf_use_s2()
    sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(s2_was_used), add = TRUE)

    # Suppress warnings about attribute variables
    basemap_clipped <- suppressWarnings(
      sf::st_intersection(sf::st_make_valid(basemap_sf), clip_poly)
    )

    if (nrow(basemap_clipped) > 0) {
      plot(sf::st_geometry(basemap_clipped),
           col = basemap_fill, border = basemap_border,
           lwd = basemap_lwd,
           add = TRUE)
    }
  }

  # Draw hexagons on top
  # Apply alpha to fill color
  fill_col <- adjustcolor(fill, alpha.f = alpha)

  plot(sf::st_geometry(hex_sf),
       col = fill_col, border = border,
       lwd = lwd,
       add = TRUE)

  invisible(NULL)
}


#' Quick world map plot
#'
#' Simple wrapper to plot the built-in world map.
#'
#' @param fill Fill color for countries
#' @param border Border color for countries
#' @param ... Additional arguments passed to plot()
#'
#' @return NULL invisibly. Creates a plot as side effect.
#'
#' @family visualization
#' @export
#' @examples
#' \dontrun{
#' # Quick world map
#' plot_world()
#'
#' # Custom colors
#' plot_world(fill = "lightblue", border = "darkblue")
#' }
plot_world <- function(fill = "gray95", border = "gray50", ...) {

  plot(sf::st_geometry(hexify_world),
       col = fill, border = border, ...)
  invisible(NULL)
}


#' Create a heatmap visualization of hexagonal grid cells
#'
#' Creates a ggplot2-based heatmap of hexagonal grid cells colored by a value
#' column. Supports continuous and discrete color scales, projection
#' transformation, and customizable styling.
#'
#' @param data Data frame from hexify() containing cell_id and cell_area columns,
#'   or an sf object with hexagon polygons. Must include a column for coloring.
#' @param value Column name (as string) to use for fill color. If NULL and data
#'   has a 'count' or 'n' column, that will be used.
#' @param basemap Optional basemap. Can be:
#'   \itemize{
#'     \item \code{NULL}: No basemap (default)
#'     \item \code{"world"}: Use built-in \code{hexify_world} map (low resolution)
#'     \item \code{"world_hires"}: Use high-resolution map from rnaturalearth (requires package)
#'     \item An sf object: User-supplied vector map
#'   }
#' @param crs Target CRS for the map projection. Can be:
#'   \itemize{
#'     \item A numeric EPSG code (e.g., 4326 for WGS84, 3035 for LAEA Europe)
#'     \item A proj4 string
#'     \item An sf crs object
#'     \item NULL to use WGS84 (EPSG:4326)
#'   }
#' @param colors Color palette for the heatmap. Can be:
#'   \itemize{
#'     \item A character vector of colors (for manual scale)
#'     \item A single RColorBrewer palette name (e.g., "YlOrRd", "Greens")
#'     \item NULL to use viridis
#'   }
#' @param breaks Numeric vector of break points for binning continuous values,
#'   or NULL for continuous scale. Use \code{Inf} and \code{-Inf} for open-ended bins.
#' @param labels Labels for the breaks (length should be one less than breaks).
#'   If NULL, labels are auto-generated.
#' @param hex_border Border color for hexagons
#' @param hex_lwd Line width for hexagon borders
#' @param hex_alpha Transparency for hexagon fill (0-1)
#' @param basemap_fill Fill color for basemap polygons
#' @param basemap_border Border color for basemap polygons
#' @param basemap_lwd Line width for basemap borders
#' @param mask_outside Logical. If TRUE and basemap is provided, mask hexagon
#'   portions that fall outside the basemap polygons.
#' @param aperture Grid aperture (default 3), used if data is from hexify()
#' @param xlim Optional x-axis limits (in target CRS units) as c(min, max)
#' @param ylim Optional y-axis limits (in target CRS units) as c(min, max)
#' @param title Plot title
#' @param legend_title Title for the color legend
#' @param na_color Color for NA values
#' @param theme_void Logical. If TRUE (default), use a minimal theme without
#'   axes, gridlines, or background.
#'
#' @return A ggplot2 object that can be further customized or saved.
#'
#' @details
#' This function provides publication-quality heatmap visualizations of
#' hexagonal grids using ggplot2. Unlike \code{\link{hexify_map}}, it returns a
#' ggplot object that can be further customized with standard ggplot2 functions.
#'
#' @section Color Scales:
#' The function supports three types of color scales:
#' \describe{
#'   \item{Continuous}{Set \code{breaks = NULL} for a continuous gradient}
#'   \item{Binned}{Provide \code{breaks} vector to bin values into categories}
#'   \item{Discrete}{If \code{value} column is a factor, discrete colors are used}
#' }
#'
#' @section Projections:
#' Common projections:
#' \describe{
#'   \item{4326}{WGS84 (unprojected lat/lon)}
#'   \item{3035}{LAEA Europe}
#'   \item{3857}{Web Mercator}
#'   \item{"+proj=robin"}{Robinson (world maps)}
#'   \item{"+proj=moll"}{Mollweide (equal-area world maps)}
#' }
#'
#' @family visualization
#' @seealso \code{\link{hexify_map}} for base R plotting,
#'   \code{\link{hexify_plot}} for simple plots
#' @export
#' @examples
#' \dontrun{
#' library(hexify)
#' library(ggplot2)
#'
#' # Sample data with counts
#' cities <- data.frame(
#'   lon = c(16.37, 2.35, -3.70, 12.5, 4.9),
#'   lat = c(48.21, 48.86, 40.42, 41.9, 52.4),
#'   count = c(100, 250, 75, 180, 300)
#' )
#' result <- hexify(cities, lon = "lon", lat = "lat", area = 5000)
#'
#' # Simple heatmap
#' hexify_heatmap(result, value = "count")
#'
#' # With world basemap and custom colors
#' hexify_heatmap(result, value = "count",
#'                basemap = "world",
#'                colors = "YlOrRd",
#'                title = "City Density")
#'
#' # Binned values with custom breaks
#' hexify_heatmap(result, value = "count",
#'                basemap = "world",
#'                breaks = c(-Inf, 100, 200, Inf),
#'                labels = c("Low", "Medium", "High"),
#'                colors = c("#fee8c8", "#fdbb84", "#e34a33"))
#'
#' # Different projection (LAEA Europe)
#' hexify_heatmap(result, value = "count",
#'                basemap = "world",
#'                crs = 3035,
#'                xlim = c(2500000, 6500000),
#'                ylim = c(1500000, 5500000))
#'
#' # Customize further with ggplot2
#' hexify_heatmap(result, value = "count", basemap = "world") +
#'   labs(caption = "Data source: Example") +
#'   theme(legend.position = "bottom")
#' }
hexify_heatmap <- function(data,
                        value = NULL,
                        basemap = NULL,
                        crs = NULL,
                        colors = NULL,
                        breaks = NULL,
                        labels = NULL,
                        hex_border = "gray30",
                        hex_lwd = 0.2,
                        hex_alpha = 0.7,
                        basemap_fill = "white",
                        basemap_border = "black",
                        basemap_lwd = 0.5,
                        mask_outside = FALSE,
                        aperture = 3L,
                        xlim = NULL,
                        ylim = NULL,
                        title = NULL,
                        legend_title = NULL,
                        na_color = "gray90",
                        theme_void = TRUE) {

  # Check ggplot2 availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for hexify_heatmap(). ",
         "Install it with: install.packages('ggplot2')")
  }

  # Prepare hex sf with extra columns merged
  hex_sf <- prepare_hex_sf(data, aperture)

  # Resolve value column
  value <- resolve_value_column(hex_sf, value)

  # Setup CRS
  crs <- if (is.null(crs)) 4326 else crs
  if (is.na(sf::st_crs(hex_sf))) sf::st_crs(hex_sf) <- 4326
  hex_sf <- sf::st_transform(hex_sf, crs)

  # Resolve and transform basemap
  basemap_sf <- resolve_basemap(basemap)
  if (!is.null(basemap_sf)) {
    basemap_sf <- sf::st_transform(basemap_sf, crs)
  }

  # Disable S2 for geometry operations (avoids edge crossing errors)
  s2_was_used <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_used), add = TRUE)

  # Create mask for areas outside basemap (if requested)
  mask_sf <- NULL
  if (mask_outside && !is.null(basemap_sf)) {
    mask_sf <- create_outside_mask(hex_sf, basemap_sf, xlim, ylim)
  }

  # Prepare fill column (apply breaks if needed)
  fill_info <- prepare_fill_column(hex_sf, value, breaks, labels)
  hex_sf <- fill_info$data
  fill_col <- fill_info$fill_col
  is_discrete <- fill_info$is_discrete

  # Set legend title
  legend_title <- if (is.null(legend_title)) value else legend_title

  # Build ggplot with layers
  p <- ggplot2::ggplot()
  if (mask_outside && !is.null(basemap_sf)) {
    p <- build_masked_layers(
      p, hex_sf, fill_col, hex_border, hex_lwd, hex_alpha,
      basemap_sf, basemap_border, basemap_lwd, mask_sf
    )
  } else {
    p <- build_standard_layers(
      p, hex_sf, fill_col, hex_border, hex_lwd, hex_alpha,
      basemap_sf, basemap_fill, basemap_border, basemap_lwd
    )
  }

  # Apply color scale using helper functions
  n_levels <- length(unique(hex_sf[[fill_col]]))
  if (is_discrete) {
    p <- apply_discrete_scale(p, colors, legend_title, na_color, n_levels)
  } else {
    p <- apply_continuous_scale(p, colors, legend_title, na_color)
  }

  # Set coordinate system with limits
  p <- if (!is.null(xlim) || !is.null(ylim)) {
    p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    p + ggplot2::coord_sf(expand = FALSE)
  }

  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  # Apply theme
  if (theme_void) {
    p <- p + ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )
  }

  p
}

