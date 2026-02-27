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
  # Handle HexData objects
  if (is_hex_data(data)) {
    g <- data@grid
    underlying_data <- data@data
    unique_cells <- unique(data@cell_id)
    hex_sf <- cell_to_sf(unique_cells, g)

    # Join extra columns from original data
    extra_cols <- setdiff(names(underlying_data), c("cell_id", "geometry"))
    if (length(extra_cols) > 0) {
      # Build data frame with cell_id for merge
      data_with_id <- cbind(underlying_data, cell_id = data@cell_id)
      cols <- c("cell_id", extra_cols)
      data_unique <- data_with_id[!duplicated(data_with_id$cell_id), cols, drop = FALSE]
      hex_sf <- merge(hex_sf, data_unique, by = "cell_id", all.x = TRUE)
    }
    return(hex_sf)
  }

  if (inherits(data, "sf")) return(data)

  if (!is.data.frame(data) || !"cell_id" %in% names(data)) {
    stop("data must be a HexData object or an sf object")
  }

  # Handle legacy data frames with cell_area column
  if ("cell_area" %in% names(data) || "cell_area_km2" %in% names(data)) {
    # Get area to determine resolution
    area <- if ("cell_area_km2" %in% names(data)) data$cell_area_km2[1]
            else if ("cell_area" %in% names(data)) data$cell_area[1]

    # Create temporary grid
    grid <- hex_grid(area_km2 = area, aperture = aperture)

    # Generate polygons
    unique_cells <- unique(data$cell_id)
    hex_sf <- cell_to_sf(unique_cells, grid)

    # Join extra columns from original data
    extra_cols <- setdiff(names(data), c("cell_id", "geometry"))
    if (length(extra_cols) > 0) {
      cols <- c("cell_id", extra_cols)
      data_unique <- data[!duplicated(data$cell_id), cols, drop = FALSE]
      hex_sf <- merge(hex_sf, data_unique, by = "cell_id", all.x = TRUE)
    }
    return(hex_sf)
  }

  stop("data must contain 'cell_area' or 'cell_area_km2' column (output from hexify()).")
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
                                basemap_lwd, mask_sf, uniform_fill = "#E69F00") {
 # Hexagons first
  if (is.null(fill_col)) {
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      fill = uniform_fill,
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )
  } else {
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      ggplot2::aes(fill = .data[[fill_col]]),
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )
  }

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
                                  basemap_border, basemap_lwd,
                                  uniform_fill = "#E69F00") {
  if (!is.null(basemap_sf)) {
    p <- p + ggplot2::geom_sf(
      data = basemap_sf,
      fill = basemap_fill,
      color = basemap_border,
      linewidth = basemap_lwd
    )
  }

  if (is.null(fill_col)) {
    p + ggplot2::geom_sf(
      data = hex_sf,
      fill = uniform_fill,
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )
  } else {
    p + ggplot2::geom_sf(
      data = hex_sf,
      ggplot2::aes(fill = .data[[fill_col]]),
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )
  }
}

#' Simple sf preparation for hexify_map (no extra column merging)
#' @noRd
prepare_hex_sf_simple <- function(data, aperture) {
  # Handle HexData objects
  if (is_hex_data(data)) {
    g <- data@grid
    unique_cells <- unique(data@cell_id)
    return(cell_to_sf(unique_cells, g))
  }

  if (inherits(data, "sf")) return(data)

  if (!is.data.frame(data) || !"cell_id" %in% names(data)) {
    stop("data must be a HexData object or an sf object")
  }

  # Handle legacy data frames with cell_area column
  if ("cell_area" %in% names(data) || "cell_area_km2" %in% names(data)) {
    # Get area to determine resolution
    area <- if ("cell_area_km2" %in% names(data)) data$cell_area_km2[1]
            else if ("cell_area" %in% names(data)) data$cell_area[1]

    # Create temporary grid
    grid <- hex_grid(area_km2 = area, aperture = aperture)

    # Generate polygons
    unique_cells <- unique(data$cell_id)
    return(cell_to_sf(unique_cells, grid))
  }

  stop("data must contain 'cell_area' or 'cell_area_km2' column (output from hexify()).")
}

#' Resolve value column name (auto-detect if NULL)
#' @noRd
resolve_value_column <- function(hex_sf, value, require = FALSE) {
  if (!is.null(value)) {
    if (!value %in% names(hex_sf)) {
      stop("Column '", value, "' not found in data. ",
           "Available columns: ", paste(names(hex_sf), collapse = ", "))
    }
    return(value)
  }

  # Auto-detect common value columns

  if ("count" %in% names(hex_sf)) return("count")
  if ("n" %in% names(hex_sf)) return("n")


  # No value column found - return NULL for uniform fill
  if (require) {
    stop("No 'value' column specified and no 'count' or 'n' column found in data")
  }
  NULL
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

stop("basemap must be 'world', an sf object, or a SpatRaster (terra)")
}

# =============================================================================
# PUBLIC FUNCTIONS
# =============================================================================

# hexify_map() removed - use plot() method on HexData objects instead


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
#' # Quick world map
#' plot_world()
#'
#' # Custom colors
#' plot_world(fill = "lightblue", border = "darkblue")
plot_world <- function(fill = "gray90", border = "gray50", ...) {

  plot(sf::st_geometry(hexify_world),
       col = fill, border = border, ...)
  invisible(NULL)
}


#' Create a ggplot2 visualization of hexagonal grid cells
#'
#' Creates a ggplot2-based visualization of hexagonal grid cells, optionally
#' colored by a value column. Supports continuous and discrete color scales,
#' projection transformation, and customizable styling.
#'
#' @param data A HexData object from \code{hexify()}, a data frame with cell_id
#'   and cell_area columns, or an sf object with hexagon polygons.
#' @param value Column name (as string) to use for fill color. If NULL, cells
#'   are drawn with a uniform fill color. If not specified but data has a
#'   'count' or 'n' column, that will be used automatically.
#' @param basemap Optional basemap. Can be:
#'   \itemize{
#'     \item \code{NULL}: No basemap (default)
#'     \item \code{"world"}: Use built-in \code{hexify_world} map (low resolution)
#'     \item \code{"world_hires"}: Use high-resolution map from rnaturalearth (requires package)
#'     \item An sf object: User-supplied vector map
#'   }
#' @param crs Target CRS for the map projection. Can be:
#'   \itemize{
#'     \item A numeric EPSG code (e.g., 4326 for 'WGS84', 3035 for LAEA Europe)
#'     \item A proj4 string
#'     \item An sf crs object
#'     \item NULL to use 'WGS84' (EPSG:4326)
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
#' hexagonal grids using ggplot2. It returns a ggplot object that can be
#' further customized with standard ggplot2 functions.
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
#'   \item{4326}{'WGS84' (unprojected lat/lon)}
#'   \item{3035}{LAEA Europe}
#'   \item{3857}{Web Mercator}
#'   \item{"+proj=robin"}{Robinson (world maps)}
#'   \item{"+proj=moll"}{Mollweide (equal-area world maps)}
#' }
#'
#' @family visualization
#' @seealso \code{\link{plot_grid}} for base R plotting,
#'   \code{\link{cell_to_sf}} to generate polygons manually
#' @export
#' @examples
#' library(hexify)
#'
#' # Sample data with counts
#' cities <- data.frame(
#'   lon = c(16.37, 2.35, -3.70, 12.5, 4.9),
#'   lat = c(48.21, 48.86, 40.42, 41.9, 52.4),
#'   count = c(100, 250, 75, 180, 300)
#' )
#' result <- hexify(cities, lon = "lon", lat = "lat", area_km2 = 5000)
#'
#' # Simple plot (uniform fill, no value mapping)
#' hexify_heatmap(result)
#'
#' \donttest{
#' library(ggplot2)
#'
#' # With world basemap
#' hexify_heatmap(result, basemap = "world")
#'
#' # Heatmap with value mapping
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
                        hex_border = "#5D4E37",
                        hex_lwd = 0.3,
                        hex_alpha = 0.7,
                        basemap_fill = "gray90",
                        basemap_border = "gray50",
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

  # Resolve value column (NULL means uniform fill)
  value <- resolve_value_column(hex_sf, value, require = FALSE)

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

  # Prepare fill column (apply breaks if needed) - skip if no value
 fill_col <- NULL
  is_discrete <- FALSE
  if (!is.null(value)) {
    fill_info <- prepare_fill_column(hex_sf, value, breaks, labels)
    hex_sf <- fill_info$data
    fill_col <- fill_info$fill_col
    is_discrete <- fill_info$is_discrete
  }

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

  # Apply color scale (only if we have a value column)
  if (!is.null(fill_col)) {
    n_levels <- length(unique(hex_sf[[fill_col]]))
    if (is_discrete) {
      p <- apply_discrete_scale(p, colors, legend_title, na_color, n_levels)
    } else {
      p <- apply_continuous_scale(p, colors, legend_title, na_color)
    }
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

