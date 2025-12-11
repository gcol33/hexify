# hex_map.R
# Visualization functions with basemap support

#' Plot hexagonal grid cells with optional basemap
#'
#' Creates a map visualization of hexagonal grid cells. Supports the built-in
#' world map or user-supplied basemaps (sf vectors or raster images).
#'
#' @param data Data frame from hexify() containing hex_id and hex_area columns,
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
#' \code{\link{hex_to_polygons}} to get an sf object and plot with ggplot2,
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
#' hex_map(result)
#'
#' # With built-in world map
#' hex_map(result, basemap = "world")
#'
#' # Custom colors
#' hex_map(result, basemap = "world",
#'         fill = "steelblue", border = "darkblue",
#'         basemap_fill = "ivory", basemap_border = "gray50")
#'
#' # With user-supplied sf basemap
#' library(rnaturalearth)
#' europe <- ne_countries(continent = "Europe", returnclass = "sf")
#' hex_map(result, basemap = europe)
#'
#' # Zoom to specific region
#' hex_map(result, basemap = "world",
#'         xlim = c(-10, 25), ylim = c(35, 55))
#' }
hex_map <- function(data,
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
 if (inherits(data, "sf")) {
    hex_sf <- data
  } else if (is.data.frame(data) && "hex_id" %in% names(data)) {
    if (!"hex_area" %in% names(data)) {
      stop("data must contain 'hex_area' column (output from hexify()). ",
           "Use hex_polygons() directly if you have pre-computed polygons.")
    }
    hex_sf <- hex_to_polygons(data, aperture = aperture, return_sf = TRUE)
  } else {
    stop("data must be a data.frame from hexify() or an sf object with polygons")
  }

  # Get bounding box for the hexagons
  hex_bbox <- sf::st_bbox(hex_sf)

  # Add buffer around data extent (10% of range)
  x_range <- hex_bbox["xmax"] - hex_bbox["xmin"]
  y_range <- hex_bbox["ymax"] - hex_bbox["ymin"]
  buffer_x <- max(x_range * 0.1, 1)  # At least 1 degree

  buffer_y <- max(y_range * 0.1, 1)

  # Determine plot limits
  if (is.null(xlim)) {
    xlim <- c(hex_bbox["xmin"] - buffer_x, hex_bbox["xmax"] + buffer_x)
  }
  if (is.null(ylim)) {
    ylim <- c(hex_bbox["ymin"] - buffer_y, hex_bbox["ymax"] + buffer_y)
  }

  # Handle basemap
  basemap_sf <- NULL
  basemap_raster <- NULL

  if (!is.null(basemap)) {
    if (is.character(basemap) && basemap == "world") {
      # Use built-in world map
      basemap_sf <- hexify_world
    } else if (inherits(basemap, "sf") || inherits(basemap, "sfc")) {
      basemap_sf <- basemap
    } else if (inherits(basemap, "SpatRaster")) {
      # terra raster
      if (!requireNamespace("terra", quietly = TRUE)) {
        stop("Package 'terra' is required for SpatRaster basemaps")
      }
      basemap_raster <- basemap
    } else if (inherits(basemap, c("RasterLayer", "RasterBrick", "RasterStack"))) {
      # raster package
      if (!requireNamespace("raster", quietly = TRUE)) {
        stop("Package 'raster' is required for Raster* basemaps")
      }
      basemap_raster <- basemap
    } else {
      stop("basemap must be 'world', an sf object, or a raster (terra/raster)")
    }
  }

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
#' @param data Data frame from hexify() containing hex_id and hex_area columns,
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
#' hexagonal grids using ggplot2. Unlike \code{\link{hex_map}}, it returns a
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
#' hex_heatmap(result, value = "count")
#'
#' # With world basemap and custom colors
#' hex_heatmap(result, value = "count",
#'             basemap = "world",
#'             colors = "YlOrRd",
#'             title = "City Density")
#'
#' # Binned values with custom breaks
#' hex_heatmap(result, value = "count",
#'             basemap = "world",
#'             breaks = c(-Inf, 100, 200, Inf),
#'             labels = c("Low", "Medium", "High"),
#'             colors = c("#fee8c8", "#fdbb84", "#e34a33"))
#'
#' # Different projection (LAEA Europe)
#' hex_heatmap(result, value = "count",
#'             basemap = "world",
#'             crs = 3035,
#'             xlim = c(2500000, 6500000),
#'             ylim = c(1500000, 5500000))
#'
#' # Customize further with ggplot2
#' hex_heatmap(result, value = "count", basemap = "world") +
#'   labs(caption = "Data source: Example") +
#'   theme(legend.position = "bottom")
#' }
hex_heatmap <- function(data,
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
    stop("Package 'ggplot2' is required for hex_heatmap(). ",
         "Install it with: install.packages('ggplot2')")
  }

  # Convert hexify output to polygons if needed
  if (inherits(data, "sf")) {
    hex_sf <- data
  } else if (is.data.frame(data) && "hex_id" %in% names(data)) {
    if (!"hex_area" %in% names(data)) {
      stop("data must contain 'hex_area' column (output from hexify()). ",
           "Use hex_polygons() directly if you have pre-computed polygons.")
    }
    # Generate polygons and merge with original data
    hex_sf <- hex_to_polygons(data, aperture = aperture, return_sf = TRUE)
    # Join original data columns (excluding duplicates)
    extra_cols <- setdiff(names(data), c("hex_id", "geometry"))
    if (length(extra_cols) > 0) {
      # Aggregate by hex_id for cases where multiple rows share same hex_id
      data_unique <- data[!duplicated(data$hex_id), c("hex_id", extra_cols), drop = FALSE]
      hex_sf <- merge(hex_sf, data_unique, by = "hex_id", all.x = TRUE)
    }
  } else {
    stop("data must be a data.frame from hexify() or an sf object with polygons")
  }

  # Auto-detect value column if not specified
  if (is.null(value)) {
    if ("count" %in% names(hex_sf)) {
      value <- "count"
    } else if ("n" %in% names(hex_sf)) {
      value <- "n"
    } else {
      stop("No 'value' column specified and no 'count' or 'n' column found in data")
    }
  }

  # Validate value column exists
  if (!value %in% names(hex_sf)) {
    stop("Column '", value, "' not found in data. ",
         "Available columns: ", paste(names(hex_sf), collapse = ", "))
  }

  # Set default CRS (WGS84)
  if (is.null(crs)) {
    crs <- 4326
  }

  # Ensure hex_sf has CRS
  if (is.na(sf::st_crs(hex_sf))) {
    sf::st_crs(hex_sf) <- 4326
  }

  # Transform to target CRS
  hex_sf <- sf::st_transform(hex_sf, crs)

  # Handle basemap
  basemap_sf <- NULL
  if (!is.null(basemap)) {
    if (is.character(basemap) && basemap == "world") {
      basemap_sf <- hexify_world
    } else if (is.character(basemap) && basemap == "world_hires") {
      # Use high-resolution map from rnaturalearth
      if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
        stop("Package 'rnaturalearth' is required for basemap = 'world_hires'. ",
             "Install with: install.packages('rnaturalearth')")
      }
      basemap_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    } else if (inherits(basemap, "sf") || inherits(basemap, "sfc")) {
      basemap_sf <- basemap
    } else {
      stop("basemap must be 'world', 'world_hires', or an sf object")
    }
    # Transform basemap to target CRS
    basemap_sf <- sf::st_transform(basemap_sf, crs)
  }

  # Disable S2 for geometry operations (avoids edge crossing errors)
  s2_was_used <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_used), add = TRUE)

  # Create mask for areas outside basemap (if requested)
  mask_sf <- NULL
  if (mask_outside && !is.null(basemap_sf)) {
    # Use xlim/ylim if provided, otherwise use hex bounding box
    if (!is.null(xlim) && !is.null(ylim)) {
      bbox_coords <- c(
        xmin = xlim[1],
        xmax = xlim[2],
        ymin = ylim[1],
        ymax = ylim[2]
      )
    } else {
      # Get bounding box of hex data with buffer
      hex_bbox <- sf::st_bbox(hex_sf)
      x_range <- hex_bbox["xmax"] - hex_bbox["xmin"]
      y_range <- hex_bbox["ymax"] - hex_bbox["ymin"]
      buffer_x <- max(x_range * 0.1, 1)
      buffer_y <- max(y_range * 0.1, 1)
      bbox_coords <- c(
        xmin = unname(hex_bbox["xmin"] - buffer_x),
        xmax = unname(hex_bbox["xmax"] + buffer_x),
        ymin = unname(hex_bbox["ymin"] - buffer_y),
        ymax = unname(hex_bbox["ymax"] + buffer_y)
      )
    }

    # Create bbox polygon in the target CRS
    bbox_poly <- sf::st_as_sfc(sf::st_bbox(bbox_coords, crs = sf::st_crs(hex_sf)))

    # Create inverse mask (areas in bbox NOT covered by basemap)
    basemap_union <- sf::st_union(sf::st_make_valid(basemap_sf))
    mask_sf <- suppressWarnings(sf::st_difference(bbox_poly, basemap_union))
  }

  # Determine if binned or continuous scale
  value_data <- hex_sf[[value]]
  is_discrete <- is.factor(value_data) || is.character(value_data)

  # Create binned factor if breaks provided
  if (!is.null(breaks) && !is_discrete) {
    if (is.null(labels)) {
      # Auto-generate labels
      n_bins <- length(breaks) - 1
      labels <- character(n_bins)
      for (i in seq_len(n_bins)) {
        low <- breaks[i]
        high <- breaks[i + 1]
        if (is.infinite(low) && low < 0) {
          labels[i] <- paste0("<", high)
        } else if (is.infinite(high)) {
          labels[i] <- paste0(">", low)
        } else {
          labels[i] <- paste0(low, "-", high)
        }
      }
    }
    hex_sf[[paste0(value, "_bin")]] <- cut(
      value_data,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE
    )
    fill_col <- paste0(value, "_bin")
    is_discrete <- TRUE
  } else {
    fill_col <- value
  }

  # Set legend title
  if (is.null(legend_title)) {
    legend_title <- value
  }

  # Build ggplot
  p <- ggplot2::ggplot()

  # Layer order depends on mask_outside setting
 if (mask_outside && !is.null(basemap_sf)) {
    # When masking: hex first, then mask, then borders on top
    # Add hexagons first
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      ggplot2::aes(fill = .data[[fill_col]]),
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )

    # Add mask layer to hide hexes outside land
    if (!is.null(mask_sf)) {
      p <- p + ggplot2::geom_sf(
        data = mask_sf,
        fill = "white",
        color = NA
      )
    }

    # Add basemap borders on top (no fill)
    p <- p + ggplot2::geom_sf(
      data = basemap_sf,
      fill = NA,
      color = basemap_border,
      linewidth = basemap_lwd
    )
  } else {
    # Standard mode: basemap underneath, hexes on top
    if (!is.null(basemap_sf)) {
      p <- p + ggplot2::geom_sf(
        data = basemap_sf,
        fill = basemap_fill,
        color = basemap_border,
        linewidth = basemap_lwd
      )
    }

    # Add hexagons on top
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      ggplot2::aes(fill = .data[[fill_col]]),
      color = hex_border,
      linewidth = hex_lwd,
      alpha = hex_alpha
    )
  }

  # Check if colors is a palette name (single string that's not a hex color)
  is_palette_name <- function(x) {
    length(x) == 1 && is.character(x) && !grepl("^#", x) && !x %in% grDevices::colors()
  }

  # Add color scale
  if (is_discrete) {
    if (!is.null(colors)) {
      if (is_palette_name(colors)) {
        # Try RColorBrewer palette
        if (requireNamespace("RColorBrewer", quietly = TRUE) &&
            colors %in% row.names(RColorBrewer::brewer.pal.info)) {
          n_colors <- length(unique(hex_sf[[fill_col]]))
          pal_colors <- RColorBrewer::brewer.pal(
            min(n_colors, RColorBrewer::brewer.pal.info[colors, "maxcolors"]),
            colors
          )
          p <- p + ggplot2::scale_fill_manual(
            values = pal_colors,
            name = legend_title,
            na.value = na_color
          )
        } else {
          # Fallback: treat as viridis option name
          p <- p + ggplot2::scale_fill_viridis_d(
            option = tolower(colors),
            name = legend_title,
            na.value = na_color
          )
        }
      } else {
        # Manual color vector
        p <- p + ggplot2::scale_fill_manual(
          values = colors,
          name = legend_title,
          na.value = na_color
        )
      }
    } else {
      # Default viridis discrete
      p <- p + ggplot2::scale_fill_viridis_d(
        name = legend_title,
        na.value = na_color
      )
    }
  } else {
    # Continuous scale
    if (!is.null(colors)) {
      if (is_palette_name(colors)) {
        # Try RColorBrewer palette via scale_fill_distiller
        if (requireNamespace("RColorBrewer", quietly = TRUE) &&
            colors %in% row.names(RColorBrewer::brewer.pal.info)) {
          p <- p + ggplot2::scale_fill_distiller(
            palette = colors,
            direction = 1,
            name = legend_title,
            na.value = na_color
          )
        } else {
          # Fallback: treat as viridis option name
          p <- p + ggplot2::scale_fill_viridis_c(
            option = tolower(colors),
            name = legend_title,
            na.value = na_color
          )
        }
      } else {
        # Manual gradient
        p <- p + ggplot2::scale_fill_gradientn(
          colors = colors,
          name = legend_title,
          na.value = na_color
        )
      }
    } else {
      # Default viridis continuous
      p <- p + ggplot2::scale_fill_viridis_c(
        name = legend_title,
        na.value = na_color
      )
    }
  }

  # Set coordinate system with limits
  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    )
  } else {
    p <- p + ggplot2::coord_sf(expand = FALSE)
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

  return(p)
}
