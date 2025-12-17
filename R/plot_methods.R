# plot_methods.R
# Plot methods for HexData objects
#
# Provides default plotting with sensible styling options.
# Users can visualize hexified data without specifying grid parameters.

# =============================================================================
# BASE R PLOT METHOD
# =============================================================================

#' Plot HexData objects
#'
#' Default plot method for HexData objects. Draws hexagonal grid cells
#' with an optional basemap.
#'
#' @param x A HexData object from \code{hexify()}
#' @param y Ignored (for S4 method compatibility)
#' @param basemap Basemap specification:
#'   \itemize{
#'     \item \code{TRUE} or \code{"world"}: Use built-in world map
#'     \item \code{FALSE} or \code{NULL}: No basemap
#'     \item sf object: Custom basemap
#'   }
#' @param clip_basemap Clip basemap to data extent (default TRUE). Clipping
#'   temporarily disables S2 spherical geometry to avoid edge-crossing errors.
#' @param basemap_fill Fill color for basemap (default "gray90")
#' @param basemap_border Border color for basemap (default "gray50")
#' @param basemap_lwd Line width for basemap borders (default 0.5)
#' @param grid_fill Fill color for grid cells (default "#E69F00" - amber/orange)
#' @param grid_border Border color for grid cells (default "#5D4E37" - dark brown)
#' @param grid_lwd Line width for cell borders (default 0.8)
#' @param grid_alpha Transparency for cell fill (0-1, default 0.7)
#' @param fill Column name for fill mapping (optional)
#' @param show_points Show original points on top of cells (default FALSE)
#' @param point_size Size of points if shown (default 1)
#' @param point_color Color of points (default "red")
#' @param crop Crop to data extent (default TRUE)
#' @param crop_expand Expansion factor for crop (default 0.1)
#' @param main Plot title
#' @param ... Additional arguments passed to base plot()
#'
#' @return Invisibly returns the HexData object
#'
#' @details
#' This function generates polygon geometries for the cells present in
#' the data and plots them. Polygons are computed on demand, not stored,
#' to minimize memory usage.
#'
#' @seealso \code{\link{hexify_ggplot}} for ggplot2 plotting
#'
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(lon = runif(100, -10, 10), lat = runif(100, 40, 50))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Basic plot
#' plot(result)
#'
#' # With world basemap
#' plot(result, basemap = TRUE)
#'
#' # Custom styling
#' plot(result, basemap = TRUE,
#'      grid_fill = "lightblue", grid_border = "darkblue",
#'      basemap_fill = "ivory")
#'
#' # Show original points
#' plot(result, basemap = TRUE, show_points = TRUE)
#' }
setMethod("plot", signature(x = "HexData", y = "missing"),
  function(x, y,
           basemap = TRUE,
           clip_basemap = TRUE,
           basemap_fill = "gray90",
           basemap_border = "gray50",
           basemap_lwd = 0.5,
           grid_fill = "#E69F00",
           grid_border = "#5D4E37",
           grid_lwd = 0.8,
           grid_alpha = 0.7,
           fill = NULL,
           show_points = FALSE,
           point_size = 1,
           point_color = "red",
           crop = TRUE,
           crop_expand = 0.1,
           main = NULL,
           ...) {

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required for plotting")
    }

    # Get grid spec
    g <- x@grid

    # Generate cell polygons
    unique_cells <- unique(x@cell_id)
    hex_sf <- cell_to_sf(unique_cells, g)

    # If fill column specified, merge data
    if (!is.null(fill)) {
      if (!fill %in% names(x@data)) {
        stop(sprintf("Column '%s' not found in data", fill))
      }
      # Create data frame with cell_id for merging
      data_with_id <- cbind(x@data, cell_id = x@cell_id)
      # Aggregate by cell if needed (take first value)
      agg_data <- data_with_id[!duplicated(data_with_id$cell_id), c("cell_id", fill)]
      hex_sf <- merge(hex_sf, agg_data, by = "cell_id", all.x = TRUE)
    }

    # Calculate bounding box
    hex_bbox <- sf::st_bbox(hex_sf)

    if (crop) {
      x_range <- hex_bbox["xmax"] - hex_bbox["xmin"]
      y_range <- hex_bbox["ymax"] - hex_bbox["ymin"]
      xlim <- c(hex_bbox["xmin"] - x_range * crop_expand,
                hex_bbox["xmax"] + x_range * crop_expand)
      ylim <- c(hex_bbox["ymin"] - y_range * crop_expand,
                hex_bbox["ymax"] + y_range * crop_expand)
    } else {
      xlim <- c(-180, 180)
      ylim <- c(-90, 90)
    }

    # Calculate aspect ratio
    mean_lat <- mean(ylim)
    asp <- 1 / cos(mean_lat * pi / 180)

    # Resolve basemap
    basemap_sf <- NULL
    if (isTRUE(basemap) || identical(basemap, "world")) {
      basemap_sf <- hexify_world
    } else if (inherits(basemap, c("sf", "sfc"))) {
      basemap_sf <- basemap
    }

    # Initialize plot
    plot(xlim, ylim, type = "n",
         xlim = xlim, ylim = ylim,
         xlab = "Longitude", ylab = "Latitude",
         asp = asp, main = main, ...)

    # Draw basemap
    if (!is.null(basemap_sf)) {
      # Disable S2 for clipping
      s2_state <- sf::sf_use_s2()
      sf::sf_use_s2(FALSE)
      on.exit(sf::sf_use_s2(s2_state), add = TRUE)

      # Clip basemap to extent
      clip_box <- sf::st_bbox(c(xmin = unname(xlim[1]),
                                 xmax = unname(xlim[2]),
                                 ymin = unname(ylim[1]),
                                 ymax = unname(ylim[2])),
                               crs = 4326)
      clip_poly <- sf::st_as_sfc(clip_box)

      basemap_clipped <- suppressWarnings(
        sf::st_intersection(sf::st_make_valid(basemap_sf), clip_poly)
      )

      if (nrow(basemap_clipped) > 0) {
        plot(sf::st_geometry(basemap_clipped),
             col = basemap_fill, border = basemap_border,
             lwd = basemap_lwd, add = TRUE)
      }
    }

    # Draw hexagons
    if (is.null(fill)) {
      fill_colors <- adjustcolor(grid_fill, alpha.f = grid_alpha)
    } else {
      # Map fill column to colors
      values <- hex_sf[[fill]]
      if (is.numeric(values)) {
        # Continuous: use viridis-like palette
        n_colors <- 100
        pal <- grDevices::colorRampPalette(
          c("#440154", "#3B528B", "#21918C", "#5DC863", "#FDE725")
        )(n_colors)
        scaled <- (values - min(values, na.rm = TRUE)) /
                  (max(values, na.rm = TRUE) - min(values, na.rm = TRUE))
        scaled[is.na(scaled)] <- 0
        idx <- pmax(1, pmin(n_colors, ceiling(scaled * n_colors)))
        fill_colors <- adjustcolor(pal[idx], alpha.f = grid_alpha)
      } else {
        # Discrete: use basic palette
        levels <- unique(values)
        pal <- grDevices::rainbow(length(levels))
        fill_colors <- adjustcolor(
          pal[match(values, levels)],
          alpha.f = grid_alpha
        )
      }
    }

    plot(sf::st_geometry(hex_sf),
         col = fill_colors, border = grid_border,
         lwd = grid_lwd, add = TRUE)

    # Draw points if requested
    if (show_points) {
      # Try to get coordinates from sf geometry or cell centers
      if (inherits(x@data, "sf")) {
        coords <- sf::st_coordinates(x@data)
        points(coords[, 1], coords[, 2],
               pch = 19, cex = point_size, col = point_color)
      } else {
        # Use cell centers as fallback
        points(x@cell_center[, "lon"], x@cell_center[, "lat"],
               pch = 19, cex = point_size, col = point_color)
      }
    }

    invisible(x)
  }
)

# =============================================================================
# GGPLOT2 PLOTTING FUNCTION
# =============================================================================

#' Create a ggplot2 visualization of HexData
#'
#' Generates a ggplot2 object for HexData, supporting fill mapping,
#' basemaps, and advanced customization.
#'
#' @param x A HexData object from \code{hexify()}
#' @param basemap Basemap specification (see \code{plot.HexData})
#' @param basemap_fill Fill color for basemap
#' @param basemap_border Border color for basemap
#' @param basemap_lwd Line width for basemap
#' @param grid_border Border color for grid cells
#' @param grid_lwd Line width for cell borders
#' @param grid_alpha Transparency for cell fill
#' @param fill Column name for fill mapping (optional)
#' @param show_points Show original points
#' @param point_size Size of points
#' @param point_color Color of points
#' @param crop Crop to data extent
#' @param crop_expand Expansion factor for crop
#' @param title Plot title
#' @param legend_title Legend title
#' @param ... Additional arguments (ignored)
#'
#' @return A ggplot object that can be further customized
#'
#' @details
#' Requires ggplot2 package. The returned object can be modified using
#' standard ggplot2 functions like \code{+}, \code{theme()}, etc.
#'
#' @seealso \code{\link{plot,HexData,missing-method}} for base R plotting
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 1000)
#'
#' # Basic plot
#' hexify_ggplot(result)
#'
#' # With fill mapping
#' result$data$count <- sample(1:100, nrow(result$data))
#' hexify_ggplot(result, fill = "count") +
#'   scale_fill_viridis_c()
#'
#' # Customize
#' hexify_ggplot(result, basemap = TRUE) +
#'   theme_minimal() +
#'   labs(title = "Hexified Data")
#' }
hexify_ggplot <- function(x,
                              basemap = TRUE,
                              basemap_fill = "gray90",
                              basemap_border = "gray50",
                              basemap_lwd = 0.3,
                              grid_border = "#5D4E37",
                              grid_lwd = 0.4,
                              grid_alpha = 0.7,
                              fill = NULL,
                              show_points = FALSE,
                              point_size = 1,
                              point_color = "red",
                              crop = TRUE,
                              crop_expand = 0.1,
                              title = NULL,
                              legend_title = NULL,
                              ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for hexify_ggplot()")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for hexify_ggplot()")
  }

  g <- x@grid

  # Generate cell polygons
  unique_cells <- unique(x@cell_id)
  hex_sf <- cell_to_sf(unique_cells, g)

  # Merge fill column if specified
  if (!is.null(fill)) {
    if (!fill %in% names(x@data)) {
      stop(sprintf("Column '%s' not found in data", fill))
    }
    # Create data frame with cell_id for merging
    data_with_id <- cbind(x@data, cell_id = x@cell_id)
    agg_data <- data_with_id[!duplicated(data_with_id$cell_id), c("cell_id", fill)]
    hex_sf <- merge(hex_sf, agg_data, by = "cell_id", all.x = TRUE)
  }

  # Calculate limits
  hex_bbox <- sf::st_bbox(hex_sf)
  if (crop) {
    x_range <- hex_bbox["xmax"] - hex_bbox["xmin"]
    y_range <- hex_bbox["ymax"] - hex_bbox["ymin"]
    xlim <- c(hex_bbox["xmin"] - x_range * crop_expand,
              hex_bbox["xmax"] + x_range * crop_expand)
    ylim <- c(hex_bbox["ymin"] - y_range * crop_expand,
              hex_bbox["ymax"] + y_range * crop_expand)
  } else {
    xlim <- NULL
    ylim <- NULL
  }

  # Resolve basemap
  basemap_sf <- NULL
  if (isTRUE(basemap) || identical(basemap, "world")) {
    basemap_sf <- hexify_world
  } else if (inherits(basemap, c("sf", "sfc"))) {
    basemap_sf <- basemap
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Add basemap
  if (!is.null(basemap_sf)) {
    p <- p + ggplot2::geom_sf(
      data = basemap_sf,
      fill = basemap_fill,
      color = basemap_border,
      linewidth = basemap_lwd
    )
  }

  # Add hexagons
  if (is.null(fill)) {
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      fill = "#E69F00",
      color = grid_border,
      linewidth = grid_lwd,
      alpha = grid_alpha
    )
  } else {
    p <- p + ggplot2::geom_sf(
      data = hex_sf,
      ggplot2::aes(fill = .data[[fill]]),
      color = grid_border,
      linewidth = grid_lwd,
      alpha = grid_alpha
    )

    # Add default color scale
    if (is.numeric(hex_sf[[fill]])) {
      p <- p + ggplot2::scale_fill_viridis_c(name = legend_title %||% fill)
    } else {
      p <- p + ggplot2::scale_fill_viridis_d(name = legend_title %||% fill)
    }
  }

  # Add points if requested
  if (show_points) {
    if (inherits(x@data, "sf")) {
      p <- p + ggplot2::geom_sf(
        data = x@data,
        color = point_color,
        size = point_size
      )
    } else {
      # Use cell centers
      pts_df <- data.frame(
        lon = x@cell_center[, "lon"],
        lat = x@cell_center[, "lat"]
      )
      p <- p + ggplot2::geom_point(
        data = pts_df,
        ggplot2::aes(x = .data$lon, y = .data$lat),
        color = point_color,
        size = point_size
      )
    }
  }

  # Set coordinates
  p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)

  # Add title
  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title)
  }

  # Apply minimal theme
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  p
}

# Helper for NULL coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a
