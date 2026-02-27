# plot_methods.R
# Plot methods for HexData objects
#
# Provides default plotting with sensible styling options.
# Users can visualize hexified data without specifying grid parameters.

# =============================================================================
# POINT SIZE AND JITTERING HELPERS
# =============================================================================

#' Calculate point size based on hexagon size
#'
#' Calculates the cex value so that a single point covers a target
#' fraction of a hexagon cell. Points may overlap when multiple
#' points are in a cell.
#'
#' @param size Size specification: numeric (direct cex), or preset name
#' @param hex_sf sf object with cell polygons
#' @param xlim,ylim Plot limits (for coordinate-to-inches conversion)
#' @return Numeric cex value
#' @noRd
resolve_point_size <- function(size, hex_sf, xlim, ylim) {
 # Coverage targets: fraction of hex area covered by ONE point
 coverage_targets <- c(
   "tiny" = 0.02,
   "small" = 0.05,
   "normal" = 0.10,
   "large" = 0.20,
   "very large" = 0.35
 )

 # If numeric, return directly
 if (is.numeric(size)) {
   return(size)
 }

 # Resolve preset name
 size_lower <- tolower(size)
 if (size_lower == "verylarge") size_lower <- "very large"
 if (size_lower == "auto") size_lower <- "normal"

 if (!size_lower %in% names(coverage_targets)) {
   warning("Unknown point_size '", size, "', using 'normal'")
   size_lower <- "normal"
 }

 target_coverage <- coverage_targets[[size_lower]]

 # Get a representative hexagon (first one)
 hex_poly <- hex_sf[1, ]

 # Get hexagon bounding box in coordinate units
 hex_bbox <- sf::st_bbox(hex_poly)
 hex_width <- hex_bbox["xmax"] - hex_bbox["xmin"]
 hex_height <- hex_bbox["ymax"] - hex_bbox["ymin"]

 # Approximate hexagon area (regular hexagon ~0.866 of bounding box)
 hex_area_approx <- hex_width * hex_height * 0.866

 # Target area for one point
 point_area <- hex_area_approx * target_coverage

 # Point radius in coordinate units
 point_radius_coord <- sqrt(point_area / pi)

 # Convert to plot units
 plot_width <- xlim[2] - xlim[1]
 plot_height <- ylim[2] - ylim[1]

 # Get device dimensions in inches
 dev_size <- grDevices::dev.size("in")
 if (any(dev_size == 0)) dev_size <- c(7, 5)

 # cex=1 gives approximately 0.2 inches diameter
 base_point_diameter_in <- 0.2

 # Convert coordinate radius to inches
 x_scale <- dev_size[1] / plot_width
 y_scale <- dev_size[2] / plot_height
 avg_scale <- (x_scale + y_scale) / 2

 point_diameter_in <- 2 * point_radius_coord * avg_scale

 # Calculate cex
 cex <- point_diameter_in / base_point_diameter_in

 # Clamp to reasonable range
 max(0.1, min(5, cex))
}

#' Jitter points within their assigned hexagon cells
#' @param cell_ids Vector of cell IDs (one per point)
#' @param hex_sf sf object with cell polygons (must have cell_id column)
#' @return Data frame with lon, lat columns for jittered positions
#' @noRd
jitter_points_in_cells <- function(cell_ids, hex_sf) {
  n_points <- length(cell_ids)

  # Pre-allocate result
  result <- data.frame(
    lon = numeric(n_points),
    lat = numeric(n_points)
  )

  # Disable S2 for sampling (simpler planar sampling in small hexagons)
  s2_state <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(invisible(suppressMessages(sf::sf_use_s2(s2_state))), add = TRUE)

  # Group points by cell
  cell_groups <- split(seq_along(cell_ids), cell_ids)

  # Determine if cell_ids are character (H3) or numeric (ISEA)
  h3_mode <- is.character(hex_sf$cell_id)

  for (cell_id_str in names(cell_groups)) {
    cell_id_val <- if (h3_mode) cell_id_str else as.numeric(cell_id_str)
    indices <- cell_groups[[cell_id_str]]
    n_in_cell <- length(indices)

    # Get the polygon for this cell
    poly_idx <- which(hex_sf$cell_id == cell_id_val)
    if (length(poly_idx) == 0) next

    poly <- hex_sf[poly_idx, ]

    # Sample random points inside the polygon
    sampled <- suppressMessages(suppressWarnings(
      sf::st_sample(poly, size = n_in_cell, type = "random")
    ))

    # If sampling failed or returned wrong count, use centroid
    if (length(sampled) == 0) {
      centroid <- suppressMessages(suppressWarnings(
        sf::st_centroid(sf::st_geometry(poly))
      ))
      coords <- sf::st_coordinates(centroid)
      result$lon[indices] <- coords[1, 1]
      result$lat[indices] <- coords[1, 2]
    } else {
      coords <- sf::st_coordinates(sampled)
      # Handle case where we got fewer points than expected
      if (nrow(coords) < n_in_cell) {
        # Pad with centroid
        centroid <- suppressMessages(suppressWarnings(
          sf::st_centroid(sf::st_geometry(poly))
        ))
        centroid_coords <- sf::st_coordinates(centroid)
        n_missing <- n_in_cell - nrow(coords)
        coords <- rbind(coords, matrix(rep(centroid_coords, n_missing),
                                        ncol = 2, byrow = TRUE))
      }
      result$lon[indices] <- coords[1:n_in_cell, 1]
      result$lat[indices] <- coords[1:n_in_cell, 2]
    }
  }

  result
}

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
#' @param show_points Show original points on top of cells (default FALSE).
#'   Points are jittered within their assigned hexagon.
#' @param point_size Size of points. Can be:
#'   \itemize{
#'     \item A number (direct cex value)
#'     \item A preset defining what fraction of a hex cell one point covers:
#'       "tiny" (~2\%), "small" (~5\%), "normal"/"auto" (~10\%),
#'       "large" (~20\%), "very large" (~35\%)
#'   }
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
#' @seealso \code{\link{hexify_heatmap}} for ggplot2 plotting
#'
#' @export
#' @examples
#' \donttest{
#' df <- data.frame(lon = runif(50, -5, 5), lat = runif(50, 45, 50))
#' result <- hexify(df, lon = "lon", lat = "lat", area_km2 = 2000)
#'
#' # Basic plot
#' plot(result, basemap = FALSE)
#'
#' # With basemap and custom styling
#' plot(result, grid_fill = "lightblue", grid_border = "darkblue")
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
           point_size = "auto",
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
      if (clip_basemap) {
        # Disable S2 for clipping to avoid spherical geometry edge-crossing errors
        s2_state <- sf::sf_use_s2()
        suppressMessages(sf::sf_use_s2(FALSE))
        on.exit(sf::sf_use_s2(s2_state), add = TRUE)

        # Clip basemap to extent
        clip_box <- sf::st_bbox(c(xmin = unname(xlim[1]),
                                   xmax = unname(xlim[2]),
                                   ymin = unname(ylim[1]),
                                   ymax = unname(ylim[2])),
                                 crs = 4326)
        clip_poly <- sf::st_as_sfc(clip_box)

        basemap_clipped <- suppressMessages(suppressWarnings(
          sf::st_intersection(sf::st_make_valid(basemap_sf), clip_poly)
        ))

        if (nrow(basemap_clipped) > 0) {
          plot(sf::st_geometry(basemap_clipped),
               col = basemap_fill, border = basemap_border,
               lwd = basemap_lwd, add = TRUE)
        }
      } else {
        # No clipping - draw full basemap
        plot(sf::st_geometry(basemap_sf),
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
      # Resolve point size based on hex cell size
      cex <- resolve_point_size(point_size, hex_sf, xlim, ylim)

      # Jitter points within their hexagon
      jittered <- jitter_points_in_cells(x@cell_id, hex_sf)

      points(jittered$lon, jittered$lat,
             pch = 19, cex = cex, col = point_color)
    }

    invisible(x)
  }
)

# Helper for NULL coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# PLOT GRID CLIPPED TO BOUNDARY
# =============================================================================

#' Plot hexagonal grid clipped to a polygon boundary
#'
#' A convenience function that creates a grid, clips it to a boundary polygon,
#' and plots the result in a single call.
#'
#' @param boundary An sf/sfc polygon to clip to (e.g., country boundary)
#' @param grid A HexGridInfo object from \code{hex_grid()}
#' @param crop If TRUE (default), cells are cropped to boundary. If FALSE,
#'   only complete hexagons within boundary are shown.
#' @param boundary_fill Fill color for the boundary polygon (default "gray95")
#' @param boundary_border Border color for boundary (default "gray40")
#' @param boundary_lwd Line width for boundary (default 0.5)
#' @param grid_fill Fill color for grid cells (default "steelblue")
#' @param grid_border Border color for grid cells (default "steelblue")
#' @param grid_lwd Line width for cell borders (default 0.3)
#' @param grid_alpha Transparency for cell fill (0-1, default 0.3)
#' @param title Plot title. If NULL (default), auto-generates title with cell area.
#' @param expand Expansion factor for plot limits (default 0.05)
#'
#' @return A ggplot object that can be further customized
#'
#' @details
#' This is a convenience wrapper around \code{grid_clip()} that handles the
#' common use case of visualizing a hexagonal grid over a geographic region.
#'
#' @seealso \code{\link{grid_clip}} for the underlying clipping function,
#'   \code{\link{hex_grid}} for grid specification
#'
#' @export
#' @examples
#' \donttest{
#' # Plot grid over France
#' france <- hexify_world[hexify_world$name == "France", ]
#' grid <- hex_grid(area_km2 = 2000)
#' plot_grid(france, grid)
#'
#' # Customize colors
#' plot_grid(france, grid,
#'           grid_fill = "coral", grid_alpha = 0.5,
#'           boundary_fill = "lightyellow")
#'
#' # Keep only complete hexagons
#' plot_grid(france, grid, crop = FALSE)
#'
#' # Add ggplot2 customizations
#' library(ggplot2)
#' plot_grid(france, grid) +
#'   labs(subtitle = "ISEA3H Discrete Global Grid") +
#'   theme_void()
#' }
plot_grid <- function(boundary,
                      grid,
                      crop = TRUE,
                      boundary_fill = "gray95",
                      boundary_border = "gray40",
                      boundary_lwd = 0.5,
                      grid_fill = "steelblue",
                      grid_border = "steelblue",
                      grid_lwd = 0.3,
                      grid_alpha = 0.3,
                      title = NULL,
                      expand = 0.05) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_grid()")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for plot_grid()")
  }

  # Extract grid spec
  g <- extract_grid(grid)

  # Clip grid to boundary
  clipped_grid <- grid_clip(boundary, g, crop = crop)

  # Get bounding box for plot limits
  bbox <- sf::st_bbox(boundary)
  x_range <- bbox["xmax"] - bbox["xmin"]
  y_range <- bbox["ymax"] - bbox["ymin"]
  xlim <- c(bbox["xmin"] - x_range * expand, bbox["xmax"] + x_range * expand)
  ylim <- c(bbox["ymin"] - y_range * expand, bbox["ymax"] + y_range * expand)

  # Auto-generate title if not provided
  if (is.null(title)) {
    title <- sprintf("Hexagonal Grid (~%.0f km\u00b2 cells)", g@area_km2)
  }

  # Build plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = boundary,
      fill = boundary_fill,
      color = boundary_border,
      linewidth = boundary_lwd
    ) +
    ggplot2::geom_sf(
      data = clipped_grid,
      fill = grDevices::adjustcolor(grid_fill, alpha.f = grid_alpha),
      color = grid_border,
      linewidth = grid_lwd
    ) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  p
}
