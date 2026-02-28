# plot_globe.R
# Globe visualization for hexified data
#
# Renders hexagonized globe with orthographic projection

# =============================================================================
# CENTER PRESETS
# =============================================================================

#' Globe center presets
#'
#' Named list of lon/lat coordinates for common globe views.
#' Used by \code{\link{plot_globe}} when center is specified as a string.
#'
#' @format Named list with elements:
#' \describe{
#'   \item{europe}{c(10, 50) - Western/Central Europe}
#'   \item{north_america}{c(-100, 45) - USA and Canada}
#'   \item{south_america}{c(-60, -15) - Full continent}
#'   \item{africa}{c(20, 5) - Central Africa}
#'   \item{asia}{c(100, 35) - China, SE Asia, Japan}
#'   \item{oceania}{c(135, -25) - Australia, NZ, Indonesia}
#'   \item{middle_east}{c(45, 25) - Arabian Peninsula, Iran, Turkey}
#'   \item{south_asia}{c(80, 20) - India, Pakistan, Bangladesh}
#'   \item{pacific}{c(-160, -10) - Polynesia, Pacific islands}
#'   \item{caribbean}{c(-70, 18) - Caribbean islands}
#'   \item{arctic}{c(0, 90) - North pole view}
#'   \item{antarctic}{c(0, -90) - South pole view}
#' }
#'
#' @export
#' @examples
#' globe_centers$europe
#' globe_centers$oceania
globe_centers <- list(
  # Continental
  europe        = c(lon = 10,   lat = 50),
  north_america = c(lon = -100, lat = 45),
  south_america = c(lon = -60,  lat = -15),
  africa        = c(lon = 20,   lat = 5),
  asia          = c(lon = 100,  lat = 35),
  oceania       = c(lon = 135,  lat = -25),

  # Regional (fills gaps)
  middle_east   = c(lon = 45,   lat = 25),
  south_asia    = c(lon = 80,   lat = 20),
  pacific       = c(lon = -160, lat = -10),
  caribbean     = c(lon = -70,  lat = 18),

  # Polar
  arctic        = c(lon = 0,    lat = 90),
  antarctic     = c(lon = 0,    lat = -90)
)

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Plot hexagonized globe
#'
#' Renders a global hexagonal grid on an orthographic projection with
#' customizable rotation, land clipping, and styling options.
#'
#' @param area Cell area in km^2 (passed to \code{\link{hex_grid}})
#' @param center Globe center: either a preset name (e.g., "europe") or
#'   numeric vector c(lon, lat). See \code{\link{globe_centers}} for presets.
#' @param clip_to_land If TRUE, clip hexagons to land boundaries
#' @param land_data Optional sf object for land boundaries. If NULL and
#'   clip_to_land is TRUE, uses rnaturalearth::ne_countries()
#' @param exclude_antarctica If TRUE, exclude Antarctica from land clipping
#' @param fill Fill color for hexagons (default "#D4B896")
#' @param border Border color for hexagons (default "grey30")
#' @param border_width Border width for hexagons (default 0.2)
#' @param ocean_fill Fill color for ocean/globe background (default "white")
#' @param ocean_border Border color for globe circle (default "grey50")
#' @param show_land If TRUE, show land boundaries (default TRUE when clipping)
#' @param land_fill Fill color for land (default NA, transparent)
#' @param land_border Border color for land boundaries (default "grey40")
#' @param land_width Border width for land boundaries (default 0.3)
#' @param use_ggplot NULL = auto-detect, TRUE = force ggplot2, FALSE = force base
#' @param return_data If TRUE, return sf objects instead of plotting
#' @param aperture Grid aperture (default 3L)
#'
#' @return If use_ggplot = TRUE: ggplot2 object (can add layers with +)
#'   If use_ggplot = FALSE: NULL invisibly (plots directly)
#'   If return_data = TRUE: list of sf objects (hexagons, land, ocean_circle, crs)
#'
#' @details
#' The function handles several technical challenges:
#' \itemize{
#'   \item Hexagons on the back side of the globe fail to transform - these are
#'     filtered out gracefully
#'   \item Invalid geometries after projection are repaired with st_buffer(0)
#'   \item Clipping is done in orthographic CRS to avoid topology errors
#' }
#'
#' @seealso \code{\link{globe_centers}} for available presets,
#'   \code{\link{grid_global}} for generating global grids without plotting
#'
#' @export
#' @examples
#' \donttest{
#' # Get data for custom plotting (no rendering)
#' data <- plot_globe(area = 100000, center = "europe", return_data = TRUE)
#' nrow(data$hexagons)
#' class(data$ocean_circle)
#'
#' # Basic usage - Europe-centered globe
#' plot_globe(area = 80000, center = "europe")
#' }
plot_globe <- function(
    area = 50000,
    center = "europe",
    clip_to_land = FALSE,
    land_data = NULL,
    exclude_antarctica = TRUE,
    fill = "#D4B896",
    border = "grey30",
    border_width = 0.2,
    ocean_fill = "white",
    ocean_border = "grey50",
    show_land = clip_to_land,
    land_fill = NA,
    land_border = "grey40",
    land_width = 0.3,
    use_ggplot = NULL,
    return_data = FALSE,
    aperture = 3L
) {
  # Check sf dependency

if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  # Resolve center
  center_coords <- resolve_center(center)

  # Auto-detect ggplot2 if not specified
  if (is.null(use_ggplot)) {
    use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  }

  # Prepare data
  data <- prepare_globe_data(
    area = area,
    center = center_coords,
    clip_to_land = clip_to_land,
    land_data = land_data,
    exclude_antarctica = exclude_antarctica,
    show_land = show_land,
    aperture = aperture
  )

  # Return data if requested
  if (return_data) {
    return(data)
  }

  # Plot with appropriate backend
  if (use_ggplot) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required when use_ggplot = TRUE")
    }
    plot_globe_ggplot(
      data = data,
      fill = fill,
      border = border,
      border_width = border_width,
      ocean_fill = ocean_fill,
      ocean_border = ocean_border,
      show_land = show_land,
      land_fill = land_fill,
      land_border = land_border,
      land_width = land_width
    )
  } else {
    plot_globe_base(
      data = data,
      fill = fill,
      border = border,
      ocean_fill = ocean_fill,
      ocean_border = ocean_border,
      show_land = show_land,
      land_fill = land_fill,
      land_border = land_border
    )
  }
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Resolve center preset or coordinates
#' @noRd
resolve_center <- function(center) {
  if (is.character(center)) {
    if (!center %in% names(globe_centers)) {
      valid <- paste(names(globe_centers), collapse = ", ")
      stop(sprintf("Unknown center preset '%s'. Valid options: %s", center, valid))
    }
    return(globe_centers[[center]])
  }

  if (is.numeric(center) && length(center) == 2) {
    # Handle named or unnamed numeric vector
    if (is.null(names(center))) {
      return(c(lon = center[1], lat = center[2]))
    }
    return(center)
  }

  stop("center must be a preset name or numeric c(lon, lat)")
}

#' Prepare globe data for plotting
#' @noRd
prepare_globe_data <- function(
    area,
    center,
    clip_to_land,
    land_data,
    exclude_antarctica,
    show_land,
    aperture
) {
  # Create orthographic CRS
  crs_string <- sprintf(
    "+proj=ortho +lat_0=%f +lon_0=%f +x_0=0 +y_0=0 +datum=WGS84",
    center["lat"], center["lon"]
  )

  # Disable S2 early - some polar geometries are invalid in spherical mode
  s2_state <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_state), add = TRUE)

  # Generate global grid — skip dateline wrapping because orthographic
  # projection handles antimeridian cells natively; wrapping splits them
  # into MULTIPOLYGON halves that create visible gaps on the globe
  grid <- hex_grid(area_km2 = area, aperture = aperture)
  global_hex <- grid_global(grid, wrap_dateline = FALSE)

  # Pre-filter hexagons by angular distance (optimization)
  # Keep cells within ~95 degrees of center (visible hemisphere + margin)
  centers_ll <- suppressWarnings(sf::st_centroid(global_hex))
  coords <- sf::st_coordinates(centers_ll)

  # Angular distance calculation (spherical)
  lon1 <- center["lon"] * pi / 180
  lat1 <- center["lat"] * pi / 180
  lon2 <- coords[, 1] * pi / 180
  lat2 <- coords[, 2] * pi / 180

  # Haversine-based angular distance
  delta_lon <- lon2 - lon1
  delta_lat <- lat2 - lat1
  a <- sin(delta_lat / 2)^2 + cos(lat1) * cos(lat2) * sin(delta_lon / 2)^2
  angular_dist <- 2 * asin(sqrt(pmin(1, a))) * 180 / pi

  # Keep hexagons within visible hemisphere (with margin for edge cells)
  visible_idx <- angular_dist < 95
  global_hex <- global_hex[visible_idx, ]

  # Transform hexagons to orthographic - handle failures gracefully
  hex_ortho <- tryCatch({
    suppressWarnings(sf::st_transform(global_hex, crs_string))
  }, error = function(e) {
    # If batch transform fails, do one by one
    transform_individually(global_hex, crs_string)
  })

  # Remove empty geometries
  hex_ortho <- hex_ortho[!sf::st_is_empty(hex_ortho), ]

  # Fix invalid geometries using st_make_valid
  # MULTIPOLYGON cells (from dateline crossing) often become invalid after
  # orthographic transform due to self-intersecting parts - st_make_valid fixes this
  # Some cells may become degenerate (2-point rings) and fail st_make_valid
  hex_ortho <- tryCatch({
    sf::st_make_valid(hex_ortho)
  }, error = function(e) {
    # Batch st_make_valid failed - process cells individually
    # Mark cells for removal if they can't be fixed
    keep_mask <- rep(TRUE, nrow(hex_ortho))

    for (i in seq_len(nrow(hex_ortho))) {
      fixed <- tryCatch({
        sf::st_make_valid(hex_ortho[i, ])
      }, error = function(e2) {
        # Try buffer as fallback
        tryCatch({
          sf::st_buffer(hex_ortho[i, ], 0)
        }, error = function(e3) {
          NULL  # Mark for removal
        })
      })

      if (is.null(fixed) || sf::st_is_empty(fixed)) {
        keep_mask[i] <- FALSE
      } else {
        hex_ortho[i, ] <- fixed
      }
    }

    hex_ortho[keep_mask, ]
  })

  # Remove empty and degenerate geometries
  # Cells at hemisphere edge get "squashed" to tiny areas - filter them out
  hex_ortho <- hex_ortho[!sf::st_is_empty(hex_ortho), ]
  areas <- tryCatch(as.numeric(sf::st_area(hex_ortho)), error = function(e) rep(NA, nrow(hex_ortho)))
  # Expected cell area is roughly area * 1e6 m² - filter cells < 1% of expected
  min_area <- area * 1e6 * 0.01  # 1% of expected area
  hex_ortho <- hex_ortho[!is.na(areas) & areas > min_area, ]

  # Create ocean circle (globe boundary)
  earth_radius <- 6371000  # meters
  ocean_circle <- create_globe_circle(earth_radius, crs_string)

  # Clip hexagons to globe circle - use st_crop for robustness
  hex_ortho <- tryCatch({
    # Get bounding box of ocean circle for cropping
    bbox <- sf::st_bbox(ocean_circle)
    cropped <- suppressWarnings(sf::st_crop(hex_ortho, bbox))
    # Then intersect with circle
    suppressWarnings(sf::st_intersection(cropped, ocean_circle))
  }, error = function(e) {
    # Fallback: just filter by centroid within circle
    tryCatch({
      centroids <- suppressWarnings(sf::st_centroid(hex_ortho))
      within <- suppressWarnings(sf::st_within(centroids, ocean_circle, sparse = FALSE))
      hex_ortho[apply(within, 1, any), ]
    }, error = function(e2) hex_ortho)
  })

  # Keep only polygon geometries
  geom_types <- sf::st_geometry_type(hex_ortho)
  hex_ortho <- hex_ortho[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]

  # Remove any empty geometries that may have resulted from clipping
  hex_ortho <- hex_ortho[!sf::st_is_empty(hex_ortho), ]

  # Handle land data
  land_ortho <- NULL
  if (clip_to_land || show_land) {
    land_ortho <- prepare_land_data(
      land_data = land_data,
      exclude_antarctica = exclude_antarctica,
      crs_string = crs_string,
      ocean_circle = ocean_circle
    )

    # Clip hexagons to land if requested
    if (clip_to_land && !is.null(land_ortho)) {
      # Filter out any malformed geometries (NA validity) before intersection
      validity <- sf::st_is_valid(hex_ortho)
      hex_ortho <- hex_ortho[!is.na(validity) & validity, ]

      hex_ortho <- tryCatch({
        result <- suppressWarnings(sf::st_intersection(hex_ortho, land_ortho))
        geom_types <- sf::st_geometry_type(result)
        result[geom_types %in% c("POLYGON", "MULTIPOLYGON"), ]
      }, error = function(e) hex_ortho)
    }
  }

  list(
    hexagons = hex_ortho,
    land = land_ortho,
    ocean_circle = ocean_circle,
    crs = crs_string,
    center = center
  )
}

#' Transform geometries individually (fallback for batch failures)
#' @noRd
transform_individually <- function(sf_obj, crs) {
  valid_idx <- logical(nrow(sf_obj))
  geoms <- vector("list", nrow(sf_obj))

  for (i in seq_len(nrow(sf_obj))) {
    result <- tryCatch({
      suppressWarnings(sf::st_transform(sf_obj[i, ], crs))
    }, error = function(e) NULL)

    if (!is.null(result) && !sf::st_is_empty(result)) {
      valid_idx[i] <- TRUE
      geoms[[i]] <- sf::st_geometry(result)[[1]]
    }
  }

  if (sum(valid_idx) == 0) {
    stop("No hexagons could be transformed to orthographic projection")
  }

  sf_obj <- sf_obj[valid_idx, ]
  sf::st_geometry(sf_obj) <- sf::st_sfc(geoms[valid_idx], crs = crs)
  sf_obj
}

#' Create globe boundary circle
#' @noRd
create_globe_circle <- function(radius, crs) {
  # Create circle in orthographic CRS
  angles <- seq(0, 2 * pi, length.out = 100)
  x <- radius * cos(angles)
  y <- radius * sin(angles)
  coords <- cbind(x, y)
  coords <- rbind(coords, coords[1, ])  # Close polygon

  circle <- sf::st_polygon(list(coords))
  sf::st_sfc(circle, crs = crs)
}

#' Prepare land data for globe
#' @noRd
prepare_land_data <- function(land_data, exclude_antarctica, crs_string, ocean_circle) {
  # Disable S2 for land data processing
  s2_state <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_state), add = TRUE)

  # Get land data
  if (is.null(land_data)) {
    if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
      warning("Package 'rnaturalearth' required for land boundaries. ",
              "Install with: install.packages('rnaturalearth')")
      return(NULL)
    }
    land_data <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  }

  # Exclude Antarctica if requested
  if (exclude_antarctica) {
    # Try common column names for country filtering
    if ("sovereignt" %in% names(land_data)) {
      land_data <- land_data[land_data$sovereignt != "Antarctica", ]
    } else if ("name" %in% names(land_data)) {
      land_data <- land_data[land_data$name != "Antarctica", ]
    } else if ("NAME" %in% names(land_data)) {
      land_data <- land_data[land_data$NAME != "Antarctica", ]
    }
  }

  # Transform individual countries to orthographic FIRST, then union
  # (Union in lat/lon then transform creates empty geometries for global polygons)
  land_ortho <- tryCatch({
    transformed <- suppressWarnings(sf::st_transform(land_data, crs_string))
    transformed <- transformed[!sf::st_is_empty(transformed), ]

    if (nrow(transformed) == 0) return(NULL)

    # Union in orthographic space - use st_combine + st_union for robustness
    land_geom <- tryCatch({
      combined <- sf::st_combine(sf::st_geometry(transformed))
      suppressWarnings(sf::st_union(combined))
    }, error = function(e) {
      # If union fails, just return combined geometries
      sf::st_combine(sf::st_geometry(transformed))
    })

    suppressWarnings(sf::st_buffer(land_geom, 0))
  }, error = function(e) NULL)

  if (is.null(land_ortho) || length(land_ortho) == 0) return(NULL)

  # Remove empty geometries
  land_ortho <- land_ortho[!sf::st_is_empty(land_ortho)]
  if (length(land_ortho) == 0) return(NULL)

  # Clip to globe circle
  land_ortho <- tryCatch({
    result <- suppressWarnings(sf::st_intersection(land_ortho, ocean_circle))
    sf::st_buffer(result, 0)
  }, error = function(e) land_ortho)

  # Fix any degenerate geometries (e.g., rings with < 3 points)
  # Use a tiny buffer (1 meter) to clean up degenerate rings, then unbuffer
  land_ortho <- tryCatch({
    buffered <- sf::st_buffer(land_ortho, 1)  # 1 meter buffer
    unbuffered <- sf::st_buffer(buffered, -1)  # Remove the buffer
    sf::st_make_valid(unbuffered)
  }, error = function(e) {
    # Fallback: just try st_make_valid
    tryCatch(sf::st_make_valid(land_ortho), error = function(e2) land_ortho)
  })

  land_ortho
}

# =============================================================================
# PLOTTING BACKENDS
# =============================================================================

#' Plot globe with ggplot2
#' @noRd
plot_globe_ggplot <- function(
    data,
    fill,
    border,
    border_width,
    ocean_fill,
    ocean_border,
    show_land,
    land_fill,
    land_border,
    land_width
) {
  ggplot2::ggplot() +
    # Ocean/globe background
    ggplot2::geom_sf(
      data = sf::st_sf(geometry = data$ocean_circle),
      fill = ocean_fill,
      color = ocean_border,
      linewidth = 0.5
    ) +
    # Land boundaries (if showing)
    {if (show_land && !is.null(data$land))
      ggplot2::geom_sf(
        data = sf::st_sf(geometry = data$land),
        fill = land_fill,
        color = land_border,
        linewidth = land_width
      )
    } +
    # Hexagons
    ggplot2::geom_sf(
      data = data$hexagons,
      fill = fill,
      color = border,
      linewidth = border_width
    ) +
    ggplot2::coord_sf(crs = data$crs) +
    ggplot2::theme_void() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
}

#' Plot globe with base R graphics
#' @noRd
plot_globe_base <- function(
    data,
    fill,
    border,
    ocean_fill,
    ocean_border,
    show_land,
    land_fill,
    land_border
) {
  # Set up plot area
  bbox <- sf::st_bbox(data$ocean_circle)
  plot(NULL,
       xlim = c(bbox["xmin"], bbox["xmax"]),
       ylim = c(bbox["ymin"], bbox["ymax"]),
       asp = 1,
       axes = FALSE,
       xlab = "",
       ylab = "")

  # Draw ocean circle
  plot(data$ocean_circle, col = ocean_fill, border = ocean_border, add = TRUE)

  # Draw land if showing
  if (show_land && !is.null(data$land)) {
    plot(data$land, col = land_fill, border = land_border, add = TRUE)
  }

  # Draw hexagons
  plot(sf::st_geometry(data$hexagons), col = fill, border = border, add = TRUE)

  invisible(NULL)
}
