# R/hex_browse.R
# Interactive leaflet-based hex maps

#' Interactive Hex Map
#'
#' Opens an interactive leaflet map with hexagonal cell polygons.
#' Cells can be colored by a data column, with popups showing cell
#' information on click.
#'
#' @param hex_data A HexData object, or a data.frame with a `cell_id` column.
#' @param grid A HexGridInfo object. Required if `hex_data` is a data.frame.
#' @param value Optional column name (character) to color cells by.
#'   If `NULL`, cells are colored uniformly.
#' @param palette Color palette name for continuous values. Default `"viridis"`.
#' @param opacity Fill opacity (0-1). Default 0.7.
#'
#' @return A `leaflet` map object (can be printed or embedded in Shiny).
#'
#' @details
#' Requires the `leaflet` package (in Suggests). The map is built using
#' [as_sf()] to generate polygon geometries, then rendered as a leaflet
#' choropleth.
#'
#' @seealso [hexify_heatmap()] for static ggplot2 maps,
#'   [as_sf()] for sf conversion
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("leaflet", quietly = TRUE)) {
#'   df <- data.frame(
#'     lon = runif(50, -5, 5),
#'     lat = runif(50, 45, 55),
#'     value = rnorm(50)
#'   )
#'   hd <- hexify(df, lon = "lon", lat = "lat", area_km2 = 500)
#'   hex_browse(hd, value = "value")
#' }
#' }
hex_browse <- function(hex_data, grid = NULL, value = NULL,
                        palette = "viridis", opacity = 0.7) {
  if (!requireNamespace("leaflet", quietly = TRUE)) {
    stop("Package 'leaflet' is required for hex_browse(). ",
         "Install with: install.packages('leaflet')")
  }

  # Get grid and cell data

  if (is_hex_data(hex_data)) {
    g <- hex_data@grid

    # Aggregate to unique cells
    df <- as.data.frame(hex_data)
    cell_ids <- unique(df$cell_id)

    # Build summary per cell
    cell_df <- data.frame(cell_id = cell_ids, stringsAsFactors = FALSE)
    first_idx <- match(cell_ids, df$cell_id)
    cell_df$cell_cen_lon <- df$cell_cen_lon[first_idx]
    cell_df$cell_cen_lat <- df$cell_cen_lat[first_idx]

    # Count points per cell
    cell_df$n_points <- as.integer(table(df$cell_id)[as.character(cell_ids)])

    # Add value column if specified
    if (!is.null(value)) {
      if (!value %in% names(df)) {
        stop(sprintf("Column '%s' not found in hex_data", value))
      }
      # Mean per cell
      agg <- tapply(df[[value]], df$cell_id, mean, na.rm = TRUE)
      cell_df[[value]] <- as.numeric(agg[as.character(cell_ids)])
    }
  } else if (is.data.frame(hex_data)) {
    if (is.null(grid)) stop("grid is required when hex_data is a data.frame")
    g <- extract_grid(grid)
    cell_df <- hex_data
    cell_ids <- cell_df$cell_id
  } else {
    stop("hex_data must be a HexData object or data.frame with cell_id column")
  }

  # Generate polygons
  hex_sf <- cell_to_sf(cell_ids, g)

  # Attach data to sf
  for (col in names(cell_df)) {
    if (col != "cell_id") {
      hex_sf[[col]] <- cell_df[[col]]
    }
  }

  # Build popup text
  popup_text <- paste0(
    "<b>Cell ID:</b> ", cell_ids, "<br>",
    "<b>Center:</b> (",
    round(cell_df$cell_cen_lon, 4), ", ",
    round(cell_df$cell_cen_lat, 4), ")<br>"
  )
  if ("n_points" %in% names(cell_df)) {
    popup_text <- paste0(popup_text, "<b>Points:</b> ", cell_df$n_points, "<br>")
  }
  if (!is.null(value) && value %in% names(cell_df)) {
    popup_text <- paste0(popup_text, "<b>", value, ":</b> ",
                          round(cell_df[[value]], 4), "<br>")
  }

  # Build leaflet map
  m <- leaflet::leaflet(hex_sf)
  m <- leaflet::addTiles(m)

  if (!is.null(value) && value %in% names(cell_df)) {
    vals <- cell_df[[value]]
    pal <- leaflet::colorNumeric(palette, domain = vals, na.color = "#808080")
    m <- leaflet::addPolygons(m,
                               fillColor = ~pal(hex_sf[[value]]),
                               fillOpacity = opacity,
                               weight = 1,
                               color = "#333333",
                               opacity = 0.5,
                               popup = popup_text)
    m <- leaflet::addLegend(m, pal = pal, values = vals, title = value)
  } else {
    m <- leaflet::addPolygons(m,
                               fillColor = "#3388ff",
                               fillOpacity = opacity,
                               weight = 1,
                               color = "#333333",
                               opacity = 0.5,
                               popup = popup_text)
  }

  m
}
