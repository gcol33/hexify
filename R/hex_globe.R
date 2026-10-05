# hex_globe.R
# An interactive globe of a grid, drawn with WebGPU

#' Interactive globe of a grid
#'
#' Draws the cells of a grid on a globe that turns under the mouse, rendered
#' on the graphics card through WebGPU. A slider folds the icosahedron the
#' grid is built on into the sphere, and the shader blends the two surfaces.
#' Values given per cell fill the cells through a colour ramp.
#'
#' The cells of an ISEA grid are found per pixel on the graphics card, by the
#' same projection and cell numbering as \code{\link{lonlat_to_cell}}, so the
#' page holds the grid's description and the values rather than the cells'
#' outlines, and a grid of any resolution draws as fast as a coarse one.
#' Borders keep their width at every zoom and fade out where the cells get
#' too small to see. The pointer shows the ID and value of the cell under it.
#' H3 cells are drawn from their outlines.
#'
#' Drag to turn the globe. Drag with Shift held to turn the view about the
#' line of sight and, in the perspective view, to tilt the camera. The mouse
#' wheel moves the perspective camera closer or zooms the orthographic view.
#' The view starts where the camera arguments, those of the grid's
#' \code{\link[=plot,HexGridInfo,missing-method]{plot}} method, place it.
#'
#' The globe needs the 'htmlwidgets' package and a viewer or browser with
#' WebGPU; without WebGPU the widget shows a notice instead.
#'
#' @inheritParams plot,HexGridInfo,missing-method
#' @param values Numeric values, one per cell in \code{cells}, filling the
#'   cells. \code{NULL} draws no fill.
#' @param cells Cell IDs to draw. \code{NULL} draws every cell of the grid.
#' @param surface Surface shown first, \code{"sphere"} or
#'   \code{"icosahedron"}; the slider folds one into the other. An H3 grid
#'   is drawn on the sphere only.
#' @param palette Colours of the ramp from low to high values, or the name of
#'   a palette of \code{\link[grDevices]{hcl.colors}}.
#' @param limits Values at the two ends of the ramp; \code{NULL} uses the
#'   range of \code{values}. Values outside are drawn in the end colours.
#' @param na_fill Fill colour of cells whose value is \code{NA}; \code{NA}
#'   leaves them unfilled.
#' @param width,height Size of the widget, as CSS units or pixels.
#' @param elementId Id of the widget's HTML element.
#'
#' @return An htmlwidget.
#'
#' @seealso \code{\link[=plot,HexGridInfo,missing-method]{plot}} for the same
#'   view as a static plot
#'
#' @export
#' @examples
#' if (requireNamespace("htmlwidgets", quietly = TRUE)) {
#'   grid <- hex_grid(resolution = 3, aperture = 3)
#'   hex_globe(grid)
#'   hex_globe(grid, surface = "icosahedron", center = "pacific")
#'
#'   cells <- seq_len(n_cells(grid))
#'   centres <- cell_to_lonlat(cells, grid)
#'   hex_globe(grid, values = centres$lat_deg, palette = "Blue-Red 3")
#' }
hex_globe <- function(x,
                      values = NULL,
                      cells = NULL,
                      surface = c("sphere", "icosahedron"),
                      center = c(lon = 15, lat = 32),
                      projection = c("orthographic", "perspective"),
                      distance = NULL,
                      tilt = 0,
                      rotation = 0,
                      fov = NULL,
                      land = TRUE,
                      ocean_fill = "#F4F6F8",
                      land_fill = "#C9D0D6",
                      land_border = "#8C969F",
                      land_lwd = 0.4,
                      grid_border = "#0072B2",
                      grid_lwd = 0.8,
                      face_edges = NULL,
                      edge_col = "black",
                      edge_lwd = 1.1,
                      palette = "viridis",
                      limits = NULL,
                      na_fill = NA,
                      step = 0.01,
                      width = NULL,
                      height = NULL,
                      elementId = NULL) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("hex_globe() needs the 'htmlwidgets' package. ",
         "Install with: install.packages('htmlwidgets')", call. = FALSE)
  }
  surface <- match.arg(surface)
  projection <- match.arg(projection)
  g <- extract_grid(x)
  face_edges <- resolve_surface(surface, face_edges, g)
  camera <- resolve_camera(projection, distance, tilt, rotation, fov)
  center <- resolve_center(center)
  land <- surface_land(land)
  arc <- step * atan(2)
  orient <- orient_arg(g)

  grid <- NULL
  fill <- NULL
  grid_lines <- NULL
  if (is_h3_grid(g)) {
    cells <- surface_cells(g, cells)
    if (!is.null(values)) {
      fill <- c(globe_mesh(h3_surface_mesh(cells, GLOBE_MESH_SPACING)),
                values = cpp_base64_buffer(ramp_position(values, cells, limits), "f32"))
    }
    grid_lines <- globe_lines(grid_surface_paths(g, cells, step))
  } else {
    grid <- globe_grid(g, cells, values, limits)
  }

  x <- list(
    shader = paste(readLines(system.file("wgsl", "globe.wgsl", package = "hexify")),
                   collapse = "\n"),
    camera = list(
      center = unname(center),
      projection = projection,
      distance = if (is.finite(camera$distance)) camera$distance,
      tilt = tilt,
      rotation = rotation,
      fov = if (!is.na(camera$fov)) camera$fov
    ),
    fold = if (surface == "sphere") 1 else 0,
    foldable = !is_h3_grid(g),
    surface = globe_mesh(cpp_globe_faces(orient, GLOBE_MESH_SPACING),
                         tri = !is.null(grid)),
    land = if (!is.null(land) && !is.na(land_fill)) {
      globe_mesh(cpp_globe_polygons(orient, sfc_polygons(land), GLOBE_MESH_SPACING))
    },
    grid = grid,
    projection = if (!is.null(grid)) cpp_globe_projection(orient),
    cells = fill,
    palette = as.vector(grDevices::col2rgb(ramp_colours(palette), alpha = TRUE)) / 255,
    grid_lines = grid_lines,
    cell_width = sqrt(4 * pi / grid_n_cells(g)),
    land_lines = if (!is.null(land) && !is.na(land_border)) {
      globe_lines(land_surface_paths(land, arc, orient))
    },
    edge_lines = if (face_edges) globe_lines(edge_surface_paths(arc, orient)),
    style = list(
      ocean_fill = globe_rgba(ocean_fill),
      land_fill = globe_rgba(land_fill),
      land_border = globe_rgba(land_border),
      land_lwd = land_lwd,
      grid_border = globe_rgba(grid_border),
      grid_lwd = grid_lwd,
      edge_col = globe_rgba(edge_col),
      edge_lwd = edge_lwd,
      na_fill = globe_rgba(na_fill)
    )
  )

  htmlwidgets::createWidget(
    "hex_globe", x, width = width, height = height, package = "hexify",
    elementId = elementId,
    sizingPolicy = htmlwidgets::sizingPolicy(
      defaultWidth = 640, defaultHeight = 640, padding = 0,
      browser.fill = TRUE
    )
  )
}

#' Save a globe as a PNG image
#'
#' Draws a globe made by \code{\link{hex_globe}} in headless Chrome, with the
#' same WebGPU renderer the widget uses, and saves the view it opens with,
#' without its controls, as a PNG file.
#'
#' Needs the 'chromote' and 'htmlwidgets' packages and a Chromium browser,
#' such as Chrome or Edge, that \code{chromote::find_chrome()} finds. The
#' browser is started with its GPU, which WebGPU draws on.
#'
#' @param widget A globe from \code{\link{hex_globe}}.
#' @param file Path of the PNG file to write.
#' @param width,height Size of the view in CSS pixels.
#' @param scale Image pixels per CSS pixel; 2 draws the view at twice the
#'   resolution.
#' @param timeout Seconds to wait for the globe to be drawn.
#'
#' @return \code{file}, invisibly.
#'
#' @seealso \code{\link[=plot,HexGridInfo,missing-method]{plot}} for a static
#'   plot that needs no browser
#'
#' @export
#' @examples
#' \dontrun{
#' grid <- hex_grid(resolution = 3, aperture = 3)
#' globe <- hex_globe(grid, surface = "icosahedron", center = "pacific")
#' hex_globe_png(globe, tempfile(fileext = ".png"), scale = 2)
#' }
hex_globe_png <- function(widget, file, width = 800, height = 800, scale = 1,
                          timeout = 30) {
  if (!inherits(widget, "hex_globe")) {
    stop("widget must be a globe from hex_globe()", call. = FALSE)
  }
  for (pkg in c("chromote", "htmlwidgets")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("hex_globe_png() needs the '", pkg, "' package. ",
           "Install with: install.packages('", pkg, "')", call. = FALSE)
    }
  }
  globe_in_chrome(widget, width, height, scale, timeout, function(session, read) {
    read("document.querySelectorAll('.hexify-globe-controls').forEach(e => e.remove()), ''")
    session$screenshot(filename = file, selector = ".hexify-globe", scale = scale,
                       show = FALSE)
  })
  invisible(file)
}

#' Open a globe in headless Chrome, with its GPU, wait until it is drawn and
#' call `f(session, read)`, where `read(js)` evaluates JavaScript in the page
#' (awaiting a promise) and returns its value
#' @noRd
globe_in_chrome <- function(widget, width, height, scale, timeout, f) {
  dir <- tempfile("hex_globe_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  html <- file.path(dir, "globe.html")
  widget$width <- width
  widget$height <- height
  htmlwidgets::saveWidget(widget, html, selfcontained = FALSE, libdir = "lib")

  args <- c(setdiff(chromote::default_chrome_args(), "--disable-gpu"),
            "--enable-unsafe-webgpu")
  browser <- chromote::Chromote$new(browser = chromote::Chrome$new(args = args))
  on.exit(browser$close(), add = TRUE)
  session <- browser$new_session(width = width, height = height)
  session$Emulation$setDeviceMetricsOverride(
    width = width, height = height, deviceScaleFactor = scale, mobile = FALSE
  )
  session$Page$navigate(paste0("file:///", normalizePath(html, winslash = "/")))

  read <- function(js) {
    r <- session$Runtime$evaluate(js, awaitPromise = TRUE, returnByValue = TRUE,
                                  timeout_ = timeout)
    if (!is.null(r$exceptionDetails)) {
      stop("in the page: ", r$exceptionDetails$text, " ",
           r$exceptionDetails$exception$description, call. = FALSE)
    }
    r$result$value
  }
  state <- ""
  deadline <- Sys.time() + timeout
  while (Sys.time() < deadline) {
    state <- read("(document.querySelector('.hexify-globe') || {dataset: {}}).dataset.state || ''")
    if (state %in% c("drawn", "failed")) break
    Sys.sleep(0.1)
  }
  if (state == "failed") {
    stop("the globe could not be drawn: ",
         read("document.querySelector('.hexify-globe').dataset.message"),
         call. = FALSE)
  }
  if (state != "drawn") {
    stop("the globe was not drawn within ", timeout, " seconds", call. = FALSE)
  }
  f(session, read)
}

#' Cell IDs of points found by the globe's shader, for testing it against
#' lonlat_to_cell(). The points go to the GPU as 32-bit floats; `xyz` are
#' their unit vectors as the GPU reads them.
#' @noRd
globe_shader_cells <- function(grid, xyz, timeout = 120) {
  widget <- hex_globe(grid, land = FALSE, face_edges = FALSE)
  chunks <- split(seq_len(nrow(xyz)), ceiling(seq_len(nrow(xyz)) / 2e5))
  out <- globe_in_chrome(widget, 64, 64, 1, timeout, function(session, read) {
    lapply(chunks, function(k) {
      b64 <- cpp_base64_buffer(as.vector(t(xyz[k, , drop = FALSE])), "f32")
      read(paste0("document.querySelector('.hexify-globe').hexGlobe.locate('", b64,
                  "').then(ids => Array.from(ids).join(','))"))
    })
  })
  as.numeric(unlist(strsplit(unlist(out), ",", fixed = TRUE)))
}

#' An ISEA grid as the widget reads it to find each pixel's cell
#'
#' The grid's frame (see cpp_globe_frame()) and the cells drawn, as their
#' sorted IDs split into high and low 32-bit words, each with its place on
#' the colour ramp when values are given. With no cells given every cell is
#' drawn, and IDs are sent only to carry values.
#' @noRd
globe_grid <- function(g, cells, values, limits) {
  mixed <- is_mixed_aperture(g@aperture)
  frame <- cpp_globe_frame(g@resolution,
                           if (mixed) 0L else aperture_to_int(g@aperture),
                           if (mixed) grid_ap_seq(g) else integer(0))
  # The shader reads a point to 32-bit float precision, a few millionths of
  # a radian; past this quad side its cells are finer than that.
  if (frame$dim >= 2^24) {
    stop("hex_globe() draws ISEA grids up to aperture 3 resolution 30, ",
         "aperture 4 resolution 23 and aperture 7 resolution 16; this grid's ",
         "cells are finer than the graphics card's 32-bit floats resolve",
         call. = FALSE)
  }
  all <- is.null(cells)
  # Values for the whole grid are read by cell ID; any other cells are found
  # among their sorted IDs.
  dense <- all && !is.null(values) && frame$n_cells <= 2^32
  if (all && !is.null(values) && !dense) cells <- seq_len(frame$n_cells)
  keys <- NULL
  if (dense) {
    check_values(values, frame$n_cells)
  } else if (!is.null(cells)) {
    if (!is.numeric(cells) || anyNA(cells) || any(cells != floor(cells)) ||
        any(cells < 1 | cells > frame$n_cells)) {
      stop("cells must be cell IDs of the grid, from 1 to ", frame$n_cells,
           call. = FALSE)
    }
    o <- order(cells)
    if (!is.null(values)) {
      check_values(values, length(cells))
      values <- values[o]
    }
    keys <- split_u64(as.numeric(cells)[o])
  }
  list(
    dim = frame$dim,
    index = frame$index,
    c = frame$c,
    generator = frame$generator,
    per_quad = split_u64(frame$per_quad),
    all = all,
    dense = dense,
    n_keys = length(keys) / 2,
    keys = if (!is.null(keys)) cpp_base64_buffer(keys, "u32"),
    values = if (!is.null(values)) cpp_base64_buffer(as.numeric(values), "f32"),
    ramp_map = if (!is.null(values)) ramp_map(values, limits)
  )
}

#' Whole numbers below 2^53 as their high and low 32-bit words, interleaved
#' @noRd
split_u64 <- function(x) {
  hi <- floor(x / 2^32)
  as.vector(rbind(hi, x - hi * 2^32))
}

#' H3 cells filled as one mesh of the sphere, item k being cell k of `cells`;
#' their edges are great-circle arcs, so each is a polygon of the sphere
#' @noRd
h3_surface_mesh <- function(cells, max_len) {
  rings <- lapply(cpp_h3_cellToBoundary(as.character(cells)), function(b) {
    b <- lonlat_ring_coords(b)
    if (any(b[1, ] != b[nrow(b), ])) b <- rbind(b, b[1, ])
    list(b)
  })
  cpp_globe_polygons(numeric(0), rings, max_len)
}

#' Coastlines placed on the faces, one path per ring
#' @noRd
land_surface_paths <- function(land, max_angle, orient) {
  rings <- unlist(sfc_polygons(land), recursive = FALSE)
  path <- rep(seq_along(rings), vapply(rings, nrow, integer(1)))
  ll <- do.call(rbind, rings)
  cpp_sphere_paths_on_faces(orient, ll[, 1], ll[, 2], path, max_angle)
}

#' The icosahedron's edges placed on the faces, one path per edge
#' @noRd
edge_surface_paths <- function(max_angle, orient) {
  solid <- icosa_solid(orient)
  V <- solid$vertices
  ends <- as.vector(t(solid$edges[, c("v1", "v2")]))
  lon <- atan2(V[ends, 2], V[ends, 1]) * 180 / pi
  lat <- asin(V[ends, 3]) * 180 / pi
  cpp_sphere_paths_on_faces(orient, lon, lat,
                            rep(seq_len(nrow(solid$edges)), each = 2L), max_angle)
}

#' Values must be numbers, one per cell
#' @noRd
check_values <- function(values, n) {
  if (!is.numeric(values) || length(values) != n) {
    stop("values must be numeric, one per cell (", n, ")", call. = FALSE)
  }
}

#' Where values sit along the colour ramp: a value v at
#' clamp((v - m[1]) * m[2] + m[3], 0, 1), the limits at its two ends, or
#' halfway when they coincide
#' @noRd
ramp_map <- function(values, limits) {
  if (is.null(limits)) {
    limits <- if (all(is.na(values))) c(0, 1) else range(values, na.rm = TRUE)
  }
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits)) {
    stop("limits must be two numbers", call. = FALSE)
  }
  span <- limits[2] - limits[1]
  if (span == 0) c(limits[1], 0, 0.5) else c(limits[1], 1 / span, 0)
}

#' Positions of values along the colour ramp: 0 to 1, -1 for NA
#' @noRd
ramp_position <- function(values, cells, limits) {
  check_values(values, length(cells))
  m <- ramp_map(values, limits)
  pos <- pmin(pmax((values - m[1]) * m[2] + m[3], 0), 1)
  pos[is.na(values)] <- -1
  pos
}

#' The 256 colours of a ramp
#' @noRd
ramp_colours <- function(palette) {
  key <- function(p) tolower(gsub("[^[:alnum:]]", "", p))
  if (length(palette) == 1L && key(palette) %in% key(grDevices::hcl.pals())) {
    return(grDevices::hcl.colors(256, palette))
  }
  if (length(palette) < 2L) {
    stop("palette must name a grDevices::hcl.colors() palette or give at ",
         "least two colours", call. = FALSE)
  }
  grDevices::colorRampPalette(palette, alpha = TRUE)(256)
}

#' A colour as red, green, blue and alpha from 0 to 1; NULL for NA
#' @noRd
globe_rgba <- function(col) {
  if (is.null(col) || is.na(col)) return(NULL)
  as.vector(grDevices::col2rgb(col, alpha = TRUE)) / 255
}

#' A mesh as the widget reads it: per vertex its icosahedron and sphere
#' positions (six 32-bit floats) and item (from 0), with `tri` also its
#' triangle coordinates on its face, and the triangles' vertex indices, each
#' as base64 text
#' @noRd
globe_mesh <- function(m, tri = FALSE) {
  pos <- rbind(matrix(m$solid, nrow = 3L), matrix(m$sphere, nrow = 3L))
  list(
    position = cpp_base64_buffer(as.vector(pos), "f32"),
    item = cpp_base64_buffer(m$item - 1L, "u32"),
    tri = if (tri) cpp_base64_buffer(m$tri, "f32"),
    index = cpp_base64_buffer(m$index, "u32"),
    n_index = length(m$index)
  )
}

#' Paths as the widget reads them: per point its icosahedron and sphere
#' positions, and the first point (from 0) of every segment between two
#' points of one path. A path with sphere positions only, as an H3 cell's,
#' takes them for both surfaces.
#' @noRd
globe_lines <- function(paths) {
  n <- nrow(paths)
  if (n < 2L) return(NULL)
  sphere <- paths[, c("sphere_x", "sphere_y", "sphere_z"), drop = FALSE]
  solid <- paths[, c("solid_x", "solid_y", "solid_z"), drop = FALSE]
  if (anyNA(solid)) solid <- sphere
  segment <- which(paths[-1L, "cell"] == paths[-n, "cell"]) - 1L
  list(
    position = cpp_base64_buffer(as.vector(t(cbind(solid, sphere))), "f32"),
    segment = cpp_base64_buffer(segment, "u32"),
    n_segment = length(segment)
  )
}
