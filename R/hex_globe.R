# hex_globe.R
# An interactive globe of a grid, drawn with WebGPU

#' Interactive globe of a grid
#'
#' Draws the cells of a grid on a globe that turns under the mouse, rendered
#' on the graphics card through WebGPU. A slider folds the solid the
#' grid is built on into the sphere, and the shader blends the two surfaces.
#' Values given per cell fill the cells through a colour ramp.
#'
#' The cells of an ISEA grid are found per pixel on the graphics card, by the
#' same projection and cell numbering as \code{\link{lonlat_to_cell}}, so the
#' page holds the grid's description and the values rather than the cells'
#' outlines, and a grid of any resolution draws as fast as a coarse one.
#' The values sit on the graphics card as a texture laid out like the grid's
#' quads, one block per quad and resolution, which the shader reads at the
#' place it finds a cell. Where the cells shrink below a pixel the fill is
#' read from the grid one resolution coarser, and so on, whose cells take the
#' mean of the values in them, so a fine grid seen whole shows its pattern
#' rather than noise. Borders keep their width at every zoom and fade out
#' where the cells get too small to see. The pointer shows the ID and value
#' of the cell under it. H3 cells are drawn from their outlines.
#'
#' Drag to turn the globe. Drag with Shift held to turn the view about the
#' line of sight and, in the perspective view, to tilt the camera. The mouse
#' wheel moves the perspective camera closer or zooms the orthographic view.
#' The view starts where the camera arguments, those of the grid's
#' \code{\link[=plot,HexGridInfo,missing-method]{plot}} method, place it.
#'
#' The globe needs the 'htmlwidgets' package and a viewer or browser with
#' WebGPU: Chrome and Edge (on Linux with some graphics cards only), Firefox
#' on Windows and macOS, Safari 26, and RStudio's viewer. Without WebGPU the
#' widget shows a notice instead. The article
#' \url{https://gillescolling.com/hexify/articles/globe.html} lists the
#' versions.
#'
#' @inheritParams plot,HexGridInfo,missing-method
#' @param values Numeric values, one per cell in \code{cells}, filling the
#'   cells. \code{NULL} draws no fill.
#' @param cells Cell IDs to draw. \code{NULL} draws every cell of the grid.
#' @param surface Surface shown first, \code{"sphere"} or \code{"solid"},
#'   the flat faces of the grid's polyhedron; the slider folds one into the
#'   other. An H3 grid
#'   is drawn on the sphere only.
#' @param palette Colours of the ramp from low to high values, or the name of
#'   a palette of \code{\link[grDevices]{hcl.colors}}.
#' @param limits Values at the two ends of the ramp; \code{NULL} uses the
#'   range of \code{values}. Values outside are drawn in the end colours.
#' @param na_fill Fill colour of cells whose value is \code{NA}; \code{NA}
#'   leaves them unfilled.
#' @param smooth Whether the fill of an ISEA grid blends between the values
#'   of neighbouring cells, linearly across the triangles of cell centres,
#'   rather than filling each cell with its own value.
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
#'   hex_globe(grid, surface = "solid", center = "pacific")
#'
#'   cells <- seq_len(n_cells(grid))
#'   centres <- cell_to_lonlat(cells, grid)
#'   hex_globe(grid, values = centres$lat_deg, palette = "Blue-Red 3")
#' }
hex_globe <- function(x,
                      values = NULL,
                      cells = NULL,
                      surface = c("sphere", "solid"),
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
                      smooth = FALSE,
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
  icosa <- icosa_arg(g)

  grid <- NULL
  fill <- NULL
  grid_lines <- NULL
  if (!isTRUE(smooth) && !isFALSE(smooth)) {
    stop("smooth must be TRUE or FALSE", call. = FALSE)
  }
  if (is_h3_grid(g)) {
    if (smooth) stop("smooth blends the cells of ISEA grids only", call. = FALSE)
    cells <- grid_cells(g, cells)
    if (!is.null(values)) {
      fill <- c(globe_mesh(h3_surface_mesh(cells, GLOBE_MESH_SPACING)),
                values = cpp_base64_buffer(ramp_position(values, cells, limits), "f32"))
    }
    grid_lines <- globe_lines(grid_surface_paths(g, cells, step))
  } else {
    grid <- globe_grid(g, cells, values, limits, smooth, globe_rgba(na_fill))
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
    lift = GLOBE_LIFT,
    foldable = !is_h3_grid(g),
    surface = globe_mesh(cpp_globe_faces(icosa, GLOBE_MESH_SPACING),
                         tri = !is.null(grid)),
    land = if (!is.null(land) && !is.na(land_fill)) {
      globe_mesh(cpp_globe_polygons(icosa, sfc_polygons(land), GLOBE_MESH_SPACING))
    },
    grid = grid,
    cells = fill,
    palette = as.vector(grDevices::col2rgb(ramp_colours(palette), alpha = TRUE)) / 255,
    grid_lines = grid_lines,
    cell_width = sqrt(4 * pi / grid_n_cells(g)),
    land_lines = if (!is.null(land) && !is.na(land_border)) {
      globe_lines(land_surface_paths(land, arc, icosa))
    },
    edge_lines = if (face_edges) globe_lines(edge_surface_paths(arc, icosa)),
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
#' Draws a globe made by \code{\link{hex_globe}} with the same WebGPU
#' shader and pipelines the widget uses, and saves the view it opens with,
#' without its controls, as a PNG file. The image is drawn on the graphics
#' card and read back from it, so it is the same on a machine without a
#' display. Pixels off the globe are transparent.
#'
#' \code{renderer = "wgpu"} draws through 'wgpu', the 'Rust' implementation
#' of WebGPU, in the 'hexglobe' package, with no browser.
#' \code{renderer = "chrome"} draws in headless Chrome, which needs the
#' 'chromote' and 'htmlwidgets' packages and a Chromium browser, such as
#' Chrome or Edge, that \code{chromote::find_chrome()} finds; the browser is
#' started with its GPU, which WebGPU draws on.
#'
#' @param widget A globe from \code{\link{hex_globe}}.
#' @param file Path of the PNG file to write.
#' @param width,height Size of the view in CSS pixels.
#' @param scale Image pixels per CSS pixel; 2 draws the view at twice the
#'   resolution.
#' @param timeout Seconds to wait for the globe to be drawn in Chrome.
#' @param renderer \code{"wgpu"} or \code{"chrome"}, as above.
#'
#' @return \code{file}, invisibly.
#'
#' @seealso \code{\link[=plot,HexGridInfo,missing-method]{plot}} for a static
#'   plot that needs no graphics card
#'
#' @export
#' @examples
#' \dontrun{
#' grid <- hex_grid(resolution = 3, aperture = 3)
#' globe <- hex_globe(grid, surface = "solid", center = "pacific")
#' hex_globe_png(globe, tempfile(fileext = ".png"), scale = 2)
#' }
hex_globe_png <- function(widget, file, width = 800, height = 800, scale = 1,
                          timeout = 30, renderer = c("wgpu", "chrome")) {
  if (!inherits(widget, "hex_globe")) {
    stop("widget must be a globe from hex_globe()", call. = FALSE)
  }
  renderer <- match.arg(renderer)
  pkgs <- if (renderer == "wgpu") "hexglobe" else c("chromote", "htmlwidgets")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install <- if (pkg == "hexglobe") {
        "remotes::install_github('gcol33/hexglobe')"
      } else {
        paste0("install.packages('", pkg, "')")
      }
      stop("hex_globe_png(renderer = \"", renderer, "\") needs the '", pkg,
           "' package. Install with: ", install, call. = FALSE)
    }
  }
  if (renderer == "wgpu") {
    image <- hexglobe::render_scene(globe_scene(widget$x, width, height, scale))
    hexglobe::write_png(image, file)
    return(invisible(file))
  }
  png <- globe_in_chrome(widget, width, height, scale, timeout, function(session, read) {
    read("document.querySelector('.hexify-globe').hexGlobe.snapshot()")
  })
  writeBin(cpp_base64_decode(png), file)
  invisible(file)
}

#' A globe as hexglobe::render_scene() draws it
#'
#' The widget's layers, in the widget's drawing order (`buildLayers()` in
#' inst/htmlwidgets/hex_globe.js), each with its uniform and buffers as raw
#' vectors, and the camera of the view the widget opens with, for an image
#' of `width` by `height` CSS pixels at `scale` pixels each. A mesh shared by
#' two layers is one raw vector, which the renderer uploads once.
#' @noRd
globe_scene <- function(x, width, height, scale) {
  w <- as.integer(round(width * scale))
  h <- as.integer(round(height * scale))
  s <- x$style
  b64 <- function(text) if (!is.null(text)) cpp_base64_decode(text)
  meshes <- list()
  mesh <- function(name) {
    if (is.null(meshes[[name]])) {
      m <- x[[name]]
      meshes[[name]] <<- list(pos = b64(m$position), item = b64(m$item),
                              tri = b64(m$tri), index = b64(m$index),
                              count = m$n_index)
    }
    meshes[[name]]
  }
  palette <- le32(x$palette, "f32")
  none <- le32(numeric(4), "f32")
  clear <- numeric(4)

  uniform <- function(color, lift, width, shaded, ramped, fade = 0) {
    le32(c(color, lift, width * scale, shaded, ramped, fade * scale, 0, 0, 0), "f32")
  }
  mesh_layer <- function(name, lift, color, kind, values = NULL) {
    if (is.null(x[[name]]) || is.null(color)) return(NULL)
    c(mesh(name), list(kind = kind,
                       uniform = uniform(color, lift, 0, TRUE, !is.null(values)),
                       values = if (is.null(values)) none else values,
                       ramp = if (is.null(values)) none else palette))
  }
  line_layer <- function(lines, lift, color, width, fade = 0) {
    if (is.null(lines) || is.null(color) || lines$n_segment == 0 || !isTRUE(width > 0)) {
      return(NULL)
    }
    list(kind = 3L, count = lines$n_segment,
         uniform = uniform(color, lift, width, FALSE, FALSE, fade),
         points = b64(lines$position), segment = b64(lines$segment))
  }
  grid_layer <- function(g) {
    color <- s$grid_border
    c(mesh("surface"), list(
      kind = 2L,
      uniform = uniform(color %||% clear, x$lift$cells,
                        if (is.null(color)) 0 else s$grid_lwd, TRUE, g$valued),
      ramp = palette,
      grid = b64(g$uniform),
      textures = lapply(g$textures, function(t) {
        list(binding = t$binding, format = t$format, dimension = t$dimension,
             size = t$size, data = b64(t$data))
      })
    ))
  }

  layers <- list(
    mesh_layer("surface", x$lift$ocean, s$ocean_fill, 0L),
    mesh_layer("land", x$lift$land, s$land_fill, 0L),
    if (!is.null(x$grid)) grid_layer(x$grid),
    if (!is.null(x$cells)) {
      mesh_layer("cells", x$lift$cells, s$na_fill %||% clear, 1L, b64(x$cells$values))
    },
    line_layer(x$grid_lines, x$lift$grid, s$grid_border, s$grid_lwd, x$cell_width),
    line_layer(x$land_lines, x$lift$coast, s$land_border, s$land_lwd),
    line_layer(x$edge_lines, x$lift$edges, s$edge_col, s$edge_lwd)
  )
  list(shader = x$shader, camera = globe_camera_uniform(x, w, h), width = w,
       height = h, layers = Filter(Negate(is.null), layers))
}

#' The Camera uniform of globe.wgsl for the view a globe opens with
#'
#' The view of surface_view() and its frame from view_frame(), as the
#' widget's `writeCamera()` lays them out, for an image of `w` by `h` pixels.
#' @noRd
globe_camera_uniform <- function(x, w, h) {
  cam <- x$camera
  persp <- cam$projection == "perspective"
  view <- surface_view(c(lon = cam$center[1], lat = cam$center[2]),
                       distance = if (persp) cam$distance %||% 3 else Inf,
                       tilt = if (persp) cam$tilt else 0,
                       rotation = cam$rotation)
  frame <- view_frame(view, if (persp && !is.null(cam$fov)) cam$fov else NA)
  eye <- view$eye
  le32(c(view$cam[, 1], 0, view$cam[, 2], 0, view$cam[, 3], 0,
         if (is.null(eye)) c(0, 0, 0, 0) else c(eye, 1),
         frame, view$scale,
         w, h,
         if (is.null(eye)) c(0, 1) else c(view$near, sqrt(sum(eye^2)) + 1.5),
         view$light, x$fold), "f32")
}

#' The Grid uniform of globe.wgsl: 1296 32-bit words, the face and edge
#' tables padded to the largest solid's and the levels to 32
#'
#' `levels` has a row per level of the table, finest first, with the columns
#' of cpp_globe_table()'s `levels`; `table` is that table, NULL for none,
#' when level 0 is the frame alone.
#' @noRd
globe_grid_uniform <- function(projection, frame, levels, table, all, valued, smooth,
                               per_quad, na_fill, ramp_map) {
  out <- raw(5184)
  put <- function(word, bytes) out[4 * word + seq_along(bytes)] <<- bytes
  put(0, le32(projection$constants[1:12], "f32"))
  put(12, le32(c(all, valued, smooth, projection$n_faces), "u32"))
  put(16, le32(c(per_quad, if (is.null(table)) 0 else nrow(levels),
                 isTRUE(table$keyed)), "u32"))
  if (!is.null(table)) put(20, le32(c(table$m, table$seed, table$buckets, 0), "u32"))
  put(24, le32(na_fill %||% numeric(4), "f32"))
  put(28, le32(ramp_map %||% numeric(4), "f32"))
  put(32, le32(projection$faces, "f32"))
  put(352, le32(projection$edges, "i32"))
  if (is.null(levels)) {
    levels <- cbind(dim = frame$dim, index = frame$index, c = frame$c,
                    ga = frame$generator[1], gb = frame$generator[2], pr = 0, pc = 0,
                    hp = 0, wp = 0, base = 0, n_cells = as.numeric(frame$n_cells))
  }
  for (k in seq_len(nrow(levels))) {
    L <- levels[k, ]
    base_hi <- floor(L[["base"]] / 2^32)
    word <- 784 + 16 * (k - 1)
    put(word, le32(L[c("dim", "index", "c", "hp")], "u32"))
    put(word + 4, le32(L[c("ga", "gb", "pr", "pc")], "i32"))
    put(word + 8, le32(c(L[["wp"]], base_hi, L[["base"]] - base_hi * 2^32, 0), "u32"))
    put(word + 12, le32(c(sqrt(4 * pi / L[["n_cells"]]), 0, 0, 0), "f32"))
  }
  out
}

#' The textures of globe.wgsl's cell table, as the widget and hexglobe make
#' them: the table's words, the perfect hash's keys and its offsets, one
#' word each where there is no table or no hash
#' @noRd
globe_textures <- function(table) {
  one <- function(words) cpp_base64_buffer(words, "u32")
  texture <- function(binding, format, dimension, size, data) {
    list(binding = binding, format = format, dimension = dimension,
         size = as.integer(size), data = data)
  }
  keyed <- isTRUE(table$keyed)
  list(
    texture(6L, "r32uint", "2d-array", table$size %||% c(1, 1, 1),
            table$values %||% one(2^32 - 1)),
    texture(7L, "rg32uint", "2d-array", if (keyed) table$size else c(1, 1, 1),
            if (keyed) table$keys else one(c(2^32 - 1, 2^32 - 1))),
    texture(8L, "r32uint", "2d", c(table$offsets_size %||% c(1, 1), 1),
            table$offsets %||% one(0))
  )
}

#' Numbers as little-endian 32-bit floats ("f32"), signed ("i32") or unsigned
#' ("u32") integers; integers are written as their two 16-bit halves, as R's
#' integers stop short of 2^31
#' @noRd
le32 <- function(x, type) {
  x <- as.numeric(x)
  if (type == "f32") return(writeBin(x, raw(), size = 4L, endian = "little"))
  x <- x %% 2^32
  writeBin(as.integer(rbind(x %% 65536, x %/% 65536)), raw(), size = 2L,
           endian = "little")
}

#' Open a globe in headless Chrome, with its GPU, wait until it is drawn and
#' call `f(session, read)`, where `read(js)` evaluates JavaScript in the page
#' (awaiting a promise) and returns its value. A browser whose WebGPU finds
#' no graphics adapter is closed and another launched, up to `launches` in
#' all: on machines with two GPUs Chrome's GPU process sometimes starts
#' without one.
#' @noRd
globe_in_chrome <- function(widget, width, height, scale, timeout, f,
                            launches = 3) {
  dir <- tempfile("hex_globe_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  html <- file.path(dir, "globe.html")
  widget$width <- width
  widget$height <- height
  htmlwidgets::saveWidget(widget, html, selfcontained = FALSE, libdir = "lib")
  url <- paste0("file:///", normalizePath(html, winslash = "/"))

  args <- c(setdiff(chromote::default_chrome_args(), "--disable-gpu"),
            "--enable-unsafe-webgpu")
  for (launch in seq_len(launches)) {
    browser <- chromote::Chromote$new(browser = chromote::Chrome$new(args = args))
    page <- tryCatch(globe_page(browser, url, width, height, scale, timeout),
                     error = function(e) {
                       browser$close()
                       stop(e)
                     })
    if (page$state != "no-adapter" || launch == launches) break
    browser$close()
  }
  on.exit(browser$close(), add = TRUE)
  if (page$state != "drawn") {
    stop("the globe could not be drawn: ", page$message, call. = FALSE)
  }
  f(page$session, page$read)
}

#' Load the page at `url` in a new session of `browser` and wait until the
#' globe is drawn. `state` is "drawn", "no-adapter" or "failed", with the
#' widget's `message` when it is not drawn.
#' @noRd
globe_page <- function(browser, url, width, height, scale, timeout) {
  session <- browser$new_session(width = width, height = height)
  session$Emulation$setDeviceMetricsOverride(
    width = width, height = height, deviceScaleFactor = scale, mobile = FALSE
  )
  session$Page$navigate(url)

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
  if (!state %in% c("drawn", "failed")) {
    stop("the globe was not drawn within ", timeout, " seconds", call. = FALSE)
  }
  message <- ""
  if (state == "failed") {
    message <- read("document.querySelector('.hexify-globe').dataset.message")
    if (identical(read("document.querySelector('.hexify-globe').dataset.reason"),
                  "no-adapter")) {
      state <- "no-adapter"
    }
  }
  list(session = session, read = read, state = state, message = message)
}

#' Cells of points found by the globe's shader, for testing it against
#' lonlat_to_cell(): a data frame of each point's cell `id`, the `word` the
#' cells' table holds for it (2^32 - 1 for none) and the `flags` of
#' spot_record() in globe.wgsl. The points go to the GPU as 32-bit floats;
#' `xyz` are their unit vectors as the GPU reads them. `...` goes to
#' hex_globe(), such as the values and cells the table holds.
#' @noRd
globe_shader_cells <- function(grid, xyz, ..., timeout = 120) {
  widget <- hex_globe(grid, land = FALSE, face_edges = FALSE, ...)
  chunks <- split(seq_len(nrow(xyz)), ceiling(seq_len(nrow(xyz)) / 2e5))
  out <- globe_in_chrome(widget, 64, 64, 1, timeout, function(session, read) {
    lapply(chunks, function(k) {
      b64 <- cpp_base64_buffer(as.vector(t(xyz[k, , drop = FALSE])), "f32")
      read(paste0("document.querySelector('.hexify-globe').hexGlobe.locate('", b64,
                  "').then(ids => Array.from(ids).join(','))"))
    })
  })
  m <- matrix(as.numeric(unlist(strsplit(unlist(out), ",", fixed = TRUE))), nrow = 3)
  data.frame(id = m[1, ], word = m[2, ], flags = m[3, ])
}

#' An ISEA grid as the widget reads it to find each pixel's cell
#'
#' The Grid uniform (see globe_grid_uniform()) and the textures of the
#' cells' table (cpp_globe_table()), base64. Values fill a table at the
#' grid's resolution and every coarser one; cells given without values fill
#' it at the grid's resolution, to mark them drawn. With no cells given every
#' cell is drawn, and with no values either there is no table.
#' @noRd
globe_grid <- function(g, cells, values, limits, smooth = FALSE, na_fill = NULL) {
  lv <- isea_levels(g@aperture, g@resolution)
  icosa <- icosa_arg(g)
  frame <- cpp_globe_frame(icosa, lv$resolution, lv$aperture, lv$ap_seq)
  # The shader reads a point to 32-bit float precision, a few millionths of
  # a radian; past this quad side its cells are finer than that.
  if (frame$dim >= 2^24) {
    stop("hex_globe() draws ISEA grids up to aperture 3 resolution 30, ",
         "aperture 4 resolution 23 and aperture 7 resolution 16; this grid's ",
         "cells are finer than the graphics card's 32-bit floats resolve",
         call. = FALSE)
  }
  all <- is.null(cells)
  valued <- !is.null(values)
  levels <- if (valued) {
    lapply(g@resolution:0, function(r) isea_levels(g@aperture, r))
  } else {
    list(lv)
  }
  table <- NULL
  map <- if (valued) ramp_map(values, limits)
  if (all && valued) {
    check_values(values, as.numeric(frame$n_cells))
    table <- cpp_globe_table(icosa, levels, bit64::integer64(0), as.numeric(values), TRUE,
                             map)
  } else if (!all) {
    cells <- as_cell_id(cells, "cells")
    if (anyNA(cells) || any(cells < 1L | cells > frame$n_cells)) {
      stop("cells must be cell IDs of the grid, from 1 to ",
           as.character(frame$n_cells), call. = FALSE)
    }
    if (anyDuplicated(cells)) stop("cells must not repeat", call. = FALSE)
    if (valued) check_values(values, length(cells))
    table <- cpp_globe_table(icosa, levels, cells,
                             if (valued) as.numeric(values) else numeric(0), FALSE,
                             map %||% numeric(0))
  }
  list(
    uniform = cpp_base64_bytes(
      globe_grid_uniform(cpp_globe_projection(icosa), frame, table$levels, table, all,
                         valued, smooth, split_u64(frame$per_quad), na_fill, map)),
    textures = globe_textures(table),
    all = all,
    valued = valued,
    levels = table$levels,
    given = table$given,
    keyed = isTRUE(table$keyed),
    ramp_map = map
  )
}

#' Non-negative integer64 values as their high and low 32-bit words,
#' interleaved, each as a double
#' @noRd
split_u64 <- function(x) {
  word <- bit64::as.integer64(2^32)
  hi <- x %/% word
  as.vector(rbind(as.numeric(hi), as.numeric(x - hi * word)))
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
land_surface_paths <- function(land, max_angle, icosa) {
  rings <- unlist(sfc_polygons(land), recursive = FALSE)
  path <- rep(seq_along(rings), vapply(rings, nrow, integer(1)))
  ll <- do.call(rbind, rings)
  cpp_sphere_paths_on_faces(icosa, ll[, 1], ll[, 2], path, max_angle)
}

#' The solid's edges placed on the faces, one path per edge
#' @noRd
edge_surface_paths <- function(max_angle, icosa) {
  solid <- icosa_solid(icosa)
  V <- solid$vertices
  ends <- as.vector(t(solid$edges[, c("v1", "v2")]))
  ll <- vec_lonlat(V[ends, ])
  cpp_sphere_paths_on_faces(icosa, ll[, 1], ll[, 2],
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
#' `clamp((v - m[1]) * m[2] + m[3], 0, 1)`, the limits at its two ends, or
#' halfway when they coincide
#' @noRd
ramp_map <- function(values, limits) {
  if (is.null(limits)) {
    limits <- suppressWarnings(c(min(values, na.rm = TRUE), max(values, na.rm = TRUE)))
    if (identical(limits, c(Inf, -Inf))) limits <- c(0, 1)
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

#' A mesh as the widget reads it: per vertex its solid and sphere
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

#' Paths as the widget reads them: per point its solid and sphere
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
