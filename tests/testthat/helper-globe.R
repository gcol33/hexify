gpu_browser <- new.env()

# Chrome is found on most machines, WebGPU only on some: the probe draws a
# small globe once per run, and the tests skip where it cannot be drawn.
skip_without_gpu_browser <- function() {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if(is.null(suppressMessages(chromote::find_chrome())), "no Chromium browser")
  if (is.null(gpu_browser$drawn)) {
    globe <- hex_globe(hex_grid(resolution = 1, aperture = 3), land = FALSE)
    gpu_browser$drawn <- tryCatch(
      hexify:::globe_in_chrome(globe, 64, 64, 1, 30, function(session, read) TRUE),
      error = function(e) conditionMessage(e)
    )
  }
  skip_if(!isTRUE(gpu_browser$drawn),
          paste("Chrome cannot draw the globe:", gpu_browser$drawn))
}

skip_without_wgpu <- function() {
  skip_on_cran()
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("hexglobe")
  skip_if(inherits(try(hexglobe::gpu_adapter(), silent = TRUE), "try-error"),
          "no graphics adapter for wgpu")
}
