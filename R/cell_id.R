#' ISEA cell IDs as exact 64-bit integers
#'
#' An ISEA cell ID numbers a cell from 1 to the grid's cell count, which passes
#' 2^53 at aperture 4 resolution 25 and aperture 7 resolution 18. A double
#' holds whole numbers exactly only below 2^53, so hexify returns ISEA cell
#' IDs as `bit64::integer64` and reads them back as given here.
#'
#' Accepted: integer64; character strings of digits, read exactly at any size;
#' whole numbers stored as double or integer below 2^53. `NA` stays `NA`.
#' @param x Cell IDs
#' @param arg The argument's name, for the error message
#' @return An integer64 vector of the same length
#' @noRd
as_cell_id <- function(x, arg = "cell_id") {
  if (bit64::is.integer64(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    out <- suppressWarnings(bit64::as.integer64(x))
    bad <- is.na(out) & !is.na(x)
    if (any(bad)) {
      stop(sprintf("`%s` must hold whole-number cell IDs; got \"%s\"",
                   arg, x[which(bad)[1L]]), call. = FALSE)
    }
    return(out)
  }
  if (is.logical(x) && all(is.na(x))) {
    return(bit64::as.integer64(rep(NA_real_, length(x))))
  }
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must hold cell IDs; got an object of class %s",
                 arg, paste(class(x), collapse = "/")), call. = FALSE)
  }
  x <- as.numeric(x)
  ok <- is.na(x)
  x_ok <- x[!ok]
  if (any(!is.finite(x_ok) | x_ok != floor(x_ok))) {
    if (all(x_ok > 0 & x_ok < 1)) {
      # The 64 bits of an integer64 read as a double are a number between 0
      # and 1 for every ID below 2^62
      stop(sprintf(paste0(
        "`%s` holds numbers between 0 and 1: integer64 cell IDs that lost ",
        "their class, as unlist(), sapply(), ifelse() and for loops drop it. ",
        "Combine a list of cell IDs with do.call(c, ids)"), arg), call. = FALSE)
    }
    stop(sprintf("`%s` must hold whole-number cell IDs", arg), call. = FALSE)
  }
  if (any(abs(x_ok) >= 2^53)) {
    stop(sprintf(paste0(
      "`%s` holds a cell ID of 2^53 or more as a double, which cannot be ",
      "read exactly; pass such IDs as bit64::integer64 or as character"), arg),
      call. = FALSE)
  }
  bit64::as.integer64(x)
}

#' Concatenate a list of cell ID vectors into one integer64 vector
#'
#' `unlist()` drops the integer64 class and returns the raw bits as doubles,
#' so the class is restored on the concatenated slots.
#' @param x List of cell ID vectors
#' @noRd
cell_id_unlist <- function(x) {
  out <- unlist(lapply(x, as_cell_id), use.names = FALSE)
  if (is.null(out)) out <- numeric(0)
  class(out) <- "integer64"
  out
}

#' Finest resolution of an ISEA aperture on a solid whose cell IDs fit in a
#' signed 64-bit integer
#'
#' Aperture 3 reaches MAX_RESOLUTION; aperture 4 stops at 29 and aperture 7 at
#' 21 on the icosahedron, where their cell counts pass 2^63 - 1. Defined for a
#' pure aperture or a family; a per-level spelling has one resolution.
#' @param aperture Aperture spelling, as a grid stores it
#' @param polyhedron The solid the grid is built on
#' @noRd
isea_max_resolution <- function(aperture, polyhedron = "icosahedron") {
  fits <- vapply(seq.int(MIN_RESOLUTION, MAX_RESOLUTION), function(r) {
    !is.na(isea_cell_count(aperture, r, polyhedron))
  }, logical(1))
  MIN_RESOLUTION + sum(cumprod(fits)) - 1L
}

#' Stop unless an ISEA grid's cell IDs fit in a signed 64-bit integer
#' @param aperture Aperture spelling, as a grid stores it
#' @param resolution Resolution
#' @param polyhedron The solid the grid is built on
#' @noRd
check_isea_resolution <- function(aperture, resolution,
                                  polyhedron = "icosahedron") {
  if (!is.na(isea_cell_count(aperture, resolution, polyhedron))) {
    return(invisible(TRUE))
  }
  msg <- sprintf(paste0(
    "Aperture %s on the %s has more cells at resolution %d than 64-bit cell ",
    "IDs can number (2^63 - 1)"), aperture, polyhedron, resolution)
  if (!is_per_level_aperture(aperture)) {
    msg <- sprintf("%s; its finest resolution is %d", msg,
                   isea_max_resolution(aperture, polyhedron))
  }
  stop(msg, call. = FALSE)
}
