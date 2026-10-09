# ogc_zirs.R
# The zone identifiers of the ISEA3H, ISEA7H, IVEA3H and IVEA7H DGGRS
# definitions registered with OGC (OGC API - DGGS, DGGRS definition examples),
# as DGGAL writes them: textZIRS ("levelRootFaceHexRowMajorSubZone") and
# uint64ZIRS.
#
# A zone is named by a root rhombus, a sub-rhombus within it and a letter. The
# root rhombi 0-9 are hexify's diamond quads, in DGGAL's staircase order: an
# even root 2k is quad k + 1 (the upper diamonds), an odd root 2k + 1 is quad
# k + 6; roots 10 and 11 are the vertex quads 0 and 11, the "North" and
# "South" poles. At an even level a quad holds P x P sub-rhombi, P = 3^(L/2)
# or 7^(L/2), whose top-left corners are the cell centres: in the quad plane
# the corner of row r and column c lies at (r e1 + c e2) / P, e1 = (1, 0),
# e2 = (-1/2, sqrt(3)/2), and its index counts row-major, r P + c. Letter A
# names that cell. At an odd level:
#   - aperture 3: the cells lie at the sub-rhombi of level L - 1, at its
#     corner (B) and at the centroids of its two triangles, (1/3, 2/3) in
#     (row, column) units (C) and (2/3, 1/3) (D);
#   - aperture 7: the seven children of the level L - 1 cell at a corner take
#     its sub-rhombus, the centred child B and the six around it C to H,
#     counter-clockwise, F lying from the parent at 19.11 degrees, arctan of
#     sqrt(3) / 5, in the quad plane. A pentagon at a quad's origin has five
#     such children: on an upper quad F, G, C, D, E from that direction on,
#     on a lower one E, F, G, C, D; a pole's children are C to G in the order
#     of their quads, 1 to 5 around the north pole and 10 down to 6 around
#     the south pole. The children of a hexagon on the first row of an upper
#     quad, or on the first column of a lower one, which DGGAL reads in the
#     neighbouring quad's frame, take the letter a sixth of a turn further
#     counter-clockwise (upper) or clockwise (lower).

#' Bits of the uint64ZIRS fields, from the least significant: letter,
#' sub-rhombus index, root rhombus, level field
#' @noRd
ZIRS_BITS <- list(`3` = c(sub = 2, ix = 51, root = 4, level = 5),
                  `7` = c(sub = 3, ix = 51, root = 4, level = 4))

#' Direction from an aperture-7 parent to its child F in the quad plane
#' @noRd
ZIRS_AP7_F_ANGLE <- atan(sqrt(3) / 5)

#' Stop unless a grid has OGC zone identifiers
#' @noRd
check_zirs_grid <- function(g) {
  if (is_h3_grid(g) || grid_polyhedron(g) != "icosahedron" ||
      !g@aperture %in% c("3", "7")) {
    stop("textZIRS and uint64ZIRS name the zones of aperture-3 and aperture-7 ",
         "grids on the icosahedron (ISEA3H, ISEA7H, IVEA3H, IVEA7H)", call. = FALSE)
  }
  if (g@aperture == "7" && g@resolution > 19L) {
    stop("textZIRS and uint64ZIRS name aperture-7 zones up to level 19; ",
         "this grid has resolution ", g@resolution, call. = FALSE)
  }
  invisible(g)
}

#' Root rhombus of a quad, and back
#' @noRd
quad_root <- function(quad) {
  ifelse(quad == 0L, 10L, ifelse(quad == 11L, 11L,
         ifelse(quad <= 5L, 2L * (quad - 1L), 2L * (quad - 6L) + 1L)))
}

#' @rdname quad_root
#' @noRd
root_quad <- function(root) {
  ifelse(root == 10L, 0L, ifelse(root == 11L, 11L,
         ifelse(root %% 2L == 0L, root %/% 2L + 1L, (root - 1L) %/% 2L + 6L)))
}

#' Sub-rhombus rows and columns per quad side at level L
#' @noRd
zirs_side <- function(aperture, level) as.numeric(aperture)^(level %/% 2L)

#' Quad-plane centres of cells in row and column units of P per quad side
#' @noRd
zirs_quad_uv <- function(cell_id, g, P) {
  q <- cpp_cell_to_quad_xy(icosa_arg(g), as_cell_id(cell_id), g@resolution,
                           as.integer(g@aperture))
  data.frame(quad = q$quad, x = q$quad_x, y = q$quad_y,
             u = (q$quad_x + q$quad_y / sqrt(3)) * P,
             v = q$quad_y * 2 / sqrt(3) * P)
}

#' The fields of the OGC zone identifier of cells
#'
#' @param cell_id Cell IDs of grid `g`
#' @param g HexGridInfo object passing check_zirs_grid()
#' @return Data frame: level, root, ix (sub-rhombus index, a whole double),
#'   sub (letter, 0 for A)
#' @noRd
zirs_fields <- function(cell_id, g) {
  cell_id <- as_cell_id(cell_id)
  L <- g@resolution
  ap <- as.integer(g@aperture)
  if (ap == 7L && L %% 2L == 1L) return(zirs_fields_ap7_odd(cell_id, g))

  P <- zirs_side(ap, L)
  q <- zirs_quad_uv(cell_id, g, P)
  pole <- q$quad %in% c(0L, 11L)
  if (L %% 2L == 0L) {
    row <- round(q$u)
    col <- round(q$v)
    sub <- rep(0L, length(cell_id))
  } else {
    row <- floor(q$u + 1e-9)
    col <- floor(q$v + 1e-9)
    fu <- q$u - row
    sub <- ifelse(fu < 1 / 6, 1L, ifelse(fu < 1 / 2, 2L, 3L))
  }
  data.frame(level = rep(L, length(cell_id)), root = quad_root(q$quad),
             ix = ifelse(pole, 0, row * P + col), sub = ifelse(pole, ifelse(L %% 2L == 0L, 0L, 1L), sub))
}

#' zirs_fields() at an odd aperture-7 level: the parent's sub-rhombus and the
#' child's letter
#' @noRd
zirs_fields_ap7_odd <- function(cell_id, g) {
  L <- g@resolution
  gp <- grid_at_resolution(g, L - 1L)
  icosa <- icosa_arg(g)
  parent <- get_parent(cell_id, g)
  pf <- zirs_fields(parent, gp)

  Pp <- zirs_side(7L, L - 1L)
  pq <- zirs_quad_uv(parent, gp, Pp)
  cq <- zirs_quad_uv(cell_id, g, Pp)
  centred <- pq$quad == cq$quad & abs(pq$x - cq$x) < 1e-9 & abs(pq$y - cq$y) < 1e-9
  sub <- rep(1L, length(cell_id))

  # Children of a pole: C onwards in the order of their quads
  north <- !centred & pf$root == 10L
  south <- !centred & pf$root == 11L
  sub[north] <- 2L + cq$quad[north] - 1L
  sub[south] <- 2L + 10L - cq$quad[south]

  # Other vertex children: counter-clockwise turn from the parent's F
  # direction, read on the sphere around the parent
  around <- which(!centred & pf$root < 10L)
  if (length(around)) {
    eps <- 1e-3 / Pp
    a <- cpp_quad_xy_to_lonlat(icosa, pq$quad[around],
                               pq$x[around] + eps * cos(ZIRS_AP7_F_ANGLE),
                               pq$y[around] + eps * sin(ZIRS_AP7_F_ANGLE))
    pc <- cell_to_lonlat(parent[around], gp)
    cc <- cell_to_lonlat(cell_id[around], g)
    P0 <- unit_vec(pc$lon_deg, pc$lat_deg, icosa)
    A <- unit_vec(a$lon_deg, a$lat_deg, icosa)
    C <- unit_vec(cc$lon_deg, cc$lat_deg, icosa)
    t1 <- A - rowSums(A * P0) * P0
    t1 <- t1 / sqrt(rowSums(t1^2))
    t2 <- cross3(P0, t1)
    d <- C - rowSums(C * P0) * P0
    turn <- atan2(rowSums(d * t2), rowSums(d * t1)) * 180 / pi
    pentagon <- pf$ix[around] == 0
    upper <- pf$root[around] %% 2L == 0L
    # DGGAL reads the children of a parent on the first row of an upper quad,
    # or the first column of a lower one, a sixth of a turn further round
    edge <- ifelse(upper & pf$ix[around] %/% Pp == 0, 1,
                   ifelse(!upper & pf$ix[around] %% Pp == 0, -1, 0))
    k6 <- (round(turn / 60) + edge) %% 6
    k5 <- round(turn / 72) %% 5
    sub[around] <- ifelse(!pentagon, c(5L, 6L, 7L, 2L, 3L, 4L)[k6 + 1L],
                          ifelse(upper, c(5L, 6L, 2L, 3L, 4L)[k5 + 1L],
                                 c(4L, 5L, 6L, 2L, 3L)[k5 + 1L]))
  }
  data.frame(level = rep(L, length(cell_id)), root = pf$root, ix = pf$ix, sub = sub)
}

#' Cells of the grid named by OGC zone identifier fields; NA where the fields
#' name no zone of it
#' @noRd
zirs_cells <- function(f, g) {
  L <- g@resolution
  ap <- as.integer(g@aperture)
  n <- nrow(f)
  out <- as_cell_id(rep(NA, n))
  P <- zirs_side(ap, L)
  max_sub <- if (L %% 2L == 0L) 0L else if (ap == 3L) 3L else 7L
  ok <- !is.na(f$root) & !is.na(f$ix) & !is.na(f$sub) & f$level == L &
    f$root <= 11L & f$sub >= as.integer(L %% 2L) & f$sub <= max_sub &
    f$ix < ifelse(f$root >= 10L, 1, P^2)
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  if (ap == 7L && L %% 2L == 1L) {
    gp <- grid_at_resolution(g, L - 1L)
    parent <- zirs_cells(transform(f, level = L - 1L, sub = 0L), gp)
    ok <- ok & !is.na(parent)
    pos <- which(ok)
    if (length(pos)) {
      parents <- unique(parent[pos])
      kids <- get_children(parents, gp)
      all_kids <- cell_id_unlist(kids)
      kid_parent <- rep(parents, lengths(kids))
      kf <- zirs_fields(all_kids, g)
      key <- paste(as.character(kid_parent), kf$sub)
      out[pos] <- all_kids[match(paste(as.character(parent[pos]), f$sub[pos]), key)]
    }
  } else {
    pole <- f$root %in% c(10L, 11L)
    row <- f$ix %/% P
    col <- f$ix %% P
    du <- c(0, 0, 1 / 3, 2 / 3)[f$sub + 1L]
    dv <- c(0, 0, 2 / 3, 1 / 3)[f$sub + 1L]
    u <- (row + du) / P
    v <- (col + dv) / P
    y <- v * sqrt(3) / 2
    x <- u - y / sqrt(3)
    quad <- root_quad(f$root)
    at <- which(ok & !pole)
    if (length(at)) {
      out[at] <- cpp_quad_xy_to_cell(icosa_arg(g), quad[at], x[at], y[at], L, ap)
    }
    out[ok & f$root == 10L] <- as_cell_id(1)
    out[ok & f$root == 11L] <- as_cell_id(grid_n_cells(g))
  }
  # A string names a zone only in its own spelling
  back <- which(!is.na(out))
  if (length(back)) {
    bf <- zirs_fields(out[back], g)
    same <- bf$root == f$root[back] & bf$ix == f$ix[back] & bf$sub == f$sub[back]
    out[back[!same]] <- NA
  }
  out
}

#' Hexadecimal digits of whole doubles below 2^53
#' @noRd
hex_digits <- function(x) {
  vapply(x, function(v) {
    if (v == 0) return("0")
    d <- character(0)
    while (v > 0) {
      d <- c(c(0:9, LETTERS[1:6])[v %% 16 + 1], d)
      v <- v %/% 16
    }
    paste(d, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

#' Whole doubles of hexadecimal digit strings, NA for anything else
#' @noRd
parse_hex <- function(s) {
  vapply(s, function(h) {
    if (is.na(h) || !grepl("^[0-9A-F]{1,13}$", h)) return(NA_real_)
    v <- 0
    for (ch in strsplit(h, "")[[1]]) v <- v * 16 + match(ch, c(0:9, LETTERS[1:6])) - 1
    v
  }, numeric(1), USE.NAMES = FALSE)
}

#' textZIRS strings of zone identifier fields
#' @noRd
zirs_text <- function(f, aperture) {
  lev <- if (aperture == "3") f$level %/% 2L else f$level
  sprintf("%s%X-%s-%s", LETTERS[lev + 1L], as.integer(f$root), hex_digits(f$ix),
          LETTERS[f$sub + 1L])
}

#' Zone identifier fields of textZIRS strings; NA fields where a string is
#' not one
#' @noRd
parse_zirs_text <- function(text, aperture) {
  m <- regmatches(text, regexec("^([A-Z])([0-9AB])-([0-9A-F]+)-([A-H])$", text))
  field <- function(k) vapply(m, function(x) if (length(x)) x[k] else NA_character_,
                              character(1))
  lev <- match(field(2), LETTERS) - 1L
  sub <- match(field(5), LETTERS) - 1L
  data.frame(level = if (aperture == "3") 2L * lev + (sub > 0L) else lev,
             root = strtoi(field(3), 16L), ix = parse_hex(field(4)), sub = sub)
}

#' uint64ZIRS values of zone identifier fields, as integer64
#' @noRd
zirs_uint64 <- function(f, aperture) {
  b <- ZIRS_BITS[[aperture]]
  lev <- if (aperture == "3") f$level %/% 2L else f$level %/% 2L
  shift <- function(k) bit64::as.integer64(2^k)
  bit64::as.integer64(lev) * shift(b[["sub"]] + b[["ix"]] + b[["root"]]) +
    bit64::as.integer64(f$root) * shift(b[["sub"]] + b[["ix"]]) +
    bit64::as.integer64(f$ix) * shift(b[["sub"]]) +
    bit64::as.integer64(f$sub)
}

#' Zone identifier fields of uint64ZIRS values
#' @noRd
parse_zirs_uint64 <- function(x, aperture) {
  b <- ZIRS_BITS[[aperture]]
  x <- bit64::as.integer64(x)
  field <- function(lo, width) {
    as.numeric((x %/% bit64::as.integer64(2^lo)) %% bit64::as.integer64(2^width))
  }
  sub <- as.integer(field(0, b[["sub"]]))
  lev <- as.integer(field(b[["sub"]] + b[["ix"]] + b[["root"]], b[["level"]]))
  bad <- is.na(x) | x < 0 |
    x >= bit64::as.integer64(2^(b[["sub"]] + b[["ix"]] + b[["root"]] + b[["level"]]))
  out <- data.frame(level = 2L * lev + (sub > 0L),
                    root = as.integer(field(b[["sub"]] + b[["ix"]], b[["root"]])),
                    ix = field(b[["sub"]], b[["ix"]]), sub = sub)
  out[bad, ] <- NA
  out
}
