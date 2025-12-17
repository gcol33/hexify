# Projection Geometry Diagrams
# Mathematical illustrations for Lambert azimuthal equal-area projection
# and ISEA icosahedral face projection

# Color palette
COL_GRAY30 <- gray(0.3)
COL_GRAY50 <- gray(0.5)
COL_RED <- "#E63946"
COL_BLUE <- "#457B9D"

# Helper functions --------------------------------------------------------

#' Draw a circular arc
#' @param x0,y0 Center coordinates
#' @param r Radius
#' @param theta1,theta2 Start and end angles in radians
#' @param n Number of segments
#' @param ... Additional graphics parameters
draw_arc <- function(x0, y0, r, theta1, theta2, n = 100, ...) {
  theta <- seq(theta1, theta2, length.out = n)
  x <- x0 + r * cos(theta)
  y <- y0 + r * sin(theta)
  lines(x, y, ...)
}

#' Draw a circle
#' @param x0,y0 Center coordinates
#' @param r Radius
#' @param n Number of segments
#' @param ... Additional graphics parameters
draw_circle <- function(x0, y0, r, n = 200, ...) {
  theta <- seq(0, 2 * pi, length.out = n)
  x <- x0 + r * cos(theta)
  y <- y0 + r * sin(theta)
  lines(x, y, ...)
}

#' Draw a filled circle
#' @param x0,y0 Center coordinates
#' @param r Radius
#' @param n Number of segments
#' @param ... Additional graphics parameters
draw_filled_circle <- function(x0, y0, r, n = 200, ...) {
  theta <- seq(0, 2 * pi, length.out = n)
  x <- x0 + r * cos(theta)
  y <- y0 + r * sin(theta)
  polygon(x, y, ...)
}

#' Draw a point with label
#' @param x,y Point coordinates
#' @param label Text label
#' @param pch Point character
#' @param col Point color
#' @param offset Label offset as c(dx, dy)
#' @param cex Character expansion
draw_labeled_point <- function(x, y, label, pch = 16, col = COL_GRAY30,
                                offset = c(0, 0.15), cex = 1) {
  points(x, y, pch = pch, col = col, cex = cex)
  text(x + offset[1], y + offset[2], label, col = col, cex = 1.2)
}

#' Draw an arrow with annotations
#' @param x0,y0 Start coordinates
#' @param x1,y1 End coordinates
#' @param col Arrow color
#' @param lwd Line width
#' @param ... Additional arrow parameters
draw_arrow <- function(x0, y0, x1, y1, col = COL_GRAY30, lwd = 1.5, ...) {
  arrows(x0, y0, x1, y1, col = col, lwd = lwd, length = 0.1, ...)
}

# Figure 1: Lambert Chord Distance ----------------------------------------

#' Create Lambert chord distance diagram
#' Shows cross-section of sphere with tangent plane projection
#'
#' @param output_path Base path for output files (without extension)
#' @param phi Angular distance in radians (default: pi/4, or 45 degrees)
create_lambert_chord_diagram <- function(output_path = NULL, phi = pi/4) {

  # Set up output
  if (!is.null(output_path)) {
    svg(paste0(output_path, ".svg"), width = 8, height = 8)
  }

  par(mar = c(1, 1, 2, 1))
  plot.new()
  plot.window(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1)
  title("Lambert Azimuthal Equal-Area Projection\nChord Distance Geometry",
        cex.main = 1.3, font.main = 1)

  # Constants
  R <- 1  # Unit sphere radius
  x0 <- 0  # Center x
  y0 <- 0  # Center y

  # Point S: tangent point (top of sphere)
  S_x <- 0
  S_y <- R

  # Point P: point on sphere at angular distance phi from S
  # In cross-section, P is at angle phi from vertical
  P_x <- R * sin(phi)
  P_y <- R * cos(phi)

  # Point P': projection of P onto tangent plane
  # Chord distance formula: d = 2R * sin(phi/2)
  # P' lies on tangent plane (y = R) along ray from O through P
  # The projection extends the line OP to intersect the tangent plane
  # Scale factor: R / (R * cos(phi)) = 1/cos(phi)
  scale <- 1 / cos(phi)
  P_prime_x <- P_x * scale
  P_prime_y <- R  # On tangent plane

  # Draw sphere (circle in cross-section)
  draw_circle(x0, y0, R, col = COL_GRAY30, lwd = 2)

  # Draw tangent plane at top
  segments(-1.5, R, 1.5, R, col = COL_GRAY30, lwd = 2, lty = 1)
  text(-1.35, R + 0.1, "Tangent plane", col = COL_GRAY30, adj = c(0, 0), cex = 1.1)

  # Draw radius to S (vertical)
  segments(x0, y0, S_x, S_y, col = COL_GRAY50, lwd = 1, lty = 2)

  # Draw radius to P
  segments(x0, y0, P_x, P_y, col = COL_GRAY50, lwd = 1, lty = 2)

  # Draw projection ray from O through P to P'
  draw_arrow(x0, y0, P_prime_x, P_prime_y, col = COL_BLUE, lwd = 2, lty = 1)

  # Draw chord distance d from S to P (along tangent plane)
  segments(S_x, S_y, P_prime_x, P_prime_y, col = COL_RED, lwd = 3)

  # Draw arc on sphere from S to P (showing angular distance)
  draw_arc(x0, y0, R, pi/2 - phi, pi/2, col = COL_GRAY50, lwd = 2)

  # Draw angular distance phi arc (small arc near center)
  arc_r <- 0.3
  draw_arc(x0, y0, arc_r, pi/2 - phi, pi/2, col = COL_GRAY50, lwd = 1.5)

  # Label angular distance
  phi_label_x <- arc_r * cos(pi/2 - phi/2) + 0.15
  phi_label_y <- arc_r * sin(pi/2 - phi/2)
  text(phi_label_x, phi_label_y, expression(phi), col = COL_GRAY50, cex = 1.3)

  # Draw points
  draw_labeled_point(x0, y0, "O", col = COL_GRAY30, offset = c(-0.15, -0.15))
  draw_labeled_point(S_x, S_y, "S", col = COL_GRAY30, offset = c(-0.15, 0.15))
  draw_labeled_point(P_x, P_y, "P", pch = 16, col = COL_RED, offset = c(0.15, 0))
  draw_labeled_point(P_prime_x, P_prime_y, "P'", pch = 16, col = COL_BLUE,
                     offset = c(0, 0.18))

  # Label chord distance d
  d_mid_x <- (S_x + P_prime_x) / 2
  d_mid_y <- S_y - 0.18
  text(d_mid_x, d_mid_y, "d", col = COL_RED, cex = 1.3, font = 2)

  # Add formula annotation
  d_val <- 2 * R * sin(phi/2)
  formula_text <- bquote(d == 2*R*sin(phi/2))
  text(0, -1.3, formula_text, col = COL_GRAY30, cex = 1.2)

  # Add grid for reference
  abline(h = 0, v = 0, col = COL_GRAY50, lwd = 0.5, lty = 3)

  if (!is.null(output_path)) {
    dev.off()

    # Also create PNG version
    png(paste0(output_path, ".png"), width = 800, height = 800, res = 100)
    par(mar = c(1, 1, 2, 1))
    plot.new()
    plot.window(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1)
    title("Lambert Azimuthal Equal-Area Projection\nChord Distance Geometry",
          cex.main = 1.3, font.main = 1)

    draw_circle(x0, y0, R, col = COL_GRAY30, lwd = 2)
    segments(-1.5, R, 1.5, R, col = COL_GRAY30, lwd = 2, lty = 1)
    text(-1.35, R + 0.1, "Tangent plane", col = COL_GRAY30, adj = c(0, 0), cex = 1.1)
    segments(x0, y0, S_x, S_y, col = COL_GRAY50, lwd = 1, lty = 2)
    segments(x0, y0, P_x, P_y, col = COL_GRAY50, lwd = 1, lty = 2)
    draw_arrow(x0, y0, P_prime_x, P_prime_y, col = COL_BLUE, lwd = 2, lty = 1)
    segments(S_x, S_y, P_prime_x, P_prime_y, col = COL_RED, lwd = 3)
    draw_arc(x0, y0, R, pi/2 - phi, pi/2, col = COL_GRAY50, lwd = 2)
    draw_arc(x0, y0, arc_r, pi/2 - phi, pi/2, col = COL_GRAY50, lwd = 1.5)
    text(phi_label_x, phi_label_y, expression(phi), col = COL_GRAY50, cex = 1.3)
    draw_labeled_point(x0, y0, "O", col = COL_GRAY30, offset = c(-0.15, -0.15))
    draw_labeled_point(S_x, S_y, "S", col = COL_GRAY30, offset = c(-0.15, 0.15))
    draw_labeled_point(P_x, P_y, "P", pch = 16, col = COL_RED, offset = c(0.15, 0))
    draw_labeled_point(P_prime_x, P_prime_y, "P'", pch = 16, col = COL_BLUE,
                       offset = c(0, 0.18))
    text(d_mid_x, d_mid_y, "d", col = COL_RED, cex = 1.3, font = 2)
    text(0, -1.3, formula_text, col = COL_GRAY30, cex = 1.2)
    abline(h = 0, v = 0, col = COL_GRAY50, lwd = 0.5, lty = 3)

    dev.off()

    message("Lambert chord diagram saved to:\n  ",
            normalizePath(paste0(output_path, ".svg")), "\n  ",
            normalizePath(paste0(output_path, ".png")))
  }

  invisible(NULL)
}

# Figure 2: Icosahedron Face Projection -----------------------------------

#' Create icosahedron face projection concept diagram
#' Shows simplified view of one triangular face with sphere arc above
#'
#' @param output_path Base path for output files (without extension)
create_icosa_face_projection_diagram <- function(output_path = NULL) {

  # Set up output
  if (!is.null(output_path)) {
    svg(paste0(output_path, ".svg"), width = 8, height = 8)
  }

  par(mar = c(1, 1, 2, 1))
  plot.new()
  plot.window(xlim = c(-2, 2), ylim = c(-1.5, 2.5), asp = 1)
  title("Icosahedron Face Projection\nSphere to Tangent Plane",
        cex.main = 1.3, font.main = 1)

  # Triangular face vertices (equilateral triangle)
  # Height of equilateral triangle with side a: h = (sqrt(3)/2) * a
  side_length <- 3
  height <- (sqrt(3)/2) * side_length

  # Center the triangle at origin, with base horizontal
  v1_x <- -side_length/2
  v1_y <- -height/3

  v2_x <- side_length/2
  v2_y <- -height/3

  v3_x <- 0
  v3_y <- 2*height/3

  # Face center (centroid)
  fc_x <- (v1_x + v2_x + v3_x) / 3
  fc_y <- (v1_y + v2_y + v3_y) / 3

  # Draw triangular face (tangent plane)
  polygon(c(v1_x, v2_x, v3_x), c(v1_y, v2_y, v3_y),
          col = rgb(0.9, 0.9, 0.9, 0.5), border = COL_GRAY30, lwd = 3)

  # Draw face center
  draw_labeled_point(fc_x, fc_y, "Face center", pch = 16, col = COL_GRAY30,
                     offset = c(0, -0.25))

  # Draw sphere arc above the face
  # Arc should curve above the triangle
  arc_height <- 1.5
  arc_center_y <- fc_y - 2  # Center below to create upward arc
  arc_radius <- arc_height + 2

  # Calculate angles for the arc that spans roughly the triangle width
  theta_range <- 0.6  # radians on each side
  theta_mid <- pi/2  # Top of circle

  # Draw the sphere arc
  draw_arc(fc_x, arc_center_y, arc_radius,
           theta_mid - theta_range, theta_mid + theta_range,
           col = COL_GRAY30, lwd = 3, n = 100)

  # Mark a point on the sphere arc
  P_theta <- theta_mid - theta_range/3
  P_sphere_x <- fc_x + arc_radius * cos(P_theta)
  P_sphere_y <- arc_center_y + arc_radius * sin(P_theta)

  # Project this point to the face (perpendicular to face plane)
  # For simplicity, project vertically down
  P_proj_x <- P_sphere_x
  P_proj_y <- fc_y + (P_sphere_x - fc_x) * tan(pi/6)  # approximate on face

  # Ensure projected point is within triangle bounds
  if (P_proj_y < v1_y) P_proj_y <- v1_y + 0.2

  # Draw projection arrow
  draw_arrow(P_sphere_x, P_sphere_y, P_proj_x, P_proj_y,
             col = COL_BLUE, lwd = 2.5)

  # Draw points
  draw_labeled_point(P_sphere_x, P_sphere_y, "P", pch = 16, col = COL_RED,
                     offset = c(0.2, 0.2))
  draw_labeled_point(P_proj_x, P_proj_y, "P'", pch = 16, col = COL_BLUE,
                     offset = c(0, -0.25))

  # Add labels for vertices
  text(v1_x - 0.15, v1_y - 0.2, "V1", col = COL_GRAY50, cex = 1.1)
  text(v2_x + 0.15, v2_y - 0.2, "V2", col = COL_GRAY50, cex = 1.1)
  text(v3_x, v3_y + 0.2, "V3", col = COL_GRAY50, cex = 1.1)

  # Add annotation for sphere surface
  text(fc_x - 1.2, P_sphere_y + 0.3, "Sphere surface",
       col = COL_GRAY30, cex = 1.1, adj = c(0, 0.5))

  # Add annotation for tangent plane
  text(v2_x + 0.3, v2_y + 0.5, "Tangent plane\n(icosa face)",
       col = COL_GRAY30, cex = 1.1, adj = c(0, 0.5))

  # Draw additional projection arrows to show the concept
  P2_theta <- theta_mid + theta_range/3
  P2_sphere_x <- fc_x + arc_radius * cos(P2_theta)
  P2_sphere_y <- arc_center_y + arc_radius * sin(P2_theta)
  P2_proj_x <- P2_sphere_x
  P2_proj_y <- fc_y - (P2_sphere_x - fc_x) * tan(pi/6)
  if (P2_proj_y < v1_y) P2_proj_y <- v1_y + 0.2

  draw_arrow(P2_sphere_x, P2_sphere_y, P2_proj_x, P2_proj_y,
             col = COL_BLUE, lwd = 1.5, lty = 2)

  # Add scale indicator
  segments(-1.8, -1.2, -1.8, -1.2 + 1, col = COL_GRAY50, lwd = 1.5)
  text(-1.8, -1.2 - 0.15, "0", col = COL_GRAY50, cex = 0.9)
  text(-1.8, -1.2 + 1 + 0.15, "R", col = COL_GRAY50, cex = 0.9)
  text(-1.4, -1.2 + 0.5, "Scale", col = COL_GRAY50, cex = 0.9, srt = 90)

  if (!is.null(output_path)) {
    dev.off()

    # Also create PNG version
    png(paste0(output_path, ".png"), width = 800, height = 800, res = 100)
    par(mar = c(1, 1, 2, 1))
    plot.new()
    plot.window(xlim = c(-2, 2), ylim = c(-1.5, 2.5), asp = 1)
    title("Icosahedron Face Projection\nSphere to Tangent Plane",
          cex.main = 1.3, font.main = 1)

    polygon(c(v1_x, v2_x, v3_x), c(v1_y, v2_y, v3_y),
            col = rgb(0.9, 0.9, 0.9, 0.5), border = COL_GRAY30, lwd = 3)
    draw_labeled_point(fc_x, fc_y, "Face center", pch = 16, col = COL_GRAY30,
                       offset = c(0, -0.25))
    draw_arc(fc_x, arc_center_y, arc_radius,
             theta_mid - theta_range, theta_mid + theta_range,
             col = COL_GRAY30, lwd = 3, n = 100)
    draw_arrow(P_sphere_x, P_sphere_y, P_proj_x, P_proj_y,
               col = COL_BLUE, lwd = 2.5)
    draw_labeled_point(P_sphere_x, P_sphere_y, "P", pch = 16, col = COL_RED,
                       offset = c(0.2, 0.2))
    draw_labeled_point(P_proj_x, P_proj_y, "P'", pch = 16, col = COL_BLUE,
                       offset = c(0, -0.25))
    text(v1_x - 0.15, v1_y - 0.2, "V1", col = COL_GRAY50, cex = 1.1)
    text(v2_x + 0.15, v2_y - 0.2, "V2", col = COL_GRAY50, cex = 1.1)
    text(v3_x, v3_y + 0.2, "V3", col = COL_GRAY50, cex = 1.1)
    text(fc_x - 1.2, P_sphere_y + 0.3, "Sphere surface",
         col = COL_GRAY30, cex = 1.1, adj = c(0, 0.5))
    text(v2_x + 0.3, v2_y + 0.5, "Tangent plane\n(icosa face)",
         col = COL_GRAY30, cex = 1.1, adj = c(0, 0.5))
    draw_arrow(P2_sphere_x, P2_sphere_y, P2_proj_x, P2_proj_y,
               col = COL_BLUE, lwd = 1.5, lty = 2)
    segments(-1.8, -1.2, -1.8, -1.2 + 1, col = COL_GRAY50, lwd = 1.5)
    text(-1.8, -1.2 - 0.15, "0", col = COL_GRAY50, cex = 0.9)
    text(-1.8, -1.2 + 1 + 0.15, "R", col = COL_GRAY50, cex = 0.9)
    text(-1.4, -1.2 + 0.5, "Scale", col = COL_GRAY50, cex = 0.9, srt = 90)

    dev.off()

    message("Icosahedron face projection diagram saved to:\n  ",
            normalizePath(paste0(output_path, ".svg")), "\n  ",
            normalizePath(paste0(output_path, ".png")))
  }

  invisible(NULL)
}

# Main execution ----------------------------------------------------------

#' Generate all projection geometry diagrams
#' @param output_dir Directory for output files
#' @export
generate_all_projection_diagrams <- function(output_dir = ".") {

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Generating projection geometry diagrams...")
  message("Output directory: ", normalizePath(output_dir))
  message("")

  # Generate Figure 1: Lambert chord distance
  message("Creating Figure 1: Lambert chord distance diagram...")
  create_lambert_chord_diagram(
    output_path = file.path(output_dir, "fig_lambert_chord_distance")
  )
  message("")

  # Generate Figure 2: Icosahedron face projection
  message("Creating Figure 2: Icosahedron face projection diagram...")
  create_icosa_face_projection_diagram(
    output_path = file.path(output_dir, "fig_icosa_face_projection")
  )
  message("")

  message("All diagrams generated successfully!")
  message("Files created:")
  message("  - fig_lambert_chord_distance.svg/.png")
  message("  - fig_icosa_face_projection.svg/.png")
}

# Example usage (commented out - uncomment to run)
# generate_all_projection_diagrams("C:/Users/Gilles Colling/Documents/dev/hexify/theory_rewrite/03_figures/projection_geometry")
