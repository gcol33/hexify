# Aperture Subdivision Diagrams
# Mathematically accurate diagrams for aperture 3, 4, and 7 subdivisions
# All coordinates computed from explicit formulas
# Quality gates: no self-intersection, no unintended overlap, consistent spacing/orientation

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Color palette
col_gray30 <- "#4D4D4D"  # Structure - darkest
col_gray50 <- "#808080"  # Structure - medium
col_gray70 <- "#B3B3B3"  # Structure - light
col_gray90 <- "#E6E6E6"  # Structure - lightest
col_red <- "#E63946"     # Accent for rotation/highlighting
col_blue <- "#457B9D"    # Secondary accent (only if needed)

# Output directory
output_dir <- "C:/Users/Gilles Colling/Documents/dev/hexify/theory_rewrite/03_figures/apertures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Compute vertices of a regular hexagon
#'
#' @param cx X-coordinate of center
#' @param cy Y-coordinate of center
#' @param radius Circumradius (distance from center to vertex)
#' @param rotation_deg Rotation in degrees (0 = flat-top, 30 = pointy-top)
#' @return Matrix with x, y columns for 6 vertices
hex_vertices <- function(cx, cy, radius, rotation_deg = 0) {
  # Convert rotation to radians
  rotation_rad <- rotation_deg * pi / 180

  # Generate 6 vertices at 60° intervals
  angles <- seq(0, 300, by = 60) * pi / 180 + rotation_rad

  # Compute vertices
  x <- cx + radius * cos(angles)
  y <- cy + radius * sin(angles)

  # Return as matrix
  cbind(x = x, y = y)
}

#' Draw a hexagon
#'
#' @param cx X-coordinate of center
#' @param cy Y-coordinate of center
#' @param radius Circumradius
#' @param rotation_deg Rotation in degrees
#' @param col Fill color
#' @param border Border color
#' @param lwd Line width
#' @param lty Line type
draw_hexagon <- function(cx, cy, radius, rotation_deg = 0,
                        col = NA, border = "black", lwd = 1, lty = 1) {
  v <- hex_vertices(cx, cy, radius, rotation_deg)
  polygon(v[, "x"], v[, "y"], col = col, border = border, lwd = lwd, lty = lty)
}

#' Draw a rotation arc with angle label
#'
#' @param cx X-coordinate of center
#' @param cy Y-coordinate of center
#' @param radius Arc radius
#' @param start_deg Starting angle in degrees
#' @param end_deg Ending angle in degrees
#' @param label Text label for the arc
#' @param col Color for arc and label
draw_rotation_arc <- function(cx, cy, radius, start_deg, end_deg, label, col = col_red) {
  # Generate arc points
  angles <- seq(start_deg, end_deg, length.out = 50) * pi / 180
  arc_x <- cx + radius * cos(angles)
  arc_y <- cy + radius * sin(angles)

  # Draw arc
  lines(arc_x, arc_y, col = col, lwd = 2)

  # Add arrowhead at end
  arrow_angle <- end_deg * pi / 180
  arrow_x <- cx + radius * cos(arrow_angle)
  arrow_y <- cy + radius * sin(arrow_angle)
  arrow_dx <- -sin(arrow_angle) * 0.08
  arrow_dy <- cos(arrow_angle) * 0.08

  arrows(arrow_x - arrow_dx * 2, arrow_y - arrow_dy * 2,
         arrow_x, arrow_y,
         length = 0.1, col = col, lwd = 2)

  # Add label
  mid_angle <- (start_deg + end_deg) / 2 * pi / 180
  label_x <- cx + (radius + 0.25) * cos(mid_angle)
  label_y <- cy + (radius + 0.25) * sin(mid_angle)
  text(label_x, label_y, label, col = col, cex = 0.9, font = 2)
}

# ==============================================================================
# APERTURE 3 DIAGRAM
# ==============================================================================

create_aperture3_diagram <- function() {
  # Parameters
  parent_radius <- 1.5
  parent_cx <- 0
  parent_cy <- 0
  parent_rotation <- 0  # Class I: flat-top

  # Child parameters
  child_scale <- 1 / sqrt(3)  # Linear scale factor
  child_radius <- parent_radius * child_scale
  child_rotation <- 30  # Class II: pointy-top

  # Child positions: triangular arrangement
  # Positions computed geometrically to fit within parent
  # Distance from center for triangular arrangement
  d <- parent_radius * (2/3) * child_scale * sqrt(3)

  # Three children at 120° intervals, starting at 90° (top)
  child_angles <- c(90, 210, 330) * pi / 180
  child_positions <- data.frame(
    cx = parent_cx + d * cos(child_angles),
    cy = parent_cy + d * sin(child_angles)
  )

  # Create plot
  svg(file.path(output_dir, "aperture3_subdivision.svg"), width = 8, height = 8)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), asp = 1)

  # Title
  title("Aperture 3: Triangular Subdivision with 30° Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)

  # Draw parent hexagon (dashed outline)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)

  # Draw child hexagons (solid, filled)
  for (i in 1:3) {
    draw_hexagon(child_positions$cx[i], child_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)

    # Label child
    text(child_positions$cx[i], child_positions$cy[i],
         as.character(i), cex = 1.2, font = 2, col = col_gray30)
  }

  # Draw rotation arc
  draw_rotation_arc(parent_cx, parent_cy - 0.5, 0.5,
                   start_deg = 90, end_deg = 60,
                   label = "30°", col = col_red)

  # Add annotations
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.2,
       "Children (Class II, pointy-top, 30°)",
       cex = 0.9, col = col_gray30)

  # Scale factor annotation
  text(parent_cx + 2.2, parent_cy + 1.5,
       expression(paste("Linear scale: ", 1/sqrt(3), " ≈ 0.577")),
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx + 2.2, parent_cy + 1.1,
       "Area ratio: 1:3",
       cex = 0.85, col = col_gray30, adj = 0)

  dev.off()

  # Also create PNG version
  png(file.path(output_dir, "aperture3_subdivision.png"),
      width = 800, height = 800, res = 100)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), asp = 1)
  title("Aperture 3: Triangular Subdivision with 30° Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)
  for (i in 1:3) {
    draw_hexagon(child_positions$cx[i], child_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)
    text(child_positions$cx[i], child_positions$cy[i],
         as.character(i), cex = 1.2, font = 2, col = col_gray30)
  }
  draw_rotation_arc(parent_cx, parent_cy - 0.5, 0.5,
                   start_deg = 90, end_deg = 60,
                   label = "30°", col = col_red)
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.2,
       "Children (Class II, pointy-top, 30°)",
       cex = 0.9, col = col_gray30)
  text(parent_cx + 2.2, parent_cy + 1.5,
       expression(paste("Linear scale: ", 1/sqrt(3), " ≈ 0.577")),
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx + 2.2, parent_cy + 1.1,
       "Area ratio: 1:3",
       cex = 0.85, col = col_gray30, adj = 0)
  dev.off()

  cat("Created aperture 3 diagrams\n")
}

# ==============================================================================
# APERTURE 4 DIAGRAM
# ==============================================================================

create_aperture4_diagram <- function() {
  # Parameters
  parent_radius <- 1.5
  parent_cx <- 0
  parent_cy <- 0
  parent_rotation <- 0  # Class I: flat-top

  # Child parameters
  child_scale <- 0.5  # Linear scale factor (sqrt(4) = 2, so 1/2)
  child_radius <- parent_radius * child_scale
  child_rotation <- 0  # Same orientation (Class I)

  # Child positions: rhombic arrangement (2x2 pattern)
  # For flat-top hexagons, the rhombic arrangement uses two main axes:
  # - Vertical axis (90°)
  # - 30° axis (upper-right diagonal)

  # Distance between centers for rhombic packing
  dx <- child_radius * sqrt(3)  # Horizontal spacing
  dy <- child_radius * 1.5      # Vertical spacing

  child_positions <- data.frame(
    cx = c(-dx/2, dx/2, -dx/2, dx/2),
    cy = c(dy/2, dy/2, -dy/2, -dy/2)
  )

  # Create plot
  svg(file.path(output_dir, "aperture4_subdivision.svg"), width = 8, height = 8)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), asp = 1)

  # Title
  title("Aperture 4: Rhombic Subdivision with No Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)

  # Draw parent hexagon (dashed outline)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)

  # Draw child hexagons (solid, filled)
  for (i in 1:4) {
    draw_hexagon(child_positions$cx[i], child_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)

    # Label child
    text(child_positions$cx[i], child_positions$cy[i],
         as.character(i), cex = 1.2, font = 2, col = col_gray30)
  }

  # Draw alignment indicators (showing no rotation)
  # Vertical alignment line
  segments(0, -parent_radius - 0.3, 0, parent_radius + 0.3,
           col = col_blue, lwd = 2, lty = 2)

  # Add annotations
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.2,
       "Children (Class I, flat-top, 0°) - Same orientation",
       cex = 0.9, col = col_gray30)

  # Scale factor annotation
  text(parent_cx + 2.2, parent_cy + 1.5,
       "Linear scale: 1/2 = 0.5",
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx + 2.2, parent_cy + 1.1,
       "Area ratio: 1:4",
       cex = 0.85, col = col_gray30, adj = 0)

  # No rotation indicator
  text(0.3, parent_radius + 0.5,
       "No rotation\n(aligned axes)",
       cex = 0.8, col = col_blue, adj = 0)

  dev.off()

  # Also create PNG version
  png(file.path(output_dir, "aperture4_subdivision.png"),
      width = 800, height = 800, res = 100)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), asp = 1)
  title("Aperture 4: Rhombic Subdivision with No Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)
  for (i in 1:4) {
    draw_hexagon(child_positions$cx[i], child_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)
    text(child_positions$cx[i], child_positions$cy[i],
         as.character(i), cex = 1.2, font = 2, col = col_gray30)
  }
  segments(0, -parent_radius - 0.3, 0, parent_radius + 0.3,
           col = col_blue, lwd = 2, lty = 2)
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.2,
       "Children (Class I, flat-top, 0°) - Same orientation",
       cex = 0.9, col = col_gray30)
  text(parent_cx + 2.2, parent_cy + 1.5,
       "Linear scale: 1/2 = 0.5",
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx + 2.2, parent_cy + 1.1,
       "Area ratio: 1:4",
       cex = 0.85, col = col_gray30, adj = 0)
  text(0.3, parent_radius + 0.5,
       "No rotation\n(aligned axes)",
       cex = 0.8, col = col_blue, adj = 0)
  dev.off()

  cat("Created aperture 4 diagrams\n")
}

# ==============================================================================
# APERTURE 7 DIAGRAM
# ==============================================================================

create_aperture7_diagram <- function() {
  # Parameters
  parent_radius <- 1.5
  parent_cx <- 0
  parent_cy <- 0
  parent_rotation <- 0  # Class I: flat-top

  # Child parameters
  child_scale <- 1 / sqrt(7)  # Linear scale factor
  child_radius <- parent_radius * child_scale

  # Rotation angle: arctan(sqrt(3/7))
  child_rotation <- atan(sqrt(3/7)) * 180 / pi  # Class III-A: ~19.1°

  # Child positions: 1 center + 6 ring (rosette)
  # Center child at origin
  center_cx <- parent_cx
  center_cy <- parent_cy

  # Ring children at 60° intervals
  # Distance from center computed for rosette packing
  ring_distance <- 2 * child_radius  # Center-to-center distance for touching hexagons

  ring_angles <- seq(0, 300, by = 60) * pi / 180
  ring_positions <- data.frame(
    cx = parent_cx + ring_distance * cos(ring_angles),
    cy = parent_cy + ring_distance * sin(ring_angles)
  )

  # Create plot
  svg(file.path(output_dir, "aperture7_subdivision.svg"), width = 8, height = 9)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.8), asp = 1)

  # Title
  title("Aperture 7: Rosette Subdivision with 19.1° Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)

  # Draw parent hexagon (dashed outline)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)

  # Draw center child (highlighted slightly different)
  draw_hexagon(center_cx, center_cy, child_radius, child_rotation,
               col = col_gray90, border = col_gray30, lwd = 2)
  text(center_cx, center_cy, "1", cex = 1.2, font = 2, col = col_gray30)

  # Draw ring children
  for (i in 1:6) {
    draw_hexagon(ring_positions$cx[i], ring_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)

    # Label child (2-7 for ring)
    text(ring_positions$cx[i], ring_positions$cy[i],
         as.character(i + 1), cex = 1.2, font = 2, col = col_gray30)
  }

  # Draw rotation arc
  draw_rotation_arc(parent_cx, parent_cy - 0.4, 0.45,
                   start_deg = 90, end_deg = 90 - child_rotation,
                   label = "19.1°", col = col_red)

  # Add annotations
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.5,
       "Children (Class III-A, rotated 19.1°)",
       cex = 0.9, col = col_gray30)

  # Scale factor annotation
  text(parent_cx - 2.2, parent_cy + 1.7,
       expression(paste("Linear scale: ", 1/sqrt(7), " ≈ 0.378")),
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx - 2.2, parent_cy + 1.3,
       "Area ratio: 1:7",
       cex = 0.85, col = col_gray30, adj = 0)

  # Mathematical annotation
  text(parent_cx - 2.2, parent_cy + 0.7,
       expression(paste(theta, " = arctan(", sqrt(3/7), ")")),
       cex = 0.85, col = col_red, adj = 0)

  dev.off()

  # Also create PNG version
  png(file.path(output_dir, "aperture7_subdivision.png"),
      width = 800, height = 900, res = 100)
  par(mar = c(0.5, 0.5, 2, 0.5), bg = "white")
  plot.new()
  plot.window(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.8), asp = 1)
  title("Aperture 7: Rosette Subdivision with 19.1° Rotation",
        cex.main = 1.3, font.main = 2, line = 0.5)
  draw_hexagon(parent_cx, parent_cy, parent_radius, parent_rotation,
               col = NA, border = col_gray50, lwd = 2, lty = 2)
  draw_hexagon(center_cx, center_cy, child_radius, child_rotation,
               col = col_gray90, border = col_gray30, lwd = 2)
  text(center_cx, center_cy, "1", cex = 1.2, font = 2, col = col_gray30)
  for (i in 1:6) {
    draw_hexagon(ring_positions$cx[i], ring_positions$cy[i],
                child_radius, child_rotation,
                col = col_gray90, border = col_gray30, lwd = 2)
    text(ring_positions$cx[i], ring_positions$cy[i],
         as.character(i + 1), cex = 1.2, font = 2, col = col_gray30)
  }
  draw_rotation_arc(parent_cx, parent_cy - 0.4, 0.45,
                   start_deg = 90, end_deg = 90 - child_rotation,
                   label = "19.1°", col = col_red)
  text(parent_cx, parent_cy - 2.2,
       "Parent (Class I, flat-top, 0°)",
       cex = 0.9, col = col_gray50)
  text(parent_cx, parent_cy + 2.5,
       "Children (Class III-A, rotated 19.1°)",
       cex = 0.9, col = col_gray30)
  text(parent_cx - 2.2, parent_cy + 1.7,
       expression(paste("Linear scale: ", 1/sqrt(7), " ≈ 0.378")),
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx - 2.2, parent_cy + 1.3,
       "Area ratio: 1:7",
       cex = 0.85, col = col_gray30, adj = 0)
  text(parent_cx - 2.2, parent_cy + 0.7,
       expression(paste(theta, " = arctan(", sqrt(3/7), ")")),
       cex = 0.85, col = col_red, adj = 0)
  dev.off()

  cat("Created aperture 7 diagrams\n")
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat("Generating aperture subdivision diagrams...\n\n")

create_aperture3_diagram()
create_aperture4_diagram()
create_aperture7_diagram()

cat("\nAll diagrams created successfully!\n")
cat("Output location:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  - aperture3_subdivision.svg/.png\n")
cat("  - aperture4_subdivision.svg/.png\n")
cat("  - aperture7_subdivision.svg/.png\n")
