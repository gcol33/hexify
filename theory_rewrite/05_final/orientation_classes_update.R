# Updated orientation-classes chunk for theory.Rmd
# Replace lines 621-651 with this code:

# Panel 1: Class I (Flat-top, 0°)
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class I (Flat-top, 0°)")
h <- hex_v(0, 0, 1.2, rotation = 0)
polygon(h$x, h$y, col = adjustcolor("#457B9D", 0.4), border = "gray30", lwd = 2)
# Highlight top edge in red
segments(h$x[1], h$y[1], h$x[2], h$y[2], col = "#E63946", lwd = 3)
text(0, -1.35, "Aperture 4 (all res)\nAperture 3 (even res)", cex = 0.8)

# Panel 2: Class II (Pointy-top, 30°)
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class II (Pointy-top, 30°)")
h <- hex_v(0, 0, 1.2, rotation = pi/6)
polygon(h$x, h$y, col = adjustcolor("#2A9D8F", 0.4), border = "gray30", lwd = 2)
# Mark top vertex with red dot
points(h$x[1], h$y[1], pch = 19, cex = 1.5, col = "#E63946")
text(0, -1.35, "Aperture 3 (odd res)", cex = 0.8)

# Panel 3: Class III (Skewed, ~19.1°)
plot(NULL, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), asp = 1,
     axes = FALSE, xlab = "", ylab = "", main = "Class III (Skewed, ~19.1°)")
# Dashed parent hexagon (flat-top for reference)
h_parent <- hex_v(0, 0, 1.2, rotation = 0)
polygon(h_parent$x, h_parent$y, col = NA, border = "gray50", lwd = 2, lty = 2)
# Solid child hexagon rotated by arctan(sqrt(3/7))
rot_angle <- atan(sqrt(3/7))
h_child <- hex_v(0, 0, 0.8, rotation = rot_angle)
polygon(h_child$x, h_child$y, col = adjustcolor("#E9C46A", 0.4),
        border = "gray30", lwd = 2)
# Draw arc showing rotation angle
arc_r <- 0.5
arc_theta <- seq(0, rot_angle, length.out = 30)
lines(arc_r * cos(arc_theta), arc_r * sin(arc_theta), col = "#E63946", lwd = 2)
# Label the arc
text(0.6, 0.15, "19.1°", cex = 0.8, col = "#E63946")
text(0, -1.35, "Aperture 7 (all res)\nRotation accumulates", cex = 0.8)

# Changes made:
# 1. Panel 1: Changed lines() to segments() for clearer top edge highlighting
# 2. Panel 2: Added explicit comment "Mark top vertex with red dot"
# 3. Panel 3: Added descriptive comments for each section
# 4. Panel 3: Increased arc resolution from 20 to 30 points for smoother curve
# 5. All panels: Changed comment headers from "Class X:" to "Panel X: Class X"
