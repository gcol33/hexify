# Newton-Raphson Iteration Visualization
#
# Demonstrates quadratic convergence of Newton-Raphson method
# for inverse projection azimuth solving.
#
# Figure shows:
# - Function curve f(Az)
# - Initial guess and subsequent iterations
# - Tangent lines showing Newton step geometry
# - Convergence to root

# Color palette
gray30 <- "#4D4D4D"
gray50 <- "#808080"
accent_red <- "#E63946"

# Define a representative function that mimics Snyder residual behavior
# Using f(x) = x^3 - 2x - 2, which has:
# - Single root near x = 1.77
# - Monotonic behavior similar to azimuth residual
# - Smooth derivative throughout domain
f <- function(x) {
  x^3 - 2*x - 2
}

# Derivative for Newton-Raphson
f_prime <- function(x) {
  3*x^2 - 2
}

# Newton-Raphson iteration
newton_step <- function(x) {
  x - f(x) / f_prime(x)
}

# Setup: compute iteration sequence
x0 <- 2.5  # Initial guess (deliberately offset for visibility)
iterations <- 4
x_vals <- numeric(iterations + 1)
x_vals[1] <- x0

for (i in 2:(iterations + 1)) {
  x_vals[i] <- newton_step(x_vals[i - 1])
}

# True root (computed to high precision)
root <- 1.769292354238632

# Setup plot domain
x_min <- 0.5
x_max <- 3.0
x_curve <- seq(x_min, x_max, length.out = 300)
y_curve <- f(x_curve)

y_min <- min(y_curve, -1)
y_max <- max(y_curve, 5)

# Create output directory
dir.create("theory_rewrite/03_figures/newton",
           showWarnings = FALSE, recursive = TRUE)

# Export as SVG
svg("theory_rewrite/03_figures/newton/newton_raphson.svg",
    width = 7, height = 5)

par(mar = c(4, 4, 1, 1))

# Plot function curve
plot(x_curve, y_curve,
     type = "l", lwd = 2, col = gray30,
     xlim = c(x_min, x_max), ylim = c(y_min, y_max),
     xlab = "Az (degrees)", ylab = "f(Az)",
     main = "", las = 1, bty = "n")

# Add y=0 line
abline(h = 0, col = gray50, lwd = 1.5, lty = 1)

# Add vertical grid for readability
grid_x <- seq(ceiling(x_min), floor(x_max), by = 0.5)
abline(v = grid_x, col = gray50, lwd = 0.5, lty = 3)

# Plot Newton-Raphson iterations
for (i in 1:iterations) {
  x_curr <- x_vals[i]
  y_curr <- f(x_curr)
  x_next <- x_vals[i + 1]

  # Draw tangent line from (x_curr, y_curr) to (x_next, 0)
  segments(x_curr, y_curr, x_next, 0,
           col = accent_red, lwd = 1.5, lty = 2)

  # Mark iteration point
  points(x_curr, y_curr, pch = 19, col = accent_red, cex = 1.2)

  # Add iteration label (offset for visibility)
  text(x_curr, y_curr + 0.4,
       labels = bquote(x[.(i-1)]),
       pos = if (i == 1) 3 else 2,
       col = accent_red, cex = 0.9)
}

# Mark final converged point on x-axis
points(x_vals[iterations + 1], 0,
       pch = 4, col = accent_red, cex = 1.5, lwd = 2)
text(x_vals[iterations + 1], -0.3,
     labels = bquote(x[.(iterations)]),
     pos = 1, col = accent_red, cex = 0.9)

# Mark true root
points(root, 0, pch = 4, col = "black", cex = 1.8, lwd = 2.5)
text(root, -0.6, labels = "Root", pos = 1, col = "black", cex = 0.9)

# Add legend
legend("topleft",
       legend = c("f(Az)", "Tangent line", "Iteration point", "Converged"),
       col = c(gray30, accent_red, accent_red, "black"),
       lty = c(1, 2, NA, NA),
       pch = c(NA, NA, 19, 4),
       lwd = c(2, 1.5, NA, 2),
       bty = "n", cex = 0.85)

dev.off()

# Export as PNG (higher resolution for web)
png("theory_rewrite/03_figures/newton/newton_raphson.png",
    width = 7, height = 5, units = "in", res = 300)

par(mar = c(4, 4, 1, 1))

# Plot function curve
plot(x_curve, y_curve,
     type = "l", lwd = 2, col = gray30,
     xlim = c(x_min, x_max), ylim = c(y_min, y_max),
     xlab = "Az (degrees)", ylab = "f(Az)",
     main = "", las = 1, bty = "n")

# Add y=0 line
abline(h = 0, col = gray50, lwd = 1.5, lty = 1)

# Add vertical grid for readability
abline(v = grid_x, col = gray50, lwd = 0.5, lty = 3)

# Plot Newton-Raphson iterations
for (i in 1:iterations) {
  x_curr <- x_vals[i]
  y_curr <- f(x_curr)
  x_next <- x_vals[i + 1]

  # Draw tangent line from (x_curr, y_curr) to (x_next, 0)
  segments(x_curr, y_curr, x_next, 0,
           col = accent_red, lwd = 1.5, lty = 2)

  # Mark iteration point
  points(x_curr, y_curr, pch = 19, col = accent_red, cex = 1.2)

  # Add iteration label (offset for visibility)
  text(x_curr, y_curr + 0.4,
       labels = bquote(x[.(i-1)]),
       pos = if (i == 1) 3 else 2,
       col = accent_red, cex = 0.9)
}

# Mark final converged point on x-axis
points(x_vals[iterations + 1], 0,
       pch = 4, col = accent_red, cex = 1.5, lwd = 2)
text(x_vals[iterations + 1], -0.3,
     labels = bquote(x[.(iterations)]),
     pos = 1, col = accent_red, cex = 0.9)

# Mark true root
points(root, 0, pch = 4, col = "black", cex = 1.8, lwd = 2.5)
text(root, -0.6, labels = "Root", pos = 1, col = "black", cex = 0.9)

# Add legend
legend("topleft",
       legend = c("f(Az)", "Tangent line", "Iteration point", "Converged"),
       col = c(gray30, accent_red, accent_red, "black"),
       lty = c(1, 2, NA, NA),
       pch = c(NA, NA, 19, 4),
       lwd = c(2, 1.5, NA, 2),
       bty = "n", cex = 0.85)

dev.off()

# Print iteration diagnostics
cat("Newton-Raphson Iteration Sequence:\n")
cat("==================================\n\n")
for (i in 1:(iterations + 1)) {
  err <- abs(x_vals[i] - root)
  correct_digits <- if (err > 0) -log10(err) else Inf
  cat(sprintf("Iteration %d: x = %.10f, |error| = %.2e, correct digits = %.1f\n",
              i - 1, x_vals[i], err, correct_digits))
}
cat(sprintf("\nTrue root: %.15f\n", root))
cat("\nQuadratic convergence: each iteration approximately doubles correct digits.\n")
