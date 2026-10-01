# --- Monte Carlo Estimation of Pi in a Square (N=2) ---

# 1. Setup Parameters
n_points <- 1000000 # Increase n for better accuracy
r <- 1              # Radius of the inscribed circle (and half the side of the square)
side <- 2           # Side length of the square (from -1 to 1)

# Theoretical Area Ratio: A_Circle / A_Square = pi/4

# --- 2. Generate Random Points (Scatter) ---

# Generate n_points random x-coordinates uniformly between -1 and 1
x <- runif(n_points, min = -r, max = r)

# Generate n_points random y-coordinates uniformly between -1 and 1
y <- runif(n_points, min = -r, max = r)

# --- 3. Check for Containment (The Circle Test) ---

# Calculate the squared distance from the origin (0, 0) for all points
# Distance^2 = x^2 + y^2
distance_sq <- x^2 + y^2

# Determine which points fall inside the circle (where Distance^2 <= r^2, which is 1^2 = 1)
points_inside <- distance_sq <= r^2

# --- 4. Count and Estimate Pi ---

# Count the total number of points inside the circle
n_in <- sum(points_inside)

# The Monte Carlo Estimate for pi is: 4 * (Points Inside / Total Points)
pi_estimate <- 4 * (n_in / n_points)

# --- 5. Output and Visualization (Plotting a small subset for clarity) ---

cat("--- Monte Carlo Pi Estimation Results ---\n")
cat("Total Random Points (n):", n_points, "\n")
cat("Points Inside Circle (n_in):", n_in, "\n")
cat("\nEstimated value of Pi:", pi_estimate, "\n")
cat("True value of Pi (difference):", pi_estimate - pi, "\n")


# ----------------------------------------------------------------------
# Optional: Plotting a subset of the data for visual confirmation
# We plot a small subset (e.g., 2000 points) because plotting 1 million is slow.
# ----------------------------------------------------------------------

n_plot <- 2000
x_plot <- x[1:n_plot]
y_plot <- y[1:n_plot]
inside_plot <- points_inside[1:n_plot]

# Set up the plot area
plot(x_plot, y_plot, 
     col = ifelse(inside_plot, "red", "blue"), # Red for inside, Blue for outside
     pch = 20,                                 # Solid dots
     cex = 0.5,                                # Smaller point size
     main = paste("Monte Carlo Pi Estimation (n =", n_plot, "points)"),
     xlab = "X coordinate",
     ylab = "Y coordinate",
     xlim = c(-1.1, 1.1),
     ylim = c(-1.1, 1.1),
     asp = 1) # Set aspect ratio to 1 for a true square/circle

# Draw the square boundary
segments(-1, -1, -1, 1, col = "black", lwd = 2)
segments(-1, 1, 1, 1, col = "black", lwd = 2)
segments(1, 1, 1, -1, col = "black", lwd = 2)
segments(1, -1, -1, -1, col = "black", lwd = 2)

# Draw the inscribed circle (using R's built-in plotting capability for circles)
# The function draw.circle() is not base R, so we use a sequence of points for the circle.
theta <- seq(0, 2 * pi, length.out = 100)
lines(r * cos(theta), r * sin(theta), col = "darkgreen", lwd = 2)

legend("topright", 
       legend = c("Points Inside Circle", "Points Outside Circle"),
       col = c("red", "blue"),
       pch = 20)