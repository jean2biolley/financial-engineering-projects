# Load necessary library for plotting
library(ggplot2)
library(ggforce) # Used for drawing circles easily, if available. 
# If not, the code below handles it manually.

set.seed(40) # Set seed for reproducibility

# 1. Setup Parameters
n_points <- 10000
plot_dim <- 7

# Generate random points in the 10x10 area
df <- data.frame(
  x = runif(n_points, min = 0, max = plot_dim),
  y = runif(n_points, min = 0, max = plot_dim)
)

# 2. Define Shapes
# Circle: Center at (5, 5), Radius 2
cx <- 2; cy <- 3; r <- 2

# Square: Center at (5, 5), Side 2
# This means x ranges from 4 to 6, and y ranges from 4 to 6
sx_min <- 5; sx_max <- 7
sy_min <- 2; sy_max <- 4

# 3. Determine where points landed
# Check if point is in Circle: (x-cx)^2 + (y-cy)^2 <= r^2
df$in_circle <- (df$x - cx)^2 + (df$y - cy)^2 <= r^2

# Check if point is in Square
df$in_square <- (df$x >= sx_min & df$x <= sx_max) & 
  (df$y >= sy_min & df$y <= sy_max)

# Create a category column for plotting colors
df$location <- "Background"
df$location[df$in_circle] <- "Circle"
df$location[df$in_square] <- "Square"
# Note: Since the square is mathematically inside the circle in this setup,
# points in the square are technically in both. 
# We label them 'Square' to make them visible.

# 4. Calculate Pi
count_circle <- sum(df$in_circle)
count_square <- sum(df$in_square)

pi_approx <- count_circle / count_square


#  Plotting
# We create a circle dataframe for drawing the boundary
theta <- seq(0, 2*pi, length.out = 100)
circle_boundary <- data.frame(
  x = cx + r * cos(theta),
  y = cy + r * sin(theta)
)

ggplot(df, aes(x = x, y = y)) +
  # Draw the scatter points
  geom_point(aes(color = location), size = 1.5, alpha = 0.6) +
  
  # Draw the Square Boundary
  geom_rect(xmin = sx_min, xmax = sx_max, ymin = sy_min, ymax = sy_max, 
            fill = NA, color = "blue", size = 1) +
  
  # Draw the Circle Boundary
  geom_path(data = circle_boundary, aes(x = x, y = y), 
            color = "red", size = 1) +
  
  # Set Plot Limits and Labels
  coord_fixed(xlim = c(0, 7), ylim = c(0, 7)) +
  theme_bw() +
  labs(
    title = paste("Monte Carlo Simulation (Pi approx:", round(pi_approx, 4), ")"),
    subtitle = "Red line = Circle (r=2), Blue box = Square (s=2)",
    x = "X Axis",
    y = "Y Axis"
  )


#  Output Results
cat("Points in Circle:", count_circle, "\n")
cat("Points in Square:", count_square, "\n")
cat("Pts in Circle divided by Pts in square:", pi_approx, "\n")





