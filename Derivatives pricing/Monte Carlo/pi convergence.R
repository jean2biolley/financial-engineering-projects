library(ggplot2)
library(dplyr)

set.seed(40) # For reproducibility

# 1. Setup Parameters
n_max <- 10000
plot_dim <- 10

# 2. Generate all 10,000 points at once
df <- data.frame(
  x = runif(n_max, min = 0, max = plot_dim),
  y = runif(n_max, min = 0, max = plot_dim)
)

# 3. Define Shapes (Same geometry as before)
# Circle: Center (5,5), Radius 2
cx <- 5; cy <- 5; r <- 2
# Square: Center (5,5), Side 2 (Limits 4 to 6)
sx_min <- 4; sx_max <- 6
sy_min <- 4; sy_max <- 6

# 4. Check inclusions
# Returns TRUE (1) or FALSE (0)
is_in_circle <- (df$x - cx)^2 + (df$y - cy)^2 <= r^2
is_in_square <- (df$x >= sx_min & df$x <= sx_max) & 
  (df$y >= sy_min & df$y <= sy_max)

# 5. Simulate the "Loop" using Cumulative Sums
# cumsum() calculates the running total. 
# It effectively shows the count at n=1, n=2, ... n=10000
running_circle_count <- cumsum(is_in_circle)
running_square_count <- cumsum(is_in_square)

# Create a results dataframe
convergence_data <- data.frame(
  n = 1:n_max,
  circle_count = running_circle_count,
  square_count = running_square_count
)

# Avoid division by zero in the very early steps if no points hit the square yet
convergence_data <- convergence_data %>%
  filter(square_count > 0) %>%
  mutate(pi_estimate = circle_count / square_count)

# 6. Plotting the Convergence
ggplot(convergence_data, aes(x = n, y = pi_estimate)) +
  # The estimated line
  geom_line(color = "#2c3e50", size = 0.8) +
  
  # The actual value of Pi (Reference line)
  geom_hline(yintercept = pi, color = "red", linetype = "dashed", size = 1) +
  
  # Annotations
  annotate("text", x = n_max, y = pi + 0.1, label = "True Pi", color = "red", hjust = 1) +
  
  # Formatting
  theme_minimal() +
  labs(
    title = "Convergence of Pi Approximation",
    subtitle = "Ratio of Points (Circle / Square) as N increases from 1 to 10,000",
    x = "Number of Points (N)",
    y = "Estimated Value of Pi"
  ) +
  # Zoom in slightly on the Y-axis to see the stabilization better
  coord_cartesian(ylim = c(2.5, 4.0))
