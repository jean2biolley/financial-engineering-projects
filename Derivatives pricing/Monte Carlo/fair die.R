# --- Exercise: Simulate fair die rolls (n=1000) using two methods ---

# 1. Define Parameters
n <- 1000        # Number of samples
outcomes <- 1:6  # Possible outcomes of the die roll (1, 2, 3, 4, 5, 6)
probs <- 1/6     # Probability of each outcome

# --- Step 1: Draw n samples from the standard uniform distribution U(0, 1) ---
u_samples <- runif(n)

# ====================================================================
## 2. Method A: Generalized Inverse PIT
# ====================================================================

# Define the cumulative probability boundaries (0, 1/6, 2/6, ..., 1)
cdf_breaks <- seq(0, 1, by = probs)

# Use cut() to implement the Generalized Inverse: 
# x_i = k if (k-1)/6 < u_i <= k/6
# 'right = TRUE' ensures the interval is (a, b], matching the discrete CDF logic.
x_pit <- as.numeric(as.character(cut(u_samples, 
                                     breaks = cdf_breaks, 
                                     labels = outcomes,
                                     right = TRUE)))

# ====================================================================
## 3. Method B: Alternative Rounding Formula
# ====================================================================

# Implement the formula: x_i = round(0.5 + 6 * u_i)
# This uses the specific continuous-to-discrete mapping provided in the exercise.
x_round <- round(0.5 + 6 * u_samples)

# ====================================================================
## 4. Verification and Comparison
# ====================================================================

cat("--- Frequencies Comparison (n =", n, ") ---\n")
cat("Theoretical Frequency (1/6):", round(1/6, 4), "\n\n")

# Check Frequencies for PIT Method
pit_frequencies <- prop.table(table(x_pit))
cat("PIT Method Frequencies:\n")
print(round(pit_frequencies, 4))

# Check Frequencies for Rounding Method
round_frequencies <- prop.table(table(x_round))
cat("\nRounding Method Frequencies:\n")
print(round(round_frequencies, 4))


# ====================================================================
## 5. Plotting Histograms
# ====================================================================

# Set up the plotting area to show two graphs side-by-side
par(mfrow = c(1, 2))

# --- Plot 1: Generalized Inverse PIT Histogram ---
hist(x_pit, 
     breaks = seq(0.5, 6.5, by = 1), # Define bins to center bars on integers 1-6
     main = "Method A: Generalized Inverse PIT",
     xlab = "Die Outcome",
     ylab = "Frequency",
     col = "skyblue",
     border = "black",
     xaxt = 'n') # Suppress default x-axis labels

axis(side = 1, at = outcomes) # Set x-axis ticks to the integer outcomes (1-6)

# Add a horizontal line for the expected frequency (n * 1/6)
abline(h = n * probs, col = "red", lty = 2, lwd = 2)
legend("topright", legend = "Expected Count (166.7)", col = "red", lty = 2)

# --- Plot 2: Alternative Rounding Method Histogram ---
hist(x_round, 
     breaks = seq(0.5, 6.5, by = 1), # Define bins to center bars on integers 1-6
     main = "Method B: Alternative Rounding Formula",
     xlab = "Die Outcome",
     ylab = "Frequency",
     col = "lightgreen",
     border = "black",
     xaxt = 'n') # Suppress default x-axis labels

axis(side = 1, at = outcomes) # Set x-axis ticks to the integer outcomes (1-6)

# Add a horizontal line for the expected frequency
abline(h = n * probs, col = "red", lty = 2, lwd = 2)

# Reset plotting layout to default
par(mfrow = c(1, 1))

