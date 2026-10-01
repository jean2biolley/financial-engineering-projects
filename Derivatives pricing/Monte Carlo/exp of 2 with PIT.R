# --- Exercise: Draw 1,000 samples from Exp(lambda=2) using the PIT ---

# 1. Define the parameters
n <- 1000       # Number of samples to generate
lambda <- 2     # Parameter of the Exponential distribution

# --- Step 2: Draw n samples from the standard uniform distribution U(0, 1) ---
# The runif() function in R generates uniform random numbers.
u_samples <- runif(n)

# --- Step 3: Apply the Inverse CDF (PIT Transform) to get Exponential samples ---
# The derived Inverse CDF for Exp(lambda) is F_X^{-1}(u) = - (1/lambda) * log(u)
x_samples <- -(1/lambda) * log(u_samples)

# --- Output and Verification (Optional but Recommended) ---

# Print the first few generated samples
cat("First 10 generated Exponential samples (lambda=2):\n")
print(head(x_samples, 10))

# Calculate the mean of the generated samples. 
# Theoretical mean of Exp(lambda=2) is 1/lambda = 1/2 = 0.5.
empirical_mean <- mean(x_samples)
cat("\nEmpirical Mean:", round(empirical_mean, 4), "\n")
cat("Theoretical Mean (1/lambda):", 1/lambda, "\n")

# A histogram can visually confirm the Exponential shape.
hist(x_samples, 
     breaks = 50, 
     main = "Histogram of Generated Exp(2) Samples via PIT", 
     xlab = "x", 
     col = "skyblue", 
     border = "white",
     freq = FALSE)

# Add the theoretical probability density function (PDF) curve for comparison.
curve(dexp(x, rate = lambda), 
      add = TRUE, 
      col = "red", 
      lwd = 2)

legend("topright", 
       legend = c("Generated Samples Density", "Theoretical PDF"), 
       col = c("skyblue", "red"), 
       lty = c(1, 1), 
       lwd = c(10, 2))
