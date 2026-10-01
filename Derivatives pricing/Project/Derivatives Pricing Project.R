# ==============================================================================
# LLSMS2225 Group Project
# ==============================================================================

# ==============================================================================
# FIGURE 1.2 - Payoff of Absolute (Straddle) vs Quadratic (Variance)
# ==============================================================================

# Clean workspace
rm(list = ls())

# Parameters
K = 30           # strike

# range (Zoomed to [22, 38] to visualize the crossing point clearly)
x = seq(22, 38, .01)
nx = length(x)

# Define Payoffs
payoff_quad = (x - K)^2
payoff_abs  = abs(x - K)

# Plotting
# We set ylim to c(0, 10) to prevent the quadratic curve from flattening the V-shape
# We set asp = 1 to ensure the Straddle appears at a 45-degree angle visually
plot(x, payoff_quad, col="blue", type="l", lwd=2, 
     xlab=expression(S[T]), ylab=expression(psi(S[T])), 
     ylim=c(0, 10), 
     asp=1, 
     las=1,
     main=bquote("Comparison: Quadratic vs Absolute (K="~.(K)*")"))

# Add Absolute Payoff (Red)
lines(x, payoff_abs, col="red", lwd=2)

# Add Reference lines (Ground and Strike) similar to reference code
points(x, rep(0,nx), col="black", type="l", lty=2)
abline(v=K, col="black", lty=2)

# Add Legend
legend("top", legend=c("Quadratic", "Absolute"), 
       col=c("blue", "red"), lwd=2, bty="n")

# Add Grid
grid()


# ==============================================================================
# Q2.2: Pricing via Numerical Integration & Visualization
# ==============================================================================
rm(list=ls()) # Clear environment

# Load libraries once at the start
library(ggplot2)
library(gridExtra) 

# 1. SETUP: Define Market Parameters (Global)
# ------------------------------------------------------------------------------
S0    <- 30      # Initial Stock Price
K     <- 30      # Strike Price
T     <- 1.0     # Time to Maturity (1 year)
r     <- 0.04    # Risk-free Rate (4%)
sigma <- 0.15    # Volatility (15%)

# 2. PRICING: Define Density & Integrands
# ------------------------------------------------------------------------------
# Black-Scholes Log-Normal Parameters
mu_log <- log(S0) + (r - 0.5 * sigma^2) * T
sd_log <- sigma * sqrt(T)

# Function for the Absolute Payoff (Straddle)
integrand_abs <- function(x) {
  payoff  <- abs(x - K)
  density <- dlnorm(x, meanlog = mu_log, sdlog = sd_log)
  return(payoff * density)
}

# Function for the Quadratic Payoff (Variance)
integrand_quad <- function(x) {
  payoff  <- (x - K)^2
  density <- dlnorm(x, meanlog = mu_log, sdlog = sd_log)
  return(payoff * density)
}

# 3. CALCULATE PRICE (Integrate & Discount)
# ------------------------------------------------------------------------------
val_abs_int  <- integrate(integrand_abs, lower = 0, upper = Inf)$value
val_quad_int <- integrate(integrand_quad, lower = 0, upper = Inf)$value

# Discount back to present value using exp(-rT)
price_abs_int  <- exp(-r * T) * val_abs_int
price_quad_int <- exp(-r * T) * val_quad_int

print(paste("Absolute Price (Integration):", round(price_abs_int, 4)))
print(paste("Quadratic Price (Integration):", round(price_quad_int, 4)))

# 4. VISUALIZATION: Plotting Logic
# ------------------------------------------------------------------------------
# Helper: Convert Standard Normal 'u' to Stock Price 'ST'
get_ST <- function(u) S0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * u)

# Unified plotting function
plot_hypothetical <- function(payoff_func, fill_col, line_col, title, y_label_col) {
  
  # A. Data Generation
  u_seq <- seq(-4, 4, length.out = 1000)
  df_smooth <- data.frame(u = u_seq, ST = get_ST(u_seq))
  df_smooth$Payoff <- payoff_func(df_smooth$ST)
  df_smooth$Prob   <- dnorm(u_seq)
  
  # B. Rectangles (Coarse Grid)
  n_rects <- 40
  u_breaks <- seq(-4, 4, length.out = n_rects + 1)
  u_centers <- (u_breaks[-1] + u_breaks[-(n_rects+1)]) / 2 # Midpoints
  
  df_rect <- data.frame(u = u_centers, ST = get_ST(u_centers))
  df_rect$Payoff <- payoff_func(df_rect$ST)
  df_rect$Width  <- diff(u_breaks)[1]
  
  # C. Dynamic Scaling (Fits Probability line to Payoff scale)
  scale_factor <- max(df_smooth$Payoff) / max(df_smooth$Prob)
  
  # D. Plotting
  ggplot() +
    # Rectangles
    geom_col(data = df_rect, aes(x = u, y = Payoff), width = df_rect$Width * 0.9,
             fill = fill_col, color = line_col, alpha = 0.6) +
    # Payoff Curve
    geom_line(data = df_smooth, aes(x = u, y = Payoff, color = "Payoff Function ($)"), size = 1.2) +
    # Probability Curve
    geom_line(data = df_smooth, aes(x = u, y = Prob * scale_factor, color = "Probability Density"),
              linetype = "dotted", size = 1.2) +
    # Scales & Theme
    scale_y_continuous(name = "Payoff Value ($)",
                       sec.axis = sec_axis(~ . / scale_factor, name = "Probability Density")) +
    scale_color_manual(values = c("Payoff Function ($)" = line_col, "Probability Density" = "black")) +
    labs(title = title, subtitle = "Note: Rectangles represent Payoff size (not weighted by probability).", 
         x = "Standard Normal Outcome (u)") +
    theme_minimal() +
    theme(legend.position = "top", 
          axis.title.y = element_text(color = y_label_col),
          axis.title.y.right = element_text(color = "black"))

}

# 5. GENERATE & DISPLAY PLOTS
# ------------------------------------------------------------------------------

# Plot 1: Absolute Payoff
p1 <- plot_hypothetical(
  payoff_func = function(S) abs(S - K),
  fill_col = "mistyrose", line_col = "red",
  title = "Straddle Payoff Profile vs. Standard Normal Shock", y_label_col = "red"
)

# Plot 2: Quadratic Payoff
p2 <- plot_hypothetical(
  payoff_func = function(S) (S - K)^2,
  fill_col = "#3182bd", line_col = "darkblue",
  title = "Quadratic Payoff Profile vs. Standard Normal Shock", y_label_col = "darkblue"
)

# Display plots individually
print(p1)
print(p2)