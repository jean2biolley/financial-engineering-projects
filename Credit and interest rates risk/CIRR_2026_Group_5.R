# =========================================================
# LLSMS2226 Credit and Interest Rates Risk - Group Project
# 
# In this code, we calibrate and price counterparty credit risk (CVA/DVA) 
# on an equity forward contract under wrong-way risk and collateral mitigants:
#
# 1. Model Calibration: Calibrating stock GBM dynamics and fitting Vasicek OU intensity parameters via OLS.
# 2. Yield Curve Bootstrapping: Bootstrapping market survival probabilities and deterministic shift extensions phi(t).
# 3. Baseline Valuation: Deriving analytical and Monte Carlo risk-free fair prices alongside unilateral/bilateral exposures.
# 4. Wrong-Way Risk Simulation: Correlating asset and hazard shocks via Cholesky decomposition to capture first-to-default dynamics.
# 5. CSA Collateral Modeling: Simulating margin call frequencies, minimum transfer amounts (MTA), and uncollateralized exposure.
# 6. Sensitivity Stress Testing: Computing price and XVA Greeks across market parameters, recovery rates, and correlation bumps.
#
# =========================================================
library(ggplot2)
library(patchwork)

# Data import 
df <- read.csv("data.csv", row.names = 1)
dates <- as.Date(rownames(df)) 
df_stock <- data.frame(Date = dates, Price = df$Underlying.stock)
# Plotting function
plot_theme <- theme_bw(base_size = 22, base_family = "serif") +
  theme(plot.title = element_text(face = "bold", color = "#1A2530", size = 24, hjust = 0.5, margin = margin(b = 15)),
        legend.position = "bottom", legend.title = element_blank(), 
        legend.text = element_text(size = 18, color = "#1A2530"), legend.key = element_blank(),
        panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.5), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "#1A2530", fill = NA, linewidth = 0.8),
        axis.title = element_text(face = "bold", color = "#1A2530", size = 20, margin = margin(t = 10, r = 10)),
        axis.text = element_text(color = "#333333", size = 18), plot.margin = margin(15, 15, 10, 10))

# Stock plot
p_stock <- ggplot(df_stock, aes(Date, Price)) + geom_line(color = "#1A2530", linewidth = 1) +
  labs(title = "Historical Evolution of the Underlying Stock", x = "Date", y = "Stock Price") + plot_theme

print(p_stock)

# Export
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/historical_stock_price.pdf", p_stock, device = "pdf", width = 12, height = 6, units = "in")

# ==========================================
# PART 1: DATA ANALYSIS AND CALIBRATION
# ==========================================
dt <- 1 / 52; R <- 0.40

# 1. Calibrate the Stock (GBM)
S <- df$Underlying.stock
log_ret <- diff(log(S))
sigma_S <- sd(log_ret) / sqrt(dt)
mu_S <- (mean(log_ret) / dt) + 0.5 * sigma_S^2
eps_S <- log_ret - mean(log_ret) # Stock shocks
s_B <- df$X5Y.Credit.Spread.Bank          
s_C <- df$X5Y.Credit.Spread.Counterparty

# 2. Calibrate Intensities (Vasicek via OLS)
lambda_B <- df$X5Y.Credit.Spread.Bank / (1 - R)
lambda_C <- df$X5Y.Credit.Spread.Counterparty / (1 - R)

calibrate_ou <- function(lam, dt) {
  mod <- lm(lam[-1] ~ lam[-length(lam)]) # y_t regressed on y_{t-1}
  k <- unname(-log(coef(mod)[2]) / dt)
  list(kappa = k, 
       theta = unname(coef(mod)[1] / (1 - coef(mod)[2])),
       eta = unname(summary(mod)$sigma * sqrt(2 * k / (1 - exp(-2 * k * dt)))),
       residuals = residuals(mod))
}

ou_B <- calibrate_ou(lambda_B, dt)
ou_C <- calibrate_ou(lambda_C, dt)

# 3. Simulation GBM & Plot

n_paths <- 1000; n_steps <- length(S); dt <- 1/52
drift_w <- (mu_S - 0.5 * sigma_S^2) * dt
vol_w <- sigma_S * sqrt(dt)

sim_paths <- matrix(0, nrow = n_steps, ncol = n_paths)
sim_paths[1, ] <- S[1]

# GBM
set.seed(5)
for (i in 2:n_steps) sim_paths[i, ] <- sim_paths[i-1, ] * exp(drift_w + vol_w * rnorm(n_paths))

quants <- t(apply(sim_paths, 1, quantile, probs = c(0.05, 0.50, 0.95)))

df_plot <- data.frame(Date = dates, Real = S, P5 = quants[,1], P50 = quants[,2], P95 = quants[,3])

# Plot
p_ice_tunnel <- ggplot(df_plot, aes(x = Date)) +
  geom_ribbon(aes(ymin = P5, ymax = P95, fill = "90% Confidence Interval"), alpha = 0.2) +
  geom_line(aes(y = P50, color = "Theoretical GBM Median"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = Real, color = "Historical Stock Price"), linewidth = 0.9) +
  scale_color_manual(name = "", values = c("Historical Stock Price" = "#1A2530", "Theoretical GBM Median" = "#005BBB")) + 
  scale_fill_manual(name = "", values = c("90% Confidence Interval" = "#7CA8D8")) + 
  labs(title = "GBM Calibration: Historical Evolution vs. Confidence Tunnel", x = "Date", y = "Stock Price") + 
  plot_theme

print(p_ice_tunnel)

# Export
ggsave("plots/gbm_ice_tunnel.pdf", p_ice_tunnel, device = "pdf", width = 12, height = 6, units = "in")

# 3. Vasicek plots

simulate_vasicek_tunnel <- function(lambda_actual, kappa, theta, eta, dt, n_paths) {
  n_steps <- length(lambda_actual)
  sim_paths <- matrix(0, nrow = n_steps, ncol = n_paths)
  sim_paths[1, ] <- lambda_actual[1]
  
  ar_coeff <- exp(-kappa * dt)
  drift_term <- theta * (1 - ar_coeff)
  vol_term <- eta * sqrt((1 - exp(-2 * kappa * dt)) / (2 * kappa))
  
  set.seed(5)
  for (i in 2:n_steps) {
    sim_paths[i, ] <- sim_paths[i-1, ] * ar_coeff + drift_term + vol_term * rnorm(n_paths)
  }
  
  quants <- t(apply(sim_paths, 1, quantile, probs = c(0.05, 0.50, 0.95)))
  return(quants)
}

# Simulate paths
quants_B <- simulate_vasicek_tunnel(lambda_B, ou_B$kappa, ou_B$theta, ou_B$eta, dt, n_paths)
quants_C <- simulate_vasicek_tunnel(lambda_C, ou_C$kappa, ou_C$theta, ou_C$eta, dt, n_paths)

# Dataframes for plotting
df_vasicek_B <- data.frame(Date = dates, Real = lambda_B, P5 = quants_B[,1], P50 = quants_B[,2], P95 = quants_B[,3])
df_vasicek_C <- data.frame(Date = dates, Real = lambda_C, P5 = quants_C[,1], P50 = quants_C[,2], P95 = quants_C[,3])

# Plot for the Bank
p_bank_tunnel <- ggplot(df_vasicek_B, aes(x = Date)) +
  geom_ribbon(aes(ymin = P5, ymax = P95, fill = "90% Confidence Interval"), alpha = 0.2) +
  geom_line(aes(y = P50, color = "Theoretical Median"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = Real, color = "Historical Intensity"), linewidth = 0.9) +
  scale_color_manual(name = "", values = c("Historical Intensity" = "#1A2530", "Theoretical Median" = "#005BBB")) +
  scale_fill_manual(name = "", values = c("90% Confidence Interval" = "#7CA8D8")) +
  labs(title = "Vasicek Calibration: Bank Default Intensity (B)", x = "Date", y = "Intensity") +
  plot_theme +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15))) # Explicitly centers the title

# Plot for the Counterparty
p_cpty_tunnel <- ggplot(df_vasicek_C, aes(x = Date)) +
  geom_ribbon(aes(ymin = P5, ymax = P95, fill = "90% Confidence Interval"), alpha = 0.2) +
  geom_line(aes(y = P50, color = "Theoretical Median"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = Real, color = "Historical Intensity"), linewidth = 0.9) +
  scale_color_manual(name = "", values = c("Historical Intensity" = "#1A2530", "Theoretical Median" = "#005BBB")) +
  scale_fill_manual(name = "", values = c("90% Confidence Interval" = "#7CA8D8")) +
  labs(title = "Vasicek Calibration: Counterparty Default Intensity (C)", x = "Date", y = "Intensity") +
  plot_theme +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15))) # Explicitly centers the title

print(p_bank_tunnel)
print(p_cpty_tunnel)

# Export as two separate PDFs (Width 12, Height 6 is the standard single-plot size)
ggsave("plots/vasicek_tunnel_bank.pdf", p_bank_tunnel, device = "pdf", width = 12, height = 6, units = "in")
ggsave("plots/vasicek_tunnel_cpty.pdf", p_cpty_tunnel, device = "pdf", width = 12, height = 6, units = "in")

# 4. Correlation matrix & calibration output

cor_matrix <- cor(cbind(Stock = eps_S, Bank = ou_B$residuals, Cpty = ou_C$residuals))

cat("--- CALIBRATION RESULTS ---\n",
    sprintf("Stock (GBM) : mu = %.4f | sigma = %.4f\n", mu_S, sigma_S),
    sprintf("Bank (OU)   : kappa = %.4f | theta = %.4f | eta = %.4f\n", ou_B$kappa, ou_B$theta, ou_B$eta),
    sprintf("Counterparty (OU)   : kappa = %.4f | theta = %.4f | eta = %.4f\n\n", ou_C$kappa, ou_C$theta, ou_C$eta),
    "--- CORRELATION MATRIX ---\n")
print(round(cor_matrix, 4))


# 5. Bootstrapping the survival probability curves

T_mat <- c(1, 3, 5, 7, 10)
r <- 0.02
LGD <- 1 - R
spreads_B <- c(15, 20, 23, 26, 30) / 10000
spreads_C <- c(55, 72, 80, 90, 105) / 10000

# CDS NPV Function
cds_npv <- function(lam_guess, m, known_lams, spread) {
  lams <- c(known_lams, lam_guess)
  PL <- 0
  DL <- 0
  G <- 1.0
  t_prev <- 0
  
  for (i in 1:m) {
    dt <- T_mat[i] - t_prev
    lam <- lams[i]
    
    chunk_val <- exp(-r * t_prev) * G * (1 - exp(-(r + lam) * dt)) / (r + lam)
    
    PL <- PL + spread * chunk_val
    DL <- DL + LGD * lam * chunk_val
    
    G <- G * exp(-lam * dt)
    t_prev <- T_mat[i]
  }
  return(DL - PL)
}

# Bootstrapping loop
bootstrap_lambdas <- function(spreads) {
  lams <- numeric(length(T_mat))
  for (m in seq_along(T_mat)) {
    known <- if (m == 1) numeric(0) else lams[1:(m-1)]
    lams[m] <- uniroot(function(x) cds_npv(x, m, known, spreads[m]), c(0, 1))$root
  }
  return(lams)
}

lambda_B <- bootstrap_lambdas(spreads_B)
lambda_C <- bootstrap_lambdas(spreads_C)

# Market Survival Probabilities G(0,T)
dt_vec <- c(T_mat[1], diff(T_mat))
G_B_mkt <- exp(-cumsum(lambda_B * dt_vec))
G_C_mkt <- exp(-cumsum(lambda_C * dt_vec))

# Output results
cat("--- BOOTSTRAPPED DEFAULT INTENSITIES & SURVIVAL PROBS ---\n",
    sprintf("Bank Lambdas : %s\n", paste(round(lambda_B, 6), collapse=" ")),
    sprintf("Cpty Lambdas : %s\n\n", paste(round(lambda_C, 6), collapse=" ")),
    sprintf("Bank G(0,T)  : %s\n", paste(round(G_B_mkt, 6), collapse=" ")),
    sprintf("Cpty G(0,T)  : %s\n\n", paste(round(G_C_mkt, 6), collapse=" ")))

# Plots of bootstraps

t_seq <- seq(0, 10, by = 0.05)

# Helper function to compute G(0,t) for piecewise constant default intensities
get_G_seq <- function(t_seq, lambdas, T_mat) {
  sapply(t_seq, function(t) {
    if (t <= 0) return(1)
    idx <- which(T_mat >= t)[1]
    if (is.na(idx)) idx <- length(T_mat)
    
    prev_int <- if (idx > 1) sum(lambdas[1:(idx-1)] * diff(c(0, T_mat[1:(idx-1)]))) else 0
    curr_part <- lambdas[idx] * (t - ifelse(idx == 1, 0, T_mat[idx-1]))
    return(exp(-(prev_int + curr_part)))
  })
}

# Helper function to get piecewise lambda(t)
get_lam_t <- function(t_seq, lambdas, T_mat) {
  sapply(t_seq, function(t) {
    idx <- which(T_mat >= t)[1]
    if (is.na(idx)) return(tail(lambdas, 1))
    return(lambdas[idx])
  })
}

# Calculate values for the time sequence
df_boot <- data.frame(
  Time = t_seq,
  G_Bank = get_G_seq(t_seq, lambda_B, T_mat),
  G_Cpty = get_G_seq(t_seq, lambda_C, T_mat),
  Lam_Bank = get_lam_t(t_seq, lambda_B, T_mat),
  Lam_Cpty = get_lam_t(t_seq, lambda_C, T_mat)
)

colors <- c("Bank" = "#005BBB", "Counterparty" = "#8B0000")

# PLOT 1: Survival Probabilities
p_surv <- ggplot(df_boot, aes(x = Time)) +
  geom_line(aes(y = G_Bank, color = "Bank"), linewidth = 1.2) +
  geom_line(aes(y = G_Cpty, color = "Counterparty"), linewidth = 1.2) +
  scale_color_manual(name = "", values = colors) +
  labs(title = expression(bold("Market Survival Probabilities ") ~ bold(G(0,t))),
       x = "Time (Years)",
       y = expression(bold("Survival Probability ") ~ bold(G(0,t)))) +
  plot_theme

# PLOT 2: Default Intensities (Step-plot)
p_haz <- ggplot(df_boot, aes(x = Time)) +
  geom_step(aes(y = Lam_Bank, color = "Bank"), linewidth = 1.2) +
  geom_step(aes(y = Lam_Cpty, color = "Counterparty"), linewidth = 1.2) +
  scale_color_manual(name = "", values = colors) +
  labs(title = expression(bold("Bootstrapped Default Intensities ") ~ bold(lambda(t))),
       x = "Time (Years)",
       y = expression(bold("Default Intensity ") ~ bold(lambda(t)))) +
  plot_theme

# 3. Combine and export
combined_plot <- p_surv + p_haz + plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(combined_plot)
ggsave("plots/bootstrap_results.pdf", combined_plot, device = "pdf", width = 16, height = 8, units = "in")


# 6. Determine the shift functions phi(t)

# Function to get the market forward rate f_mkt(0,t) 
get_f_mkt <- function(t, lambdas, T_mat) {
  if (t <= 0) return(lambdas[1])
  idx <- which(T_mat >= t)[1]
  return(if (is.na(idx)) tail(lambdas, 1) else lambdas[idx])
}

# Function to calculate the Vasicek base forward rate f_y(0,t)
get_f_y <- function(t, kappa, theta, eta) {
  term1 <- theta * (1 - exp(-kappa * t))
  term2 <- (eta^2 / (2 * kappa^2)) * (1 - exp(-kappa * t))^2
  return(term1 - term2)
}

# The final Shift Function phi(t)
phi_t <- function(t, lambdas, T_mat, kappa, theta, eta) {
  return(get_f_mkt(t, lambdas, T_mat) - get_f_y(t, kappa, theta, eta))
}


# 7. Plots

# historical intensities 
lambda_B_hist <- s_B / LGD
lambda_C_hist <- s_C / LGD

colors <- c("Bank" = "#005BBB", "Counterparty" = "#8B0000")

# --- PLOT 1: Historical Default Intensities ---
df_hist <- data.frame(Weeks = 1:length(lambda_B_hist), Bank = lambda_B_hist, Cpty = lambda_C_hist)

p_hist <- ggplot(df_hist, aes(x = Weeks)) +
  geom_line(aes(y = Bank, color = "Bank"), linewidth = 1.2) +
  geom_line(aes(y = Cpty, color = "Counterparty"), linewidth = 1.2) +
  scale_color_manual(name = "", values = colors) +
  labs(title = "Historical Default Intensities", 
       x = "Time (Weeks)", 
       y = expression(bold("Intensity (") * bold(lambda) * bold(")"))) + 
  plot_theme

print(p_hist)
ggsave("plots/1_historical_intensities.pdf", p_hist, device = "pdf", width = 12, height = 7, units = "in")

# --- PLOT 2: Correlation Scatter Plot ---
df_scatter <- data.frame(Stock_Res = eps_S, Cpty_Res = ou_C$residuals)

p_scatter <- ggplot(df_scatter, aes(x = Stock_Res, y = Cpty_Res)) +
  geom_point(color = "#8B0000", alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", color = "#1A2530", linetype = "dashed", linewidth = 1.2, se = FALSE) +
  labs(title = "Stock Returns vs. Counterparty Intensity Shocks", x = "Stock Returns (Residuals)", y = "Counterparty Intensity (Residuals)") +
  plot_theme

print(p_scatter)
ggsave("plots/2_correlation_scatter.pdf", p_scatter, device = "pdf", width = 12, height = 7, units = "in")


# --- PLOT 3: The Shift Function phi(t) over 10 years ---
t_seq <- seq(0, 10, by = 0.1)
df_phi <- data.frame(
  Time = t_seq,
  Bank = sapply(t_seq, function(t) phi_t(t, lambda_B, T_mat, ou_B$kappa, ou_B$theta, ou_B$eta)),
  Cpty = sapply(t_seq, function(t) phi_t(t, lambda_C, T_mat, ou_C$kappa, ou_C$theta, ou_C$eta))
)

p_phi_10y <- ggplot(df_phi, aes(x = Time)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 1) +
  geom_line(aes(y = Bank, color = "Bank"), linewidth = 1.2) +
  geom_line(aes(y = Cpty, color = "Counterparty"), linewidth = 1.2) +
  scale_color_manual(name = "", values = colors) +
  labs(title = expression(bold("Deterministic Shift Function ") ~ bold(phi(t))), 
       x = "Time (Years)", 
       y = expression(bold("Shift Value ") ~ bold(phi(t)))) +
  plot_theme

print(p_phi_10y)
ggsave("plots/4_deterministic_shift_10y.pdf", p_phi_10y, device = "pdf", width = 12, height = 7, units = "in")

# --- PLOT 4: RESIDUAL DIAGNOSTICS (QQ PLOTS) ---
df_qq <- data.frame(
  Residuals = c(eps_S, ou_B$residuals, ou_C$residuals),
  Model = factor(rep(c("Stock (GBM)", "Bank (Vasicek)", "Counterparty (Vasicek)"), 
                     times = c(length(eps_S), length(ou_B$residuals), length(ou_C$residuals))))
)

df_qq$Model <- factor(df_qq$Model, levels = c("Stock (GBM)", "Bank (Vasicek)", "Counterparty (Vasicek)"))

p_qq <- ggplot(df_qq, aes(sample = Residuals)) +
  stat_qq(shape = 21, color = "#005BBB", fill = "#7CA8D8", size = 3, alpha = 0.7) +
  stat_qq_line(color = "#8B0000", linewidth = 1.2, linetype = "dashed") +
  facet_wrap(~ Model, scales = "free") +
  labs(title = "Residual Diagnostics: Normal Q-Q Plots", 
       x = "Theoretical Quantiles", 
       y = "Sample Quantiles") +
  plot_theme +
  theme(
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5, margin = margin(b = 25)),
    axis.title.x = element_text(face = "bold", size = 22, margin = margin(t = 15)),
    axis.title.y = element_text(face = "bold", size = 22, margin = margin(r = 15)),
    axis.text = element_text(size = 18, color = "#333333"),
    strip.background = element_rect(fill = "#F4F6F8", color = "#1A2530", linewidth = 1),
    strip.text = element_text(face = "bold", size = 20, margin = margin(t = 12, b = 12)),
    panel.border = element_rect(color = "#1A2530", fill = NA, linewidth = 1),
    plot.margin = margin(20, 20, 20, 20)
  )

print(p_qq)

ggsave("plots/3_residual_qq_plots.pdf", p_qq, device = "pdf", width = 16, height = 7, units = "in")

# ==========================================================
# PART 2 : FAIR PRICE, EPE/ENE, CVA & DVA
# ==========================================================

# 1. Forward parameters
S0 <- tail(S, 1)
K <- 150
Time_Maturity <- 5
n_steps <- 260
dt_sim <- Time_Maturity / n_steps

# 2. Fair Price Analytic and Monte Carlo

V0_analytic <- S0 - K * exp(-r * Time_Maturity)

drift_Q <- (r - 0.5 * sigma_S^2) * dt_sim
vol_Q <- sigma_S * sqrt(dt_sim)

sim_S <- matrix(0, nrow = n_steps + 1, ncol = n_paths)
sim_S[1, ] <- S0

set.seed(42)
for (t in 2:(n_steps + 1)) {
  Z_half <- rnorm(n_paths / 2)
  Z <- c(Z_half, -Z_half)
  sim_S[t, ] <- sim_S[t - 1, ] * exp(drift_Q + vol_Q * Z)
}

V0_MC <- mean(sim_S[n_steps + 1, ] - K) * exp(-r * Time_Maturity)

# 3. EPE and ENE Profiles

EPE <- numeric(n_steps + 1)
ENE <- numeric(n_steps + 1)
time_grid <- seq(0, Time_Maturity, length.out = n_steps + 1) 

for (t in 1:(n_steps + 1)) {
  tau <- Time_Maturity - time_grid[t]
  Vt <- sim_S[t, ] - K * exp(-r * tau)
  
  EPE[t] <- mean(pmax(Vt, 0))
  ENE[t] <- mean(pmin(Vt, 0))
}

# --- PLOT 5: EPE AND ENE EXPOSURE PROFILES ---
df_expo <- data.frame(Time = time_grid, EPE = EPE, ENE = ENE)
colors_expo <- c("Expected Positive Exposure (EPE)" = "#005BBB", "Expected Negative Exposure (ENE)" = "#8B0000")

p_expo <- ggplot(df_expo, aes(x = Time)) +
  geom_hline(yintercept = 0, color = "#1A2530", linewidth = 1) +
  geom_ribbon(aes(ymin = 0, ymax = EPE, fill = "Expected Positive Exposure (EPE)"), alpha = 0.1) +
  geom_ribbon(aes(ymin = ENE, ymax = 0, fill = "Expected Negative Exposure (ENE)"), alpha = 0.1) +
  geom_line(aes(y = EPE, color = "Expected Positive Exposure (EPE)"), linewidth = 1.2) +
  geom_line(aes(y = ENE, color = "Expected Negative Exposure (ENE)"), linewidth = 1.2) +
  scale_color_manual(name = "", values = colors_expo) +
  scale_fill_manual(name = "", values = colors_expo, guide = "none") + 
  labs(title = "Forward Contract Exposure Profiles (EPE & ENE)", x = "Time (Years)", y = "Expected Exposure (MtM)") +
  plot_theme

print(p_expo)
ggsave("plots/5_exposure_profiles.pdf", p_expo, device = "pdf", width = 12, height = 7, units = "in")

# 4. Interpolate Market Survival Probabilities G(0,t)

get_G_vec <- function(t_vec, lambdas, T_mat) {
  sapply(t_vec, function(t) {
    if (t <= 0) return(1)
    idx <- which(T_mat >= t)[1]
    if (is.na(idx)) idx <- length(T_mat)
    prev_int <- if(idx > 1) sum(lambdas[1:(idx-1)] * diff(c(0, T_mat[1:(idx-1)]))) else 0
    curr_part <- lambdas[idx] * (t - ifelse(idx==1, 0, T_mat[idx-1]))
    exp(-(prev_int + curr_part))
  })
}

G_B_grid <- get_G_vec(time_grid, lambda_B, T_mat) 
G_C_grid <- get_G_vec(time_grid, lambda_C, T_mat) 

PD_B <- -diff(G_B_grid)
PD_C <- -diff(G_C_grid)

# 5. CVA & DVA by Monte Carlo

discount_factors <- exp(-r * time_grid[-1])
EPE_disc_MC <- EPE[-1] * discount_factors
ENE_disc_MC <- abs(ENE[-1]) * discount_factors

CVA_MC <- LGD * sum(EPE_disc_MC * PD_C)
DVA_MC <- LGD * sum(ENE_disc_MC * PD_B)


# 6. CVA & DVA Analytic (Black-Scholes)

d1 <- function(S, K, r, sigma, t) { (log(S/K) + (r + 0.5 * sigma^2) * t) / (sigma * sqrt(t)) }
d2 <- function(S, K, r, sigma, t) { d1(S, K, r, sigma, t) - sigma * sqrt(t) }

BS_Call <- function(S, K, r, sigma, t) {
  if (t == 0) return(max(S - K, 0))
  S * pnorm(d1(S, K, r, sigma, t)) - K * exp(-r * t) * pnorm(d2(S, K, r, sigma, t))
}
BS_Put <- function(S, K, r, sigma, t) {
  if (t == 0) return(max(K - S, 0))
  K * exp(-r * t) * pnorm(-d2(S, K, r, sigma, t)) - S * pnorm(-d1(S, K, r, sigma, t))
}

EPE_disc_analytic <- numeric(n_steps)
ENE_disc_analytic <- numeric(n_steps)

for (i in 1:n_steps) {
  t_step <- time_grid[i+1]
  K_opt <- K * exp(-r * (Time_Maturity - t_step))
  
  EPE_disc_analytic[i] <- BS_Call(S0, K_opt, r, sigma_S, t_step)
  ENE_disc_analytic[i] <- BS_Put(S0, K_opt, r, sigma_S, t_step)
}

CVA_analytic <- LGD * sum(EPE_disc_analytic * PD_C)
DVA_analytic <- LGD * sum(ENE_disc_analytic * PD_B)

# 7. First-to-Default Adjustment

# Multiply default probability by the survival probability of the other entity
G_B_prev <- head(G_B_grid, -1) 
G_C_prev <- head(G_C_grid, -1)

PD_C_FtD <- PD_C * G_B_prev
PD_B_FtD <- PD_B * G_C_prev

CVA_analytic_FtD <- LGD * sum(EPE_disc_analytic * PD_C_FtD)
DVA_analytic_FtD <- LGD * sum(ENE_disc_analytic * PD_B_FtD)

CVA_MC_FtD <- LGD * sum(EPE_disc_MC * PD_C_FtD)
DVA_MC_FtD <- LGD * sum(ENE_disc_MC * PD_B_FtD)


# 8. Print All Results 

Adjusted_Price_Analytic <- V0_analytic - CVA_analytic + DVA_analytic
Adjusted_Price_MC <- V0_MC - CVA_MC + DVA_MC
Adjusted_Price_FtD <- V0_analytic - CVA_analytic_FtD + DVA_analytic_FtD

cat("--- 1. FAIR PRICE (RISK-FREE) ---\n")
cat("Analytical Price :", round(V0_analytic, 4), "\n")
cat("Monte Carlo Price:", round(V0_MC, 4), "\n\n")

cat("--- 2. CVA & DVA (UNDER INDEPENDENCE) ---\n")
cat(sprintf("Analytic     | CVA: %.4f | DVA: %.4f | Adj Price: %.4f\n", CVA_analytic, DVA_analytic, Adjusted_Price_Analytic))
cat(sprintf("Monte Carlo  | CVA: %.4f | DVA: %.4f | Adj Price: %.4f\n\n", CVA_MC, DVA_MC, Adjusted_Price_MC))

cat("--- 3. CVA & DVA : FIRST-TO-DEFAULT (UNDER INDEPENDENCE) ---\n")
cat(sprintf("Analytic FtD | CVA: %.4f | DVA: %.4f | Adj Price: %.4f\n", CVA_analytic_FtD, DVA_analytic_FtD, Adjusted_Price_FtD))
cat(sprintf("MC FtD       | CVA: %.4f | DVA: %.4f | Adj Price: %.4f\n\n", CVA_MC_FtD, DVA_MC_FtD, V0_MC - CVA_MC_FtD + DVA_MC_FtD))




# ==========================================================
# PART 3.1: WRONG-WAY RISK (CORRELATION & FIRST-TO-DEFAULT)
# ==========================================================

# 1. Setup & pre-computation
L_chol <- t(chol(cor_matrix)) 
n_paths_corr <- 20000  

sim_S_corr <- matrix(0, nrow = n_steps + 1, ncol = n_paths_corr)
sim_yB     <- matrix(0, nrow = n_steps + 1, ncol = n_paths_corr)
sim_yC     <- matrix(0, nrow = n_steps + 1, ncol = n_paths_corr)

sim_S_corr[1, ] <- S0
sim_yB[1, ]     <- 0   
sim_yC[1, ]     <- 0

phi_B_vec <- sapply(time_grid, function(t) phi_t(t, lambda_B, T_mat, ou_B$kappa, ou_B$theta, ou_B$eta))
phi_C_vec <- sapply(time_grid, function(t) phi_t(t, lambda_C, T_mat, ou_C$kappa, ou_C$theta, ou_C$eta))

# 2. Simulation loop
set.seed(42)
for (t in 2:(n_steps + 1)) {
  
  Z_indep <- matrix(rnorm(3 * n_paths_corr), nrow = 3, ncol = n_paths_corr)
  Z_corr <- L_chol %*% Z_indep 
  
  sim_S_corr[t, ] <- sim_S_corr[t-1, ] * exp(drift_Q + vol_Q * Z_corr[1, ])
  sim_yB[t, ] <- sim_yB[t-1, ] + ou_B$kappa * (ou_B$theta - sim_yB[t-1, ]) * dt + ou_B$eta * sqrt(dt) * Z_corr[2, ]
  sim_yC[t, ] <- sim_yC[t-1, ] + ou_C$kappa * (ou_C$theta - sim_yC[t-1, ]) * dt + ou_C$eta * sqrt(dt) * Z_corr[3, ]
}

# 3. Full default intensities (floored at zero) & integrated hazards
sim_lamB <- pmax(sweep(sim_yB, 1, phi_B_vec, "+"), 0)
sim_lamC <- pmax(sweep(sim_yC, 1, phi_C_vec, "+"), 0)

int_lamB <- apply(sim_lamB[1:n_steps, ], 2, cumsum) * dt
int_lamC <- apply(sim_lamC[1:n_steps, ], 2, cumsum) * dt

# 4. Default time simulation & Pathwise XVA
E_B <- rexp(n_paths_corr)
E_C <- rexp(n_paths_corr)

CVA_paths <- numeric(n_paths_corr)
DVA_paths <- numeric(n_paths_corr)

for (p in 1:n_paths_corr) {
  
  idx_B <- which(int_lamB[, p] > E_B[p])[1]
  idx_C <- which(int_lamC[, p] > E_C[p])[1]
  
  tau_B <- ifelse(is.na(idx_B), Inf, time_grid[idx_B + 1])
  tau_C <- ifelse(is.na(idx_C), Inf, time_grid[idx_C + 1])
  tau_min <- min(tau_B, tau_C)
  
  if (tau_min <= Time_Maturity) {
    
    # Counterparty defaults first or simultaneously (CVA)
    if (tau_C <= tau_B) {
      Vt <- sim_S_corr[idx_C + 1, p] - K * exp(-r * (Time_Maturity - tau_C))
      if (Vt > 0) CVA_paths[p] <- LGD * Vt * exp(-r * tau_C)
    }
    
    # Bank defaults first or simultaneously (DVA)
    if (tau_B <= tau_C) {
      Vt <- sim_S_corr[idx_B + 1, p] - K * exp(-r * (Time_Maturity - tau_B))
      if (Vt < 0) DVA_paths[p] <- LGD * abs(Vt) * exp(-r * tau_B)
    }
  }
}

CVA_Final_WWR <- mean(CVA_paths)
DVA_Final_WWR <- mean(DVA_paths)
Adjusted_Price_WWR <- V0_analytic - CVA_Final_WWR + DVA_Final_WWR

# 5. Output
cat("\n--- PART 3.1 RESULTS: CORRELATION AND WRONG-WAY RISK ---\n",
    sprintf("CVA (Independent) : %.4f  -->  CVA (Correlated WWR) : %.4f\n", CVA_analytic_FtD, CVA_Final_WWR),
    sprintf("DVA (Independent) : %.4f  -->  DVA (Correlated WWR) : %.4f\n\n", DVA_analytic_FtD, DVA_Final_WWR),
    sprintf("Adjusted Price (Independent FtD) : %.4f\n", Adjusted_Price_FtD),
    sprintf("Adjusted Price (Correlated WWR)  : %.4f\n", Adjusted_Price_WWR))

# ==========================================================
# PART 3.2: COLLATERAL AND CREDIT SUPPORT ANNEX (CSA)
# ==========================================================

# 1. Function to compute CVA and DVA with Collateral
# Unilateral CSA: Counterparty posts collateral for bank's positive exposure.
# Bank does not post collateral.
compute_collateral_XVA <- function(freq_margin_input = 2, MTA_input = 5) {
  
  sim_C_tmp <- matrix(0, nrow = n_steps + 1, ncol = n_paths_corr)
  
  # 1.1 Build collateral account
  for (t in 2:(n_steps + 1)) {
    Vt <- sim_S_corr[t, ] - K * exp(-r * (Time_Maturity - time_grid[t]))
    C_target <- pmax(Vt, 0) # Target collateral covers positive exposure
    
    # Base case: carry over previous collateral for all paths
    sim_C_tmp[t, ] <- sim_C_tmp[t - 1, ]
    
    # Margin call dates
    if ((t - 1) %% freq_margin_input == 0) {
      margin_call <- C_target - sim_C_tmp[t - 1, ]
      update_mask <- abs(margin_call) > MTA_input
      
      # Overwrite collateral ONLY for paths exceeding the MTA
      sim_C_tmp[t, update_mask] <- C_target[update_mask]
    }
  }
  
  # 1.2 Recompute CVA and DVA with collateral
  CVA_paths_tmp <- numeric(n_paths_corr)
  DVA_paths_tmp <- numeric(n_paths_corr)
  
  for (p in 1:n_paths_corr) {
    idx_B <- which(int_lamB[, p] > E_B[p])[1]
    idx_C <- which(int_lamC[, p] > E_C[p])[1]
    
    tau_B <- ifelse(is.na(idx_B), Inf, time_grid[idx_B + 1])
    tau_C <- ifelse(is.na(idx_C), Inf, time_grid[idx_C + 1])
    tau_min <- min(tau_B, tau_C)
    
    if (tau_min <= Time_Maturity) {
      
      # Counterparty defaults first or simultaneously (CVA)
      if (tau_C <= tau_B) {
        Vt <- sim_S_corr[idx_C + 1, p] - K * exp(-r * (Time_Maturity - tau_C))
        exposure_collat <- max(Vt - sim_C_tmp[idx_C + 1, p], 0)
        
        if (exposure_collat > 0) CVA_paths_tmp[p] <- LGD * exposure_collat * exp(-r * tau_C)
      }
      
      # Bank defaults first or simultaneously (DVA - no collateral posted by bank)
      if (tau_B <= tau_C) {
        Vt <- sim_S_corr[idx_B + 1, p] - K * exp(-r * (Time_Maturity - tau_B))
        if (Vt < 0) DVA_paths_tmp[p] <- LGD * abs(Vt) * exp(-r * tau_B)
      }
    }
  }
  
  CVA_tmp <- mean(CVA_paths_tmp)
  DVA_tmp <- mean(DVA_paths_tmp)
  
  return(c(CVA = CVA_tmp, DVA = DVA_tmp, Adjusted_Price = V0_analytic - CVA_tmp + DVA_tmp))
}


# 2. Main CSA case

res_collat_main <- compute_collateral_XVA(freq_margin_input = 2, MTA_input = 5)

cat("\n--- PART 3.2 RESULTS: COLLATERAL EFFECT ---\n",
    sprintf("CVA without collateral, with WWR : %.4f\n", CVA_Final_WWR),
    sprintf("DVA without collateral, with WWR : %.4f\n", DVA_Final_WWR),
    sprintf("Adjusted price without collat    : %.4f\n\n", Adjusted_Price_WWR),
    sprintf("CVA with collateral, MTA = 5     : %.4f\n", res_collat_main["CVA"]),
    sprintf("DVA with collateral, MTA = 5     : %.4f\n", res_collat_main["DVA"]),
    sprintf("Adjusted price with collateral   : %.4f\n", res_collat_main["Adjusted_Price"]))


# 3. Impact of collateral frequency and MTA

collateral_sensi <- rbind(
  "Weekly collateral, MTA = 5"    = compute_collateral_XVA(freq_margin_input = 1, MTA_input = 5),
  "Two-week collateral, MTA = 5"  = res_collat_main, # Reusing the calculation from above
  "Monthly collateral, MTA = 5"   = compute_collateral_XVA(freq_margin_input = 4, MTA_input = 5),
  "Two-week collateral, MTA = 0"  = compute_collateral_XVA(freq_margin_input = 2, MTA_input = 0),
  "Two-week collateral, MTA = 10" = compute_collateral_XVA(freq_margin_input = 2, MTA_input = 10)
)

collateral_sensi_df <- as.data.frame(round(collateral_sensi, 4))
print(collateral_sensi_df)

# ==========================================================
# PART 3.3: SENSITIVITY ANALYSIS
# ==========================================================

# 0. Helper Functions
# Restores positive-definiteness to a bumped correlation matrix
make_valid_correlation <- function(C) {
  C <- (C + t(C)) / 2
  diag(C) <- 1
  eig <- eigen(C)
  eig$values <- pmax(eig$values, 1e-8)
  C_new <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  D_inv <- diag(1 / sqrt(diag(C_new)))
  C_new <- D_inv %*% C_new %*% D_inv
  C_new <- (C_new + t(C_new)) / 2
  diag(C_new) <- 1
  return(C_new)
}

# Re-bootstrap default intensities dynamically when 'r' or 'LGD' are bumped
bootstrap_lambdas_sensi <- function(spreads, maturities, r_bumped, LGD_bumped) {
  lams_out <- numeric(length(maturities))
  for (m in seq_along(maturities)) {
    known <- if (m == 1) numeric(0) else lams_out[1:(m - 1)]
    
    # Self-contained NPV root function to avoid global variable conflicts
    root_func <- function(lam_guess) {
      lams <- c(known, lam_guess)
      PL <- 0; DL <- 0; G <- 1.0; t_prev <- 0
      for (i in 1:m) {
        dt_chunk <- maturities[i] - t_prev
        chunk_val <- exp(-r_bumped * t_prev) * G * (1 - exp(-(r_bumped + lams[i]) * dt_chunk)) / (r_bumped + lams[i])
        PL <- PL + spreads[m] * chunk_val
        DL <- DL + LGD_bumped * lams[i] * chunk_val
        G <- G * exp(-lams[i] * dt_chunk)
        t_prev <- maturities[i]
      }
      return(DL - PL)
    }
    lams_out[m] <- uniroot(root_func, c(0, 1))$root
  }
  return(lams_out)
}

# 1. Full Sensitivity Pricer
compute_Full_Sensi <- function(
    K_bump=0, T_bump=0, r_bump=0, mu_S_bump=0, sigma_S_bump=0,
    rho_SB_bump=0, rho_SC_bump=0, rho_BC_bump=0,
    kappa_B_bump=0, theta_B_bump=0, eta_B_bump=0,
    kappa_C_bump=0, theta_C_bump=0, eta_C_bump=0,
    R_B_bump=0, R_C_bump=0, n_paths_sensi=20000) {
  
  # a. Apply Bumps
  K_sim <- K * (1 + K_bump); T_sim <- Time_Maturity + T_bump; r_sim <- r + r_bump
  sigma_sim <- sigma_S * (1 + sigma_S_bump)
  
  LGD_B_sim <- 1 - min(max(R + R_B_bump, 0.0001), 0.9999)
  LGD_C_sim <- 1 - min(max(R + R_C_bump, 0.0001), 0.9999)
  
  kappa_B_sim <- max(ou_B$kappa * (1 + kappa_B_bump), 1e-8)
  theta_B_sim <- ou_B$theta * (1 + theta_B_bump)
  eta_B_sim   <- max(ou_B$eta * (1 + eta_B_bump), 1e-8)
  
  kappa_C_sim <- max(ou_C$kappa * (1 + kappa_C_bump), 1e-8)
  theta_C_sim <- ou_C$theta * (1 + theta_C_bump)
  eta_C_sim   <- max(ou_C$eta * (1 + eta_C_bump), 1e-8)
  
  # b. Re-bootstrap default intensities
  lamB_mkt <- bootstrap_lambdas_sensi(spreads_B, T_mat, r_sim, LGD_B_sim)
  lamC_mkt <- bootstrap_lambdas_sensi(spreads_C, T_mat, r_sim, LGD_C_sim)
  
  # c. Dynamic Time Grid & Shifts
  n_steps_sim <- round(T_sim * 52)
  dt_sim <- T_sim / n_steps_sim
  time_grid_sim <- seq(0, T_sim, length.out = n_steps_sim + 1)
  
  phi_B_vec <- sapply(time_grid_sim, function(t) phi_t(t, lamB_mkt, T_mat, kappa_B_sim, theta_B_sim, eta_B_sim))
  phi_C_vec <- sapply(time_grid_sim, function(t) phi_t(t, lamC_mkt, T_mat, kappa_C_sim, theta_C_sim, eta_C_sim))
  
  # d. Bumped Correlation Matrix
  cor_sim <- cor_matrix
  cor_sim[1, 2] <- cor_sim[2, 1] <- min(max(cor_sim[1, 2] + rho_SB_bump, -0.99), 0.99)
  cor_sim[1, 3] <- cor_sim[3, 1] <- min(max(cor_sim[1, 3] + rho_SC_bump, -0.99), 0.99)
  cor_sim[2, 3] <- cor_sim[3, 2] <- min(max(cor_sim[2, 3] + rho_BC_bump, -0.99), 0.99)
  L_chol_sim <- t(chol(make_valid_correlation(cor_sim)))
  
  # e. Stochastic Engine
  drift_Q_sim <- (r_sim - 0.5 * sigma_sim^2) * dt_sim
  vol_Q_sim <- sigma_sim * sqrt(dt_sim)
  
  S_sim <- matrix(0, nrow = n_steps_sim + 1, ncol = n_paths_sensi)
  yB_sim <- matrix(0, nrow = n_steps_sim + 1, ncol = n_paths_sensi)
  yC_sim <- matrix(0, nrow = n_steps_sim + 1, ncol = n_paths_sensi)
  S_sim[1, ] <- S0
  
  set.seed(42)
  for (t in 2:(n_steps_sim + 1)) {
    Z_corr <- L_chol_sim %*% matrix(rnorm(3 * n_paths_sensi), nrow = 3)
    
    S_sim[t, ] <- S_sim[t - 1, ] * exp(drift_Q_sim + vol_Q_sim * Z_corr[1, ])
    yB_sim[t, ] <- yB_sim[t-1, ] + kappa_B_sim * (theta_B_sim - yB_sim[t-1, ]) * dt_sim + eta_B_sim * sqrt(dt_sim) * Z_corr[2, ]
    yC_sim[t, ] <- yC_sim[t-1, ] + kappa_C_sim * (theta_C_sim - yC_sim[t-1, ]) * dt_sim + eta_C_sim * sqrt(dt_sim) * Z_corr[3, ]
  }
  
  # f. Default Times & Exposures
  int_lamB <- apply(pmax(sweep(yB_sim[1:n_steps_sim, ], 1, head(phi_B_vec, -1), "+"), 0), 2, cumsum) * dt_sim
  int_lamC <- apply(pmax(sweep(yC_sim[1:n_steps_sim, ], 1, head(phi_C_vec, -1), "+"), 0), 2, cumsum) * dt_sim
  
  E_B <- rexp(n_paths_sensi); E_C <- rexp(n_paths_sensi)
  
  CVA_paths <- numeric(n_paths_sensi)
  DVA_paths <- numeric(n_paths_sensi)
  
  for (p in 1:n_paths_sensi) {
    idx_B <- which(int_lamB[, p] > E_B[p])[1]
    idx_C <- which(int_lamC[, p] > E_C[p])[1]
    
    tau_B <- ifelse(is.na(idx_B), Inf, time_grid_sim[idx_B + 1])
    tau_C <- ifelse(is.na(idx_C), Inf, time_grid_sim[idx_C + 1])
    tau_min <- min(tau_B, tau_C)
    
    if (tau_min <= T_sim) {
      if (tau_C <= tau_B) { # CVA
        Vt <- S_sim[idx_C + 1, p] - K_sim * exp(-r_sim * (T_sim - tau_C))
        if (Vt > 0) CVA_paths[p] <- LGD_C_sim * Vt * exp(-r_sim * tau_C)
      }
      if (tau_B <= tau_C) { # DVA
        Vt <- S_sim[idx_B + 1, p] - K_sim * exp(-r_sim * (T_sim - tau_B))
        if (Vt < 0) DVA_paths[p] <- LGD_B_sim * abs(Vt) * exp(-r_sim * tau_B)
      }
    }
  }
  
  # g. Aggregation
  V0_clean_sim <- S0 - K_sim * exp(-r_sim * T_sim)
  CVA_sim <- mean(CVA_paths)
  DVA_sim <- mean(DVA_paths)
  
  return(c(Clean_Price = V0_clean_sim, CVA = CVA_sim, DVA = DVA_sim, Adjusted_Price = V0_clean_sim - CVA_sim + DVA_sim))
}

# 2. Run Scenarios & Build Table

results_sensi <- rbind(
  "Base case"                        = compute_Full_Sensi(),
  "Strike K +5%"                     = compute_Full_Sensi(K_bump = 0.05),
  "Maturity T +1 year"               = compute_Full_Sensi(T_bump = 1.00),
  "Risk-free rate r +100 bps"        = compute_Full_Sensi(r_bump = 0.01),
  "Correlation rho(S,B) +0.10"       = compute_Full_Sensi(rho_SB_bump = 0.10),
  "Correlation rho(S,C) +0.10"       = compute_Full_Sensi(rho_SC_bump = 0.10),
  "Correlation rho(B,C) +0.10"       = compute_Full_Sensi(rho_BC_bump = 0.10),
  "Stock drift mu_S +10%"            = compute_Full_Sensi(mu_S_bump = 0.10),
  "Stock volatility sigma_S +10%"    = compute_Full_Sensi(sigma_S_bump = 0.10),
  "Bank kappa_B +10%"                = compute_Full_Sensi(kappa_B_bump = 0.10),
  "Bank theta_B +10%"                = compute_Full_Sensi(theta_B_bump = 0.10),
  "Bank eta_B +10%"                  = compute_Full_Sensi(eta_B_bump = 0.10),
  "Counterparty kappa_C +10%"        = compute_Full_Sensi(kappa_C_bump = 0.10),
  "Counterparty theta_C +10%"        = compute_Full_Sensi(theta_C_bump = 0.10),
  "Counterparty eta_C +10%"          = compute_Full_Sensi(eta_C_bump = 0.10),
  "Bank recovery R_B +10 pp"         = compute_Full_Sensi(R_B_bump = 0.10),
  "Counterparty recovery R_C +10 pp" = compute_Full_Sensi(R_C_bump = 0.10)
)

results_df <- as.data.frame(results_sensi)

# Calculate Deltas relative to Base Case
base_vals <- results_df["Base case", ]
results_df$Delta_CVA_pct <- sprintf("%+.2f%%", (results_df$CVA / base_vals$CVA - 1) * 100)
results_df$Delta_DVA_pct <- sprintf("%+.2f%%", (results_df$DVA / base_vals$DVA - 1) * 100)
results_df$Delta_Adj_Price_pct <- sprintf("%+.2f%%", (results_df$Adjusted_Price / base_vals$Adjusted_Price - 1) * 100)

# Round raw prices for clean display
results_df[, 1:4] <- round(results_df[, 1:4], 4)

print(results_df)