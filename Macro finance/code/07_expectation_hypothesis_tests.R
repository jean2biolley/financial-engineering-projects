# Author: Jean de Biolley #

# The following code tests the Expectations Hypothesis of the term structure of interest rates
# using Fama-McBeth regressions and moving-block bootstrap resampling.


library(readxl)
library(dplyr)
library(lmtest)
library(sandwich)

zcb_raw <- read_excel("raw data/LW_monthly.xlsx", col_names = FALSE, col_types = "text")

# Locate the row containing the text "1 m"
anchor_idx <- which(apply(zcb_raw, 1, function(row) any(grepl("1 m", row))))[1]

# Convert the anchor row to character headers and assign "date" to the first element
headers <- as.character(zcb_raw[anchor_idx, ])
headers[1] <- "date"

# Subset everything below the anchor row and apply the extracted headers
zcb <- zcb_raw[(anchor_idx + 1):nrow(zcb_raw), ]
colnames(zcb) <- headers

# Filter out empty rows based on the date column and cast all variables to numeric
zcb <- zcb %>%
  filter(!is.na(date)) %>%
  mutate(across(everything(), as.numeric))


# Set up regression and bootstrap configurations
n_values <- c(2, 3, 4, 6, 9, 12, 24, 36, 48)
m_values <- c(1, 2, 3, 4, 6, 12)
B <- 1000
block_size <- 12

results <- data.frame()

for (n in n_values) {
  for (m in m_values) {
    if (n <= m) next

    # Calculate the remaining maturity after m periods
    rem_mat <- n - m

    # Skip iteration if maturity values exceed data dimensions
    if (n > ncol(zcb) | m > ncol(zcb) | rem_mat > ncol(zcb)) next

    # Pull out yield series as numeric vectors using element selection
    y_n <- as.numeric(zcb[[n]])
    y_m <- as.numeric(zcb[[m]])
    y_rem <- as.numeric(zcb[[rem_mat]])

    # Construct the dependent variable representing long yield changes
    dep_var <- as.numeric(dplyr::lead(y_rem, m) - y_n)

    # Construct the independent variable representing the scaled yield spread
    indep_var <- as.numeric((m / rem_mat) * (y_n - y_m))

    # Pack variables into a dataframe and drop terminal rows containing missing lead elements
    data_reg <- data.frame(Y = dep_var, X = indep_var)
    data_reg <- na.omit(data_reg)

    if (nrow(data_reg) == 0) next

    # Estimate the primary baseline OLS model
    fit <- lm(Y ~ X, data = data_reg)
    alpha1_hat <- coef(fit)["X"]

    # Compute Newey-West standard errors adjusted for moving average terms
    nw_vcov <- NeweyWest(fit, lag = m - 1, prewhite = FALSE, adjust = TRUE)
    se_alpha1_hat <- sqrt(diag(nw_vcov))["X"]

    # Calculate the empirical t-statistic testing the null hypothesis of equality to 1
    t_stat_hat <- (alpha1_hat - 1) / se_alpha1_hat

    # Initialize containers for the moving block bootstrap distributions
    n_obs <- nrow(data_reg)
    alpha1_boot <- numeric(B)
    t_boot <- numeric(B)

    # Enforce reproducibility for the randomized bootstrap resampling loop
    set.seed(123)
    for (b in 1:B) {
      # Resample random historical blocks with replacement
      start_indices <- sample(1:(n_obs - block_size + 1), ceiling(n_obs / block_size), replace = TRUE)
      boot_idx <- unlist(lapply(start_indices, function(i) i:(i + block_size - 1)))
      boot_idx <- boot_idx[1:n_obs]

      boot_data <- data_reg[boot_idx, ]

      # Estimate the OLS parameters on the bootstrapped data block
      boot_fit <- lm(Y ~ X, data = boot_data)
      a1_b <- coef(boot_fit)["X"]

      # Extract bootstrap standard errors while catching potential matrix singularities
      nw_vcov_b <- tryCatch(
        {
          NeweyWest(boot_fit, lag = m - 1, prewhite = FALSE, adjust = TRUE)
        },
        error = function(e) matrix(NA, 2, 2)
      )

      se_a1_b <- sqrt(diag(nw_vcov_b))["X"]

      alpha1_boot[b] <- a1_b

      # Shift the bootstrap test statistics to center around the empirical point estimate
      t_boot[b] <- (a1_b - alpha1_hat) / se_a1_b
    }

    # Remove records reflecting failed bootstrap iterations
    valid_boot <- !is.na(t_boot) & !is.na(alpha1_boot)
    alpha1_boot <- alpha1_boot[valid_boot]
    t_boot <- t_boot[valid_boot]

    # Check the proportion of bootstrap estimates extending further away from unity
    if (alpha1_hat < 1) {
      frac_dir <- mean(alpha1_boot < alpha1_hat)
    } else {
      frac_dir <- mean(alpha1_boot > alpha1_hat)
    }

    # Check the proportion of bootstrap t-statistics exceeding the empirical boundary
    frac_t <- mean(abs(t_boot) > abs(t_stat_hat))

    # Append calculated indicators to the final collection dataframe
    results <- rbind(results, data.frame(
      n = n,
      m = m,
      alpha1_hat = alpha1_hat,
      SE_alpha1 = se_alpha1_hat,
      Frac_dir = frac_dir,
      Frac_t = frac_t
    ))
  }
}

# Standardize output tables and map significance indicators
results$p_val <- 2 * (1 - pnorm(abs((results$alpha1_hat - 1) / results$SE_alpha1)))
results$stars <- ifelse(results$p_val < 0.01, "***",
  ifelse(results$p_val < 0.05, "**",
    ifelse(results$p_val < 0.1, "*", "")
  )
)

results$alpha1_hat <- round(results$alpha1_hat, 3)
results$alpha1_str <- paste0(results$alpha1_hat, results$stars)
results$SE_alpha1 <- round(results$SE_alpha1, 4)
results$Frac_dir <- round(results$Frac_dir, 3)
results$Frac_t <- round(results$Frac_t, 3)

final_table <- results[, c("n", "m", "alpha1_str", "SE_alpha1", "Frac_dir", "Frac_t")]
colnames(final_table) <- c("n", "m", "alpha1_hat", "SE(alpha1_hat)", "Frac_dir", "Frac_t")

cat("\n--- Test 1 Results: Long Yield Changes Over Life of Short Bond ---\n")
print(final_table, row.names = FALSE)


# Reset structural data collection for the second experiment layout
results_t2 <- data.frame()

for (n in n_values) {
  for (m in m_values) {
    # Isolate sample matches mapping to integer multiple maturities
    if (n <= m || n %% m != 0) next

    if (n > ncol(zcb) | m > ncol(zcb)) next

    k <- n / m

    y_n <- as.numeric(zcb[[n]])
    y_m <- as.numeric(zcb[[m]])

    # Accumulate leading values to evaluate perfect-foresight short yield paths
    future_sum <- y_m
    for (i in 1:(k - 1)) {
      future_sum <- future_sum + as.numeric(dplyr::lead(zcb[[m]], i * m))
    }
    dep_var <- (future_sum / k) - y_m

    # Measure the baseline current spread structure
    indep_var <- y_n - y_m

    data_reg <- data.frame(Y = dep_var, X = indep_var)
    data_reg <- na.omit(data_reg)

    if (nrow(data_reg) == 0) next

    # Fit the second baseline OLS structure
    fit <- lm(Y ~ X, data = data_reg)
    beta1_hat <- coef(fit)["X"]

    # Estimate long-horizon Newey-West standard errors to clean overlap biases
    nw_vcov <- NeweyWest(fit, lag = n - m, prewhite = FALSE, adjust = TRUE)
    se_beta1_hat <- sqrt(diag(nw_vcov))["X"]

    t_stat_hat <- (beta1_hat - 1) / se_beta1_hat

    n_obs <- nrow(data_reg)
    beta1_boot <- numeric(B)
    t_boot <- numeric(B)

    set.seed(123)
    for (b in 1:B) {
      start_indices <- sample(1:(n_obs - block_size + 1), ceiling(n_obs / block_size), replace = TRUE)
      boot_idx <- unlist(lapply(start_indices, function(i) i:(i + block_size - 1)))
      boot_idx <- boot_idx[1:n_obs]

      boot_data <- data_reg[boot_idx, ]

      boot_fit <- lm(Y ~ X, data = boot_data)
      b1_b <- coef(boot_fit)["X"]

      nw_vcov_b <- tryCatch(
        {
          NeweyWest(boot_fit, lag = n - m, prewhite = FALSE, adjust = TRUE)
        },
        error = function(e) matrix(NA, 2, 2)
      )

      se_b1_b <- sqrt(diag(nw_vcov_b))["X"]

      beta1_boot[b] <- b1_b
      t_boot[b] <- (b1_b - beta1_hat) / se_b1_b
    }

    valid_boot <- !is.na(t_boot) & !is.na(beta1_boot)
    beta1_boot <- beta1_boot[valid_boot]
    t_boot <- t_boot[valid_boot]

    if (beta1_hat < 1) {
      frac_dir <- mean(beta1_boot < beta1_hat)
    } else {
      frac_dir <- mean(beta1_boot > beta1_hat)
    }

    frac_t <- mean(abs(t_boot) > abs(t_stat_hat))

    results_t2 <- rbind(results_t2, data.frame(
      n = n,
      m = m,
      beta1_hat = beta1_hat,
      SE_beta1 = se_beta1_hat,
      Frac_dir = frac_dir,
      Frac_t = frac_t
    ))
  }
}

# Clean and report summary analytics for the second testing framework
results_t2$p_val <- 2 * (1 - pnorm(abs((results_t2$beta1_hat - 1) / results_t2$SE_beta1)))
results_t2$stars <- ifelse(results_t2$p_val < 0.01, "***",
  ifelse(results_t2$p_val < 0.05, "**",
    ifelse(results_t2$p_val < 0.1, "*", "")
  )
)

results_t2$beta1_hat <- round(results_t2$beta1_hat, 3)
results_t2$beta1_str <- paste0(results_t2$beta1_hat, results_t2$stars)
results_t2$SE_beta1 <- round(results_t2$SE_beta1, 4)
results_t2$Frac_dir <- round(results_t2$Frac_dir, 3)
results_t2$Frac_t <- round(results_t2$Frac_t, 3)

final_table_t2 <- results_t2[, c("n", "m", "beta1_str", "SE_beta1", "Frac_dir", "Frac_t")]
colnames(final_table_t2) <- c("n", "m", "beta1_hat", "SE(beta1_hat)", "Frac_dir", "Frac_t")

cat("\n--- Test 2 Results: Short Yield Changes Over Life of Long Bond ---\n")
print(final_table_t2, row.names = FALSE)
