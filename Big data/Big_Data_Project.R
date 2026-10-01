# In this code, we perform the following:
# 1. Data Processing: Cleaning and standardizing 8 multi-asset datasets (10, 25, 48, 100 portfolios).
# 2. Efficient Frontiers: Tracing Markowitz hyperbolas, GMVP, and 1/N baselines using `quadprog`.
# 3. Regularization: Applying Ridge (L2), Lasso (L1), and Elastic Net penalties via cross-validation.
# 4. Spectral Shrinkage: Filtering covariance noise using PCA truncation and Stein shrinkage.
# 5. Tail Risk Optimization: Optimizing skewness, excess kurtosis, and Modified VaR using `fastICA`.
# 6. Out-of-Sample Evaluation: Analyzing turnover, weight stability, Sharpe ratios, and cumulative wealth.

library(tidyverse)
library(stringr)
library(lubridate)
library(MASS)       
library(quadprog)  
library(patchwork)

all_files <- c(
  "10_Industry_Portfolios.csv",
  "10_Industry_Portfolios_daily.csv",
  "25_Portfolios_5x5.csv",
  "25_Portfolios_5x5_Daily.csv",
  "48_Industry_Portfolios.csv",
  "48_Industry_Portfolios_daily.csv",
  "100_Portfolios_10x10.csv",
  "100_Portfolios_10x10_Daily.csv",
  "25_Portfolios_ME_Prior_12_2.csv",        
  "25_Portfolios_ME_Prior_12_2_Daily.csv"  
)

dir.create("Clean_Datasets", showWarnings = FALSE)

process_french_data <- function(file_name) {
  input_file_path  <- file.path("Datasets", file_name)
  output_file_path <- file.path("Clean_Datasets", file_name)
  
  if (!file.exists(input_file_path)) {
    warning(paste("File not found, skipping:", input_file_path))
    return(NULL)
  }
  
  is_daily <- str_detect(tolower(file_name), "daily")
  
  if (is_daily) {
    start_date <- 19700101
    end_date   <- 20251231
    date_regex <- "^\\d{8}$" 
  } else {
    start_date <- 197001
    end_date   <- 202512
    date_regex <- "^\\d{6}$" 
  }
  
  all_lines <- readLines(input_file_path, warn = FALSE)
  
  anchor_string <- if (is_daily) {
    "Average Value Weighted Returns -- Daily"
  } else {
    "Average Value Weighted Returns -- Monthly"
  }
  
  anchor_idx <- which(str_detect(all_lines, fixed(anchor_string)))
  if (length(anchor_idx) == 0) {
    warning(paste("Anchor not found in:", file_name, "- Skipping."))
    return(NULL)
  }
  start_idx <- anchor_idx[1] + 1 
  
  is_blank <- str_trim(all_lines) == ""
  blank_lines_after_anchor <- which(is_blank & seq_along(all_lines) > start_idx)
  
  end_idx <- if (length(blank_lines_after_anchor) == 0) length(all_lines) else min(blank_lines_after_anchor) - 1
  table_lines <- all_lines[start_idx:end_idx]
  
  df <- read_csv(paste(table_lines, collapse = "\n"), show_col_types = FALSE)
  
  clean_df <- df %>%
    rename(date = 1) %>%
    mutate(date = str_trim(as.character(date))) %>%
    filter(str_detect(date, date_regex)) %>%
    mutate(date = as.integer(date)) %>%
    filter(date >= start_date & date <= end_date) %>%
    mutate(across(-date, ~ {
      val <- as.numeric(.x)
      val[val <= -99] <- 0
      return(val / 100)  # Raw decimal percentages format
    }))
  
  write_csv(clean_df, output_file_path)
}

cat("\nProcessing raw CSV files...\n")
walk(all_files, process_french_data)

X10_monthly   <- read_csv("Clean_Datasets/10_Industry_Portfolios.csv", show_col_types = FALSE)
X10_daily     <- read_csv("Clean_Datasets/10_Industry_Portfolios_daily.csv", show_col_types = FALSE)
X25_monthly   <- read_csv("Clean_Datasets/25_Portfolios_5x5.csv", show_col_types = FALSE)
X25_daily     <- read_csv("Clean_Datasets/25_Portfolios_5x5_Daily.csv", show_col_types = FALSE)
X48_monthly   <- read_csv("Clean_Datasets/48_Industry_Portfolios.csv", show_col_types = FALSE)
X48_daily     <- read_csv("Clean_Datasets/48_Industry_Portfolios_daily.csv", show_col_types = FALSE)
X100_monthly  <- read_csv("Clean_Datasets/100_Portfolios_10x10.csv", show_col_types = FALSE)
X100_daily    <- read_csv("Clean_Datasets/100_Portfolios_10x10_Daily.csv", show_col_types = FALSE)

Mom25_monthly <- read_csv("Clean_Datasets/25_Portfolios_ME_Prior_12_2.csv", show_col_types = FALSE)
Mom25_daily   <- read_csv("Clean_Datasets/25_Portfolios_ME_Prior_12_2_Daily.csv", show_col_types = FALSE)

cat("\nAll 8 core portfolio return datasets loaded into workspace (Ignoring Risk-Free Rate)!\n")

calc_weights_ew <- function(data_mat) {
  N <- ncol(data_mat)
  rep(1 / N, N)
}

calc_weights_mvp <- function(data_mat) {
  N <- ncol(data_mat)
  cov_mat <- cov(data_mat)
  inv_cov <- ginv(cov_mat) 
  ones <- rep(1, N)
  
  numerator <- inv_cov %*% ones
  denominator <- as.numeric(t(ones) %*% inv_cov %*% ones)
  return(as.numeric(numerator / denominator))
}

calc_weights_mv <- function(data_mat) {
  N <- ncol(data_mat)
  mu <- colMeans(data_mat)
  cov_mat <- cov(data_mat)
  inv_cov <- ginv(cov_mat) 
  ones <- rep(1, N)
  
  numerator <- inv_cov %*% mu
  denominator <- as.numeric(t(ones) %*% inv_cov %*% mu)
  return(as.numeric(numerator / denominator))
}

calc_metrics <- function(returns_vector, strategy_name, portfolio_name, freq_multiplier) {
  ann_mean <- mean(returns_vector) * freq_multiplier
  ann_vol  <- sd(returns_vector) * sqrt(freq_multiplier)
  sharpe   <- if(ann_vol > 0) ann_mean / ann_vol else NA
  
  data.frame(
    Portfolio = portfolio_name,
    Strategy  = strategy_name,
    Ann_Mean  = round(ann_mean, 4),
    Ann_Vol   = round(ann_vol, 4),
    Sharpe    = round(sharpe, 4)
  )
}

calc_turnover <- function(w_mat) {
  if(nrow(w_mat) <= 1) return(0)
  diffs <- diff(w_mat)
  mean(rowSums(abs(diffs)))
}

evaluate_portfolio <- function(data_df, portfolio_name, frequency = "monthly") {
  
  if (frequency == "monthly") {
    freq_mult <- 12; window_size <- 120; roll_step <- 6
  } else {
    freq_mult <- 252; window_size <- 2520; roll_step <- 126
  }
  
  ret_data <- as.matrix(data_df[, -1])
  T_total  <- nrow(ret_data)
  N <- ncol(ret_data)
  
  w_is_ew  <- calc_weights_ew(ret_data)
  w_is_mvp <- calc_weights_mvp(ret_data)
  w_is_mv  <- calc_weights_mv(ret_data)
  
  is_ret_ew  <- ret_data %*% w_is_ew
  is_ret_mvp <- ret_data %*% w_is_mvp
  is_ret_mv  <- ret_data %*% w_is_mv
  
  oos_ret_ew  <- c()
  oos_ret_mvp <- c()
  oos_ret_mv  <- c()
  
  weights_mvp <- matrix(nrow=0, ncol=N)
  weights_mv  <- matrix(nrow=0, ncol=N)
  
  for (t_end in seq(window_size, T_total - 1, by = roll_step)) {
    t_start    <- t_end - window_size + 1
    train_data <- ret_data[t_start:t_end, , drop = FALSE]
    
    w_oos_ew  <- calc_weights_ew(train_data)
    w_oos_mvp <- calc_weights_mvp(train_data)
    w_oos_mv  <- calc_weights_mv(train_data)
    
    weights_mvp <- rbind(weights_mvp, w_oos_mvp)
    weights_mv  <- rbind(weights_mv, w_oos_mv)
    
    test_start <- t_end + 1
    test_end   <- min(t_end + roll_step, T_total)
    test_data  <- ret_data[test_start:test_end, , drop = FALSE]
    
    oos_ret_ew  <- c(oos_ret_ew,  test_data %*% w_oos_ew)
    oos_ret_mvp <- c(oos_ret_mvp, test_data %*% w_oos_mvp)
    oos_ret_mv  <- c(oos_ret_mv,  test_data %*% w_oos_mv)
  }
  
  res_is <- bind_rows(
    calc_metrics(is_ret_ew,  "IS: 1/N", portfolio_name, freq_mult),
    calc_metrics(is_ret_mvp, "IS: Min-Variance", portfolio_name, freq_mult),
    calc_metrics(is_ret_mv,  "IS: Mean-Variance", portfolio_name, freq_mult)
  )
  
  res_oos <- bind_rows(
    calc_metrics(oos_ret_ew,  "OOS: 1/N", portfolio_name, freq_mult),
    calc_metrics(oos_ret_mvp, "OOS: Min-Variance", portfolio_name, freq_mult),
    calc_metrics(oos_ret_mv,  "OOS: Mean-Variance", portfolio_name, freq_mult)
  )
  
  turnover_res <- data.frame(
    Portfolio = portfolio_name,
    EW_Turnover = 0,
    MinVar_Turnover = round(calc_turnover(weights_mvp), 4),
    MeanVar_Turnover = round(calc_turnover(weights_mv), 4)
  )
  
  return(list(Results = bind_rows(res_is, res_oos), 
              Turnover = turnover_res, 
              Weights_MV = weights_mv))
}

cat("\nEvaluating Baseline Strategies across data tracks...\n")

res_10_m  <- evaluate_portfolio(X10_monthly, "10 Ind Monthly", "monthly")
res_10_d  <- evaluate_portfolio(X10_daily, "10 Ind Daily", "daily")
res_25_m  <- evaluate_portfolio(X25_monthly, "25 Size Monthly", "monthly")
res_25_d  <- evaluate_portfolio(X25_daily, "25 Size Daily", "daily")
res_48_m  <- evaluate_portfolio(X48_monthly, "48 Ind Monthly", "monthly")
res_48_d  <- evaluate_portfolio(X48_daily, "48 Ind Daily", "daily")
res_100_m <- evaluate_portfolio(X100_monthly, "100 Size Monthly", "monthly")
res_100_d <- evaluate_portfolio(X100_daily, "100 Size Daily", "daily")

baseline_perf <- bind_rows(res_10_m$Results, res_10_d$Results, res_25_m$Results, res_25_d$Results,
                           res_48_m$Results, res_48_d$Results, res_100_m$Results, res_100_d$Results)
print(baseline_perf)

cat("\nComputing true Markowitz hyperbola Efficient Frontiers...\n")

library(quadprog)
library(ggplot2)
library(dplyr)

compute_efficient_frontier <- function(data_df, frequency, universe_name, n_points = 150) {
  ret_data <- as.matrix(data_df[, -1])
  N <- ncol(ret_data)
  
  freq_mult <- if (frequency == "monthly") 12 else 252
  
  mu     <- colMeans(ret_data) * freq_mult
  sigma  <- cov(ret_data) * freq_mult
  
  if (rcond(sigma) < 1e-12) {
    sigma <- sigma + 1e-6 * diag(diag(sigma))
  }
  
  ones <- rep(1, N)
  inv_sigma <- solve(sigma)
  w_gmv <- as.numeric((inv_sigma %*% ones) / as.numeric(t(ones) %*% inv_sigma %*% ones))
  
  gmv_return <- as.numeric(t(w_gmv) %*% mu)
  gmv_risk   <- sqrt(as.numeric(t(w_gmv) %*% sigma %*% w_gmv))
  
  w_ew <- rep(1 / N, N)
  ew_return <- as.numeric(t(w_ew) %*% mu)
  ew_risk   <- sqrt(as.numeric(t(w_ew) %*% sigma %*% w_ew))
  
  return_spread <- max(mu) - min(mu)
  min_target <- gmv_return - (return_spread * 1.5)
  max_target <- gmv_return + (return_spread * 2.0)
  target_returns <- seq(min_target, max_target, length.out = n_points)
  
  frontier_risks <- numeric(n_points)
  Amat <- cbind(ones, mu)
  
  for (i in seq_along(target_returns)) {
    bvec <- c(1, target_returns[i])
    qp_solution <- tryCatch({
      solve.QP(Dmat = sigma, dvec = rep(0, N), Amat = Amat, bvec = bvec, meq = 2)
    }, error = function(e) { NULL })
    
    if (!is.null(qp_solution)) {
      frontier_risks[i] <- sqrt(2 * qp_solution$value)
    } else {
      frontier_risks[i] <- NA
    }
  }
  
  frontier_df <- data.frame(
    Risk = frontier_risks, Return = target_returns, Universe = universe_name, Freq = tools::toTitleCase(frequency)
  ) %>% 
    filter(!is.na(Risk)) %>%
    mutate(Efficiency = ifelse(Return >= gmv_return, "Efficient", "Inefficient"))
  
  gmv_point <- data.frame(Risk = gmv_risk, Return = gmv_return, Universe = universe_name, Freq = tools::toTitleCase(frequency))
  ew_point  <- data.frame(Risk = ew_risk, Return = ew_return, Universe = universe_name, Freq = tools::toTitleCase(frequency))
  asset_points <- data.frame(Risk = sqrt(diag(sigma)), Return = mu, Universe = universe_name, Freq = tools::toTitleCase(frequency))
  
  return(list(Frontier = frontier_df, GMV = gmv_point, EW = ew_point, Assets = asset_points))
}

target_tracks <- list(
  list(df = X10_monthly,  f = "monthly", n = "10 Industry"),
  list(df = X10_daily,    f = "daily",   n = "10 Industry"),
  list(df = X25_monthly,  f = "monthly", n = "25 Portfolios"),
  list(df = X25_daily,    f = "daily",   n = "25 Portfolios"),
  list(df = X48_monthly,  f = "monthly", n = "48 Industry"),
  list(df = X48_daily,    f = "daily",   n = "48 Industry"),
  list(df = X100_monthly, f = "monthly", n = "100 Portfolios"),
  list(df = X100_daily,   f = "daily",   n = "100 Portfolios")
)

for (track in target_tracks) {
  res <- compute_efficient_frontier(track$df, track$f, track$n)
  
  p <- ggplot() +

    geom_path(data = filter(res$Frontier, Efficiency == "Efficient"), aes(x = Risk, y = Return), color = "#1B3A61", linewidth = 1.2) +

    geom_path(data = filter(res$Frontier, Efficiency == "Inefficient"), aes(x = Risk, y = Return), color = "#7F8C8D", linewidth = 0.8, linetype = "dotted") +

    geom_point(data = res$Assets, aes(x = Risk, y = Return), color = "#B0BEC5", alpha = 0.6, size = 2) +
    
    geom_point(data = res$GMV, aes(x = Risk, y = Return), fill = "#7B1113", color = "white", size = 5.5, shape = 23, stroke = 0.8) + 
    geom_point(data = res$EW, aes(x = Risk, y = Return), fill = "#0F5257", color = "white", size = 5.5, shape = 22, stroke = 0.8) +   
    
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste(track$n, " (", tools::toTitleCase(track$f), ")", sep = ""),
      x = "Annualized Volatility (Risk)",
      y = "Annualized Expected Return"
    ) +

    theme_bw(base_size = 20) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 26, color = "#2C3E50"), 
      axis.title = element_text(size = 20, face = "bold", color = "#2C3E50"), 
      axis.text = element_text(size = 18, color = "#34495E"), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#BDC3C7", linewidth = 1),
      text = element_text(family = "serif") 
    )
  
  clean_filename <- paste0("Efficient_Frontier_", gsub(" ", "_", track$n), "_", tools::toTitleCase(track$f), ".png")
  ggsave(clean_filename, plot = p, width = 8, height = 5.5, dpi = 300)
}

cat("\nAll 8 institutional-grade Efficient Frontier plots saved successfully!\n")

cat("\nGenerating Out-of-Sample portfolio weight stability boxplots...\n")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

extract_oos_weights <- function(data_df, frequency) {
  if (frequency == "monthly") {
    window_size <- 120; roll_step <- 6
  } else {
    window_size <- 2520; roll_step <- 126
  }
  
  local_mat <- as.matrix(data_df[, -1])
  T_total    <- nrow(local_mat)
  N          <- ncol(local_mat)
  asset_names <- colnames(local_mat)
  
  weights_ew_list  <- list()
  weights_mvp_list <- list()
  weights_mv_list  <- list()
  window_idx       <- 1
  
  for (t_end in seq(window_size, T_total - 1, by = roll_step)) {
    t_start    <- t_end - window_size + 1
    train_data <- local_mat[t_start:t_end, , drop = FALSE]

    w_ew  <- as.numeric(calc_weights_ew(train_data))
    w_mvp <- as.numeric(calc_weights_mvp(train_data))
    w_mv  <- as.numeric(calc_weights_mv(train_data))
    
    weights_ew_list[[window_idx]]  <- data.frame(Window = window_idx, Strategy = "Equally Weighted", Asset = asset_names, Weight = w_ew)
    weights_mvp_list[[window_idx]] <- data.frame(Window = window_idx, Strategy = "Minimum Variance", Asset = asset_names, Weight = w_mvp)
    weights_mv_list[[window_idx]]  <- data.frame(Window = window_idx, Strategy = "Mean-Variance", Asset = asset_names, Weight = w_mv)
    
    window_idx <- window_idx + 1
  }
  
  all_weights <- bind_rows(weights_ew_list, weights_mvp_list, weights_mv_list) %>%
    mutate(Freq = tools::toTitleCase(frequency))
  
  return(all_weights)
}

dataset_list <- list(
  "10 Industry"   = list(m_df = X10_monthly,  d_df = X10_daily,   type = "detailed"),
  "25 Portfolios" = list(m_df = X25_monthly,  d_df = X25_daily,   type = "summary"),
  "48 Industry"   = list(m_df = X48_monthly,  d_df = X48_daily,   type = "summary"),
  "100 Portfolios"= list(m_df = X100_monthly, d_df = X100_daily,  type = "summary")
)

for (name in names(dataset_list)) {
  cat(paste("Processing weight distributions for:", name, "...\n"))
  
  df_w_daily   <- extract_oos_weights(dataset_list[[name]]$d_df, "daily")
  df_w_monthly <- extract_oos_weights(dataset_list[[name]]$m_df, "monthly")
  
  combined_w <- bind_rows(df_w_daily, df_w_monthly)
  combined_w$Strategy <- factor(combined_w$Strategy, levels = c("Equally Weighted", "Mean-Variance", "Minimum Variance"))
  
  plot_generator <- function(filtered_df, subtitle_text) {
    p <- ggplot(filtered_df, aes(fill = Strategy))
    
    p <- p + 
      facet_wrap(~Strategy, scales = "free_x", ncol = 3) + 
      scale_fill_manual(values = c("Equally Weighted" = "#8DA0CB", "Mean-Variance" = "#FC8D62", "Minimum Variance" = "#66C2A5")) + 
      coord_cartesian(ylim = c(-1.5, 2.5)) + 
      labs(title = paste(name, "-", subtitle_text), y = "Portfolio Allocation Weight") +
      theme_bw(base_size = 20) + 
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 24), 
        axis.title.y = element_text(size = 20, face = "bold"),            
        axis.text.y = element_text(size = 16),
        strip.background = element_rect(fill = "#F0F0F0", color = "#CCCCCC"), 
        strip.text = element_text(face = "bold", color = "#333333", size = 16),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(family = "serif")
      )
    
    if (dataset_list[[name]]$type == "detailed") {
      p <- p + geom_boxplot(aes(x = Asset, y = Weight), outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8, color = "#2C3E50", linewidth = 0.3) +
        labs(x = "Industry") + 
        theme(
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 13)
        ) 
    } else {
      p <- p + geom_boxplot(aes(x = Strategy, y = Weight), outlier.size = 0.5, outlier.alpha = 0.3, alpha = 0.8, color = "#2C3E50", linewidth = 0.3) +
        labs(x = NULL) + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
    }
    
    return(p)
  }
  
  p_d <- plot_generator(filter(combined_w, Freq == "Daily"), "Daily OOS Weights")
  p_m <- plot_generator(filter(combined_w, Freq == "Monthly"), "Monthly OOS Weights")
  
  ggsave(paste0("Weight_Stability_", gsub(" ", "_", name), "_Daily.png"), plot = p_d, width = 8, height = 5, dpi = 300)
  ggsave(paste0("Weight_Stability_", gsub(" ", "_", name), "_Monthly.png"), plot = p_m, width = 8, height = 5, dpi = 300)
}

cat("\nAll individual stability boxplots successfully generated and saved!\n")

cat("\nCompiling all 8 Global Minimum-Variance portfolio weights into a unified DataFrame...\n")

library(dplyr)
library(purrr)

extract_gmv_weights_df <- function(data_df, frequency, universe_name) {
  ret_data <- as.matrix(data_df[, -1])
  N <- ncol(ret_data)
  asset_names <- colnames(ret_data)
  
  freq_mult <- if (frequency == "monthly") 12 else 252
  sigma     <- cov(ret_data) * freq_mult
  
  if (rcond(sigma) < 1e-12) {
    sigma <- sigma + 1e-6 * diag(diag(sigma))
  }
  
  ones <- rep(1, N)
  inv_sigma <- solve(sigma)
  w_gmv <- as.numeric((inv_sigma %*% ones) / as.numeric(t(ones) %*% inv_sigma %*% ones))

  data.frame(
    Universe      = universe_name,
    Frequency     = tools::toTitleCase(frequency),
    Asset_ID      = seq_len(N),
    Asset_Name    = asset_names,
    GMV_Weight    = round(w_gmv, 6) 
  )
}

gmv_tracks <- list(
  list(df = X10_monthly,  f = "monthly", n = "10 Industry"),
  list(df = X10_daily,    f = "daily",   n = "10 Industry"),
  list(df = X25_monthly,  f = "monthly", n = "25 Portfolios"),
  list(df = X25_daily,    f = "daily",   n = "25 Portfolios"),
  list(df = X48_monthly,  f = "monthly", n = "48 Industry"),
  list(df = X48_daily,    f = "daily",   n = "48 Industry"),
  list(df = X100_monthly, f = "monthly", n = "100 Portfolios"),
  list(df = X100_daily,   f = "daily",   n = "100 Portfolios")
)

gmw_weights <- map_df(gmv_tracks, function(track) {
  extract_gmv_weights_df(track$df, track$f, track$n)
})


calc_weights_l2_q23 <- function(data_mat, gamma) {
  N <- ncol(data_mat)
  sigma_hat <- cov(data_mat)
  sigma_tilde <- sigma_hat + gamma * diag(N)
  ones <- rep(1, N)
  
  inv_sigma_tilde <- tryCatch(
    solve(sigma_tilde),
    error = function(e) MASS::ginv(sigma_tilde)
  )
  
  numerator <- inv_sigma_tilde %*% ones
  denominator <- as.numeric(t(ones) %*% inv_sigma_tilde %*% ones)
  
  as.numeric(numerator / denominator)
}

make_gamma_grid_q23 <- function(data_mat, n_points = 12) {
  sigma_hat <- cov(data_mat)
  scale_gamma <- mean(diag(sigma_hat))
  scale_gamma <- max(scale_gamma, 1e-8)
  
  gamma_grid <- scale_gamma * 10^seq(-4, 1, length.out = n_points)
  unique(sort(gamma_grid))
}

make_block_folds_q23 <- function(T_train, k_folds) {
  cut_points <- floor(seq(0, T_train, length.out = k_folds + 1))
  folds <- vector("list", k_folds)
  
  for (k in seq_len(k_folds)) {
    start_idx <- cut_points[k] + 1
    end_idx <- cut_points[k + 1]
    folds[[k]] <- start_idx:end_idx
  }
  
  folds
}

cv_score_gamma_l2_q23 <- function(train_data, gamma, k_folds = NULL) {
  T_train <- nrow(train_data)
  N <- ncol(train_data)
  
  if (is.null(k_folds)) {
    k_folds <- ifelse(T_train <= 240, 5, 10)
  }
  
  folds <- make_block_folds_q23(T_train, k_folds)
  cv_returns <- c()
  
  for (fold_idx in seq_along(folds)) {
    val_idx <- folds[[fold_idx]]
    train_idx <- setdiff(seq_len(T_train), val_idx)

    
    train_fold <- train_data[train_idx, , drop = FALSE]
    val_fold <- train_data[val_idx, , drop = FALSE]
    
    w_cv <- calc_weights_l2_q23(train_fold, gamma)
    cv_returns <- c(cv_returns, as.numeric(val_fold %*% w_cv))
  }
  
  if (length(cv_returns) <= 1) return(Inf)
  
  var(cv_returns)
}

find_optimal_gamma_q23 <- function(train_data, gamma_grid, k_folds = NULL) {
  cv_scores <- sapply(gamma_grid, function(g) {
    cv_score_gamma_l2_q23(train_data, gamma = g, k_folds = k_folds)
  })
  
  best_idx <- which.min(cv_scores)
  
  list(
    best_gamma = gamma_grid[best_idx],
    cv_scores = data.frame(
      Gamma = gamma_grid,
      CV_Variance = cv_scores
    )
  )
}

evaluate_l2_fixed_gamma_q23 <- function(data_df, portfolio_name, gamma, frequency = "monthly") {
  
  if (frequency == "monthly") {
    freq_mult <- 12
    window_size <- 120
    roll_step <- 6
  } else {
    freq_mult <- 252
    window_size <- 2520
    roll_step <- 126
  }
  
  ret_data <- as.matrix(data_df[, -1])
  T_total <- nrow(ret_data)
  N <- ncol(ret_data)
  
  oos_ret_l2 <- c()
  weights_l2 <- matrix(nrow = 0, ncol = N)
  weights_mv <- matrix(nrow = 0, ncol = N)
  
  for (t_end in seq(window_size, T_total - 1, by = roll_step)) {
    t_start <- t_end - window_size + 1
    train_data <- ret_data[t_start:t_end, , drop = FALSE]
    
    w_l2 <- calc_weights_l2_q23(train_data, gamma)
    w_mv <- calc_weights_mvp(train_data)
    
    weights_l2 <- rbind(weights_l2, w_l2)
    weights_mv <- rbind(weights_mv, w_mv)
    
    test_start <- t_end + 1
    test_end <- min(t_end + roll_step, T_total)
    test_data <- ret_data[test_start:test_end, , drop = FALSE]
    
    oos_ret_l2 <- c(oos_ret_l2, as.numeric(test_data %*% w_l2))
  }
  
  perf_res <- calc_metrics(
    returns_vector = oos_ret_l2,
    strategy_name = paste0("OOS: L2 MinVar (gamma=", signif(gamma, 3), ")"),
    portfolio_name = portfolio_name,
    freq_multiplier = freq_mult
  )
  
  weight_comp <- data.frame(
    Portfolio = portfolio_name,
    Gamma = gamma,
    AvgAbsWeight_L2 = mean(abs(weights_l2)),
    MaxAbsWeight_L2 = mean(apply(abs(weights_l2), 1, max)),
    MinWeight_L2 = mean(apply(weights_l2, 1, min)),
    MaxWeight_L2 = mean(apply(weights_l2, 1, max)),
    AvgAbsWeight_MV = mean(abs(weights_mv)),
    MaxAbsWeight_MV = mean(apply(abs(weights_mv), 1, max)),
    MinWeight_MV = mean(apply(weights_mv, 1, min)),
    MaxWeight_MV = mean(apply(weights_mv, 1, max))
  )
  
  list(
    Results = perf_res,
    Weight_Comparison = weight_comp,
    OOS_Returns = oos_ret_l2,
    Weights_L2 = weights_l2,
    Weights_MV = weights_mv
  )
}

evaluate_dynamic_l2_q23 <- function(data_df, portfolio_name, frequency = "monthly", k_folds = NULL) {
  
  if (frequency == "monthly") {
    freq_mult <- 12
    window_size <- 120
    roll_step <- 6
  } else {
    freq_mult <- 252
    window_size <- 2520
    roll_step <- 126
  }
  
  ret_data <- as.matrix(data_df[, -1])
  T_total <- nrow(ret_data)
  N <- ncol(ret_data)
  
  oos_ret_dyn <- c()
  weights_dyn <- matrix(nrow = 0, ncol = N)
  weights_mv <- matrix(nrow = 0, ncol = N)
  chosen_gammas <- c()
  cv_tables <- list()
  window_counter <- 1
  
  for (t_end in seq(window_size, T_total - 1, by = roll_step)) {
    t_start <- t_end - window_size + 1
    train_data <- ret_data[t_start:t_end, , drop = FALSE]
    
    gamma_grid <- make_gamma_grid_q23(train_data)
    gamma_selection <- find_optimal_gamma_q23(
      train_data = train_data,
      gamma_grid = gamma_grid,
      k_folds = k_folds
    )
    
    opt_gamma <- gamma_selection$best_gamma
    chosen_gammas <- c(chosen_gammas, opt_gamma)
    
    cv_table_window <- gamma_selection$cv_scores
    cv_table_window$Window <- window_counter
    cv_table_window$Portfolio <- portfolio_name
    cv_tables[[window_counter]] <- cv_table_window
    
    w_dyn <- calc_weights_l2_q23(train_data, opt_gamma)
    w_mv <- calc_weights_mvp(train_data)
    
    weights_dyn <- rbind(weights_dyn, w_dyn)
    weights_mv <- rbind(weights_mv, w_mv)
    
    test_start <- t_end + 1
    test_end <- min(t_end + roll_step, T_total)
    test_data <- ret_data[test_start:test_end, , drop = FALSE]
    
    oos_ret_dyn <- c(oos_ret_dyn, as.numeric(test_data %*% w_dyn))
    
    window_counter <- window_counter + 1
  }
  
  perf_res <- calc_metrics(
    returns_vector = oos_ret_dyn,
    strategy_name = "OOS: L2 MinVar (dynamic gamma by CV)",
    portfolio_name = portfolio_name,
    freq_multiplier = freq_mult
  )
  
  gamma_summary <- data.frame(
    Portfolio = portfolio_name,
    Avg_Gamma = mean(chosen_gammas),
    Median_Gamma = median(chosen_gammas),
    Min_Gamma = min(chosen_gammas),
    Max_Gamma = max(chosen_gammas)
  )
  
  weight_comp <- data.frame(
    Portfolio = portfolio_name,
    AvgAbsWeight_L2 = mean(abs(weights_dyn)),
    MaxAbsWeight_L2 = mean(apply(abs(weights_dyn), 1, max)),
    MinWeight_L2 = mean(apply(weights_dyn, 1, min)),
    MaxWeight_L2 = mean(apply(weights_dyn, 1, max)),
    AvgAbsWeight_MV = mean(abs(weights_mv)),
    MaxAbsWeight_MV = mean(apply(abs(weights_mv), 1, max)),
    MinWeight_MV = mean(apply(weights_mv, 1, min)),
    MaxWeight_MV = mean(apply(weights_mv, 1, max))
  )
  
  list(
    Results = perf_res,
    Gamma_Summary = gamma_summary,
    Weight_Comparison = weight_comp,
    Chosen_Gammas = chosen_gammas,
    CV_Table = bind_rows(cv_tables),
    OOS_Returns = oos_ret_dyn,
    Weights_L2 = weights_dyn,
    Weights_MV = weights_mv
  )
}

cat("\nGenerating L2 Regularization (Gamma) Shrinkage Plots...\n")
library(ggplot2)
library(dplyr)
library(tidyr)

q2_perf_monthly <- q2_performance_table %>% filter(grepl("Monthly", Portfolio))
q2_weights_monthly <- q2_weight_comparison_table %>% filter(grepl("Monthly", Portfolio))

p_sharpe <- ggplot(q2_perf_monthly, aes(x = Gamma, y = Sharpe, color = Portfolio)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10() + 
  labs(
    title = "OOS Sharpe Ratio vs. L2 Penalty (Monthly)",
    x = "Gamma (Log Scale)",
    y = "Annualized Sharpe Ratio"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 19),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Gamma_vs_Sharpe_Monthly.png", plot = p_sharpe, width = 10, height = 6, dpi = 300)

shrinkage_data <- q2_weights_monthly %>% filter(Portfolio == "100 Size Monthly")

p_shrinkage <- ggplot(shrinkage_data, aes(x = Gamma)) +

  geom_hline(aes(yintercept = MaxWeight_MV), linetype = "dashed", color = "#7B1113", linewidth = 1) +
  geom_hline(aes(yintercept = MinWeight_MV), linetype = "dashed", color = "#1B3A61", linewidth = 1) +
  
  geom_line(aes(y = MaxWeight_L2, color = "Max Long Position"), linewidth = 1.5) +
  geom_line(aes(y = MinWeight_L2, color = "Max Short Position"), linewidth = 1.5) +
  
  scale_x_log10() +
  scale_color_manual(values = c("Max Long Position" = "#7B1113", "Max Short Position" = "#1B3A61")) +
  labs(
    title = "Weight Shrinkage Effect (100 Portfolios - Monthly)",
    x = "Gamma (Log Scale)",
    y = "Portfolio Weight Allocation"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 19),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Gamma_Shrinkage_100_Monthly.png", plot = p_shrinkage, width = 10, height = 6, dpi = 300)

cat("\nL2 Regularization plots saved successfully to Q2_Q3_Output/ folder!\n")

q2_perf_daily <- q2_performance_table %>% filter(grepl("Daily", Portfolio))
q2_weights_daily <- q2_weight_comparison_table %>% filter(grepl("Daily", Portfolio))

p_sharpe_daily <- ggplot(q2_perf_daily, aes(x = Gamma, y = Sharpe, color = Portfolio)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_log10() + 
  labs(
    title = "OOS Sharpe Ratio vs. L2 Penalty (Daily)",
    x = "Gamma (Log Scale)",
    y = "Annualized Sharpe Ratio"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 19),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Gamma_vs_Sharpe_Daily.png", plot = p_sharpe_daily, width = 10, height = 6, dpi = 300)

shrinkage_data_daily <- q2_weights_daily %>% filter(Portfolio == "100 Size Daily")

p_shrinkage_daily <- ggplot(shrinkage_data_daily, aes(x = Gamma)) +

  geom_hline(aes(yintercept = MaxWeight_MV), linetype = "dashed", color = "#7B1113", linewidth = 1) +
  geom_hline(aes(yintercept = MinWeight_MV), linetype = "dashed", color = "#1B3A61", linewidth = 1) +
  
  geom_line(aes(y = MaxWeight_L2, color = "Max Long Position"), linewidth = 1.5) +
  geom_line(aes(y = MinWeight_L2, color = "Max Short Position"), linewidth = 1.5) +
  
  scale_x_log10() +
  scale_color_manual(values = c("Max Long Position" = "#7B1113", "Max Short Position" = "#1B3A61")) +
  labs(
    title = "Weight Shrinkage Effect (100 Portfolios - Daily)",
    x = "Gamma (Log Scale)",
    y = "Portfolio Weight Allocation"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 25),
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 19),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Gamma_Shrinkage_100_Daily.png", plot = p_shrinkage_daily, width = 10, height = 6, dpi = 300)

cat("\nDaily L2 Regularization plots saved successfully!\n")

cat("\nRunning Question 2: fixed gamma analysis...\n")

datasets_q23 <- list(
  list(data = X10_monthly,  name = "10 Ind Monthly",   freq = "monthly"),
  list(data = X10_daily,    name = "10 Ind Daily",     freq = "daily"),
  list(data = X25_monthly,  name = "25 Size Monthly",  freq = "monthly"),
  list(data = X25_daily,    name = "25 Size Daily",    freq = "daily"),
  list(data = X48_monthly,  name = "48 Ind Monthly",   freq = "monthly"),
  list(data = X48_daily,    name = "48 Ind Daily",     freq = "daily"),
  list(data = X100_monthly, name = "100 Size Monthly", freq = "monthly"),
  list(data = X100_daily,   name = "100 Size Daily",   freq = "daily")
)

q2_perf_all <- list()
q2_weights_all <- list()

for (ds in datasets_q23) {
  ret_mat <- as.matrix(ds$data[, -1])
  gamma_grid_ds <- make_gamma_grid_q23(ret_mat)
  
  for (g in gamma_grid_ds) {
    res_q2 <- evaluate_l2_fixed_gamma_q23(
      data_df = ds$data,
      portfolio_name = ds$name,
      gamma = g,
      frequency = ds$freq
    )
    
    q2_perf_all[[length(q2_perf_all) + 1]] <- cbind(
      data.frame(Gamma = g),
      res_q2$Results
    )
    
    q2_weights_all[[length(q2_weights_all) + 1]] <- res_q2$Weight_Comparison
  }
}

q2_performance_table <- bind_rows(q2_perf_all)
q2_weight_comparison_table <- bind_rows(q2_weights_all)

cat("\nQuestion 2 finished.\n")
print(q2_performance_table)
print(q2_weight_comparison_table)

cat("\nRunning Question 3: dynamic gamma chosen by cross-validation...\n")

q3_results_all <- list()
q3_gamma_all <- list()
q3_weights_all <- list()
q3_cv_all <- list()

for (ds in datasets_q23) {
  res_q3 <- evaluate_dynamic_l2_q23(
    data_df = ds$data,
    portfolio_name = ds$name,
    frequency = ds$freq,
    k_folds = NULL
  )
  
  q3_results_all[[length(q3_results_all) + 1]] <- res_q3$Results
  q3_gamma_all[[length(q3_gamma_all) + 1]] <- res_q3$Gamma_Summary
  q3_weights_all[[length(q3_weights_all) + 1]] <- res_q3$Weight_Comparison
  q3_cv_all[[length(q3_cv_all) + 1]] <- res_q3$CV_Table
}

q3_performance_table <- bind_rows(q3_results_all)
q3_gamma_summary_table <- bind_rows(q3_gamma_all)
q3_weight_comparison_table <- bind_rows(q3_weights_all)
q3_cv_table <- bind_rows(q3_cv_all)

cat("\nQuestion 3 finished.\n")
print(q3_performance_table)
print(q3_gamma_summary_table)
print(q3_weight_comparison_table)

dir.create("Q2_Q3_Output", showWarnings = FALSE)

write_csv(q2_performance_table, "Q2_Q3_Output/q2_performance_table.csv")
write_csv(q2_weight_comparison_table, "Q2_Q3_Output/q2_weight_comparison_table.csv")

write_csv(q3_performance_table, "Q2_Q3_Output/q3_performance_table.csv")
write_csv(q3_gamma_summary_table, "Q2_Q3_Output/q3_gamma_summary_table.csv")
write_csv(q3_weight_comparison_table, "Q2_Q3_Output/q3_weight_comparison_table.csv")
write_csv(q3_cv_table, "Q2_Q3_Output/q3_cv_table.csv")

cat("\nGenerating Q3 Cross-Validation Plots...\n")

cv_example_data <- q3_cv_table %>% 
  filter(Portfolio == "100 Size Monthly", Window == 1)

p_cv_curve <- ggplot(cv_example_data, aes(x = Gamma, y = CV_Variance)) +
  geom_line(color = "#1B3A61", linewidth = 1.2) +
  geom_point(color = "#7B1113", size = 4) +
  scale_x_log10() +
  labs(
    title = "Cross-Validation Score vs. Gamma (100 Size Monthly - Window 1)",
    x = "Gamma (Log Scale)",
    y = "CV OOS Variance"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    axis.title = element_text(face = "bold"),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Q3_CV_Curve_Example.png", plot = p_cv_curve, width = 10, height = 6, dpi = 300)

gamma_evolution_data <- q3_cv_table %>%
  group_by(Portfolio, Window) %>%
  slice(which.min(CV_Variance)) %>% 
  ungroup() %>%
  filter(grepl("Monthly", Portfolio)) 

p_gamma_time <- ggplot(gamma_evolution_data, aes(x = Window, y = Gamma, color = Portfolio)) +
  geom_step(linewidth = 1) + 
  scale_y_log10() +
  labs(
    title = "Evolution of Optimal Gamma over Time (Monthly Datasets)",
    x = "Rebalancing Window Index",
    y = "Chosen Optimal Gamma (Log Scale)"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(family = "serif")
  )

ggsave("Q2_Q3_Output/Q3_Dynamic_Gamma_Time.png", plot = p_gamma_time, width = 10, height = 6, dpi = 300)

cat("\nQ3 plots saved successfully to Q2_Q3_Output/ folder!\n")

add_months <- function(date, n_months) { seq(date, by = paste(n_months, "months"), length.out = 2)[2] }
add_years <- function(date, n_years) { seq(date, by = paste(n_years, "years"), length.out = 2)[2] }

make_calendar_rolling_windows <- function(dates, train_years = 10, test_months = 6, step_months = 6) {
  first_month <- as.Date(format(min(dates), "%Y-%m-01"))
  last_month_end <- add_months(as.Date(format(max(dates), "%Y-%m-01")), 1) - 1
  start_dates <- seq(first_month, max(dates), by = paste(step_months, "months"))
  windows <- list()
  
  for (i in seq_along(start_dates)) {
    train_start <- start_dates[i]
    train_end <- add_years(train_start, train_years) - 1
    test_start <- train_end + 1
    test_end <- add_months(test_start, test_months) - 1
    
    if (test_end > last_month_end) break
    
    windows[[length(windows) + 1]] <- data.frame(
      window_id = length(windows) + 1, train_start = train_start,
      train_end = train_end, test_start = test_start, test_end = test_end
    )
  }
  do.call(rbind, windows)
}

make_dynamic_gamma_grid <- function(returns_mat) {
  sigma <- stats::cov(returns_mat)
  sigma[!is.finite(sigma)] <- 0
  sigma <- (sigma + t(sigma)) / 2
  variance_scale <- median(diag(sigma), na.rm = TRUE)
  if (!is.finite(variance_scale) || variance_scale <= 0) variance_scale <- 1e-4
  
  nonzero_grid <- variance_scale * 10^seq(-4, 2, length.out = 8) 
  sort(unique(c(0, nonzero_grid)))
}

stabilize_covariance <- function(returns_mat) {
  sigma <- stats::cov(returns_mat)
  sigma[!is.finite(sigma)] <- 0
  sigma <- (sigma + t(sigma)) / 2
  scale <- mean(diag(sigma), na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- 1e-6
  sigma + diag(scale * 1e-8, nrow(sigma))
}

solve_penalized_minvar <- function(sigma, gamma, penalty_type, alpha = 0.5) {
  n_assets <- ncol(sigma)
  l1_weight <- if (penalty_type == "Lasso") 1 else alpha
  l2_weight <- if (penalty_type == "Lasso") 0 else (1 - alpha)
  
  sigma_scale <- mean(diag(sigma), na.rm = TRUE)
  if (!is.finite(sigma_scale) || sigma_scale <= 0) sigma_scale <- 1e-6
  
  n_variables <- 2 * n_assets
  dvec <- c(rep(0, n_assets), rep(-gamma * l1_weight, n_assets))

  equality_constraint <- c(rep(1, n_assets), rep(0, n_assets))
  upper_abs_constraints <- rbind(-diag(n_assets), diag(n_assets))
  lower_abs_constraints <- rbind(diag(n_assets), diag(n_assets))
  amat <- cbind(equality_constraint, upper_abs_constraints, lower_abs_constraints)
  bvec <- c(1, rep(0, 2 * n_assets))
  
  base_dmat <- matrix(0, n_variables, n_variables)
  base_dmat[1:n_assets, 1:n_assets] <- 2 * (sigma + gamma * l2_weight * diag(n_assets))
  
  for (multiplier in c(1, 1e3, 1e6, 1e9)) {
    eps <- sigma_scale * 1e-10 * multiplier
    dmat <- base_dmat + diag(eps, n_variables)
    solution <- tryCatch(solve.QP(Dmat = dmat, dvec = dvec, Amat = amat, bvec = bvec, meq = 1), error = function(e) NULL)
    if (!is.null(solution)) {
      weights <- solution$solution[seq_len(n_assets)]
      return(weights / sum(weights))
    }
  }
  stop("Quadprog solver failure.")
}

find_optimal_gamma_q4 <- function(train_data, gamma_grid, penalty_type, seed_val) {
  T_train <- nrow(train_data)
  n_assets <- ncol(train_data)
  k_folds <- 10
  
  set.seed(seed_val)
  shuffled_indices <- sample(seq_len(T_train))
  fold_ids <- integer(T_train)
  fold_ids[shuffled_indices] <- rep(seq_len(k_folds), length.out = T_train)
  
  best_gamma <- gamma_grid[1]
  best_avg_var <- Inf 
  
  for (g in gamma_grid) {
    fold_variances <- c()
    for (i in seq_len(k_folds)) {
      val_idx <- which(fold_ids == i)
      train_idx <- which(fold_ids != i)
      
      if (length(val_idx) < 2 || length(train_idx) < (n_assets + 1)) next
      
      tt_data <- train_data[train_idx, , drop = FALSE]
      val_data <- train_data[val_idx, , drop = FALSE]
      
      sigma_tt <- stabilize_covariance(tt_data)
      w_val <- tryCatch(solve_penalized_minvar(sigma_tt, g, penalty_type), error = function(e) NULL)
      if(is.null(w_val)) next
      
      val_returns <- val_data %*% w_val
      fold_variances <- c(fold_variances, stats::var(val_returns, na.rm = TRUE))
    }
    
    avg_var <- mean(fold_variances, na.rm = TRUE)
    if (!is.na(avg_var) && avg_var < best_avg_var) {
      best_avg_var <- avg_var
      best_gamma <- g
    }
  }
  return(best_gamma)
}

evaluate_upgraded_penalties <- function(raw_returns_df, portfolio_name, penalty_type, frequency = "monthly") {
  
  if (frequency == "monthly") {
    dates <- as.Date(paste0(raw_returns_df$date, "01"), format = "%Y%m%d")
    freq_m <- 12; t_years <- 10; t_months <- 6; s_months <- 6
  } else {
    dates <- as.Date(as.character(raw_returns_df$date), format = "%Y%m%d")
    freq_m <- 252; t_years <- 10; t_months <- 6; s_months <- 6
  }
  
  ret_data <- as.matrix(raw_returns_df[, -1])
  N <- ncol(ret_data)
  
  outer_windows <- make_calendar_rolling_windows(dates, train_years = t_years, test_months = t_months, step_months = s_months)
  oos_ret <- c()
  weights_mat <- matrix(nrow = 0, ncol = N)
  
  cat(sprintf("\nRunning %s Engine on %s (%s frequency)...\n", penalty_type, portfolio_name, frequency))
  
  for (i in seq_len(nrow(outer_windows))) {
    train_idx <- which(dates >= outer_windows$train_start[i] & dates <= outer_windows$train_end[i])
    test_idx  <- which(dates >= outer_windows$test_start[i]  & dates <= outer_windows$test_end[i])
    
    if(length(train_idx) == 0 || length(test_idx) == 0) next
    train_data <- ret_data[train_idx, , drop = FALSE]
    
    gamma_grid <- make_dynamic_gamma_grid(train_data)
    opt_gamma <- find_optimal_gamma_q4(train_data, gamma_grid, penalty_type, seed_val = 12345 + i)
    
    sigma_full <- stabilize_covariance(train_data)
    w_oos  <- solve_penalized_minvar(sigma_full, opt_gamma, penalty_type)
    weights_mat  <- rbind(weights_mat, w_oos)
    
    test_data  <- ret_data[test_idx, , drop = FALSE]
    oos_ret <- c(oos_ret, test_data %*% w_oos)
  }
  
  res_oos <- calc_metrics(oos_ret, paste("OOS: Upgraded", penalty_type), portfolio_name, freq_multiplier = freq_m)
  turnover <- round(calc_turnover(weights_mat), 4)
  
  return(list(
    Results = res_oos, 
    Turnover = turnover, 
    Weights = weights_mat,      # Added this
    OOS_Returns = oos_ret       # Added this
  ))
}

upgraded_lasso <- evaluate_upgraded_penalties(X100_monthly, "100 Size Monthly", "Lasso", "monthly")
upgraded_enet  <- evaluate_upgraded_penalties(X100_monthly, "100 Size Monthly", "ElasticNet", "monthly")

cat("\n--- FINAL MASTER SUMMARY REPORT (RAW RETURNS FRAMEWORK) ---\n")
final_comparison_table <- bind_rows(results_q3$Results, upgraded_lasso$Results, upgraded_enet$Results)
print(final_comparison_table)

cat("\nGenerating Q4 Visualizations...\n")
library(ggplot2)
library(tidyr)
library(dplyr)

cat("Re-evaluating L2 Baseline for 100 Size Monthly to extract raw returns...\n")
baseline_l2 <- evaluate_dynamic_l2_q23(X100_monthly, "100 Size Monthly", frequency = "monthly")

min_length <- min(length(baseline_l2$OOS_Returns), length(upgraded_lasso$OOS_Returns))

ret_l2 <- tail(baseline_l2$OOS_Returns, min_length) 
ret_lasso <- tail(upgraded_lasso$OOS_Returns, min_length)
ret_enet <- tail(upgraded_enet$OOS_Returns, min_length)


cum_l2 <- cumprod(1 + ret_l2)
cum_lasso <- cumprod(1 + ret_lasso)
cum_enet <- cumprod(1 + ret_enet)

dates_plot <- seq(1, min_length)

plot_data_ret <- data.frame(
  Time = rep(dates_plot, 3),
  Cumulative_Return = c(cum_l2, cum_lasso, cum_enet),
  Strategy = rep(c("L2 (Ridge)", "L1 (Lasso)", "Elastic Net"), each = min_length)
)

p_cum_ret <- ggplot(plot_data_ret, aes(x = Time, y = Cumulative_Return, color = Strategy)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("L2 (Ridge)" = "#7B1113", "L1 (Lasso)" = "#1B3A61", "Elastic Net" = "#E3B505")) +
  labs(
    title = "Out-of-Sample Cumulative Returns (100 Size Monthly)",
    x = "Months Out-of-Sample",
    y = "Cumulative Wealth"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(family = "serif")
  )

ggsave("Q4_Cumulative_Returns.png", plot = p_cum_ret, width = 10, height = 6, dpi = 300)

sparsity_l2 <- apply(baseline_l2$Weights_L2, 1, function(x) sum(abs(x) < 1e-4) / length(x))
sparsity_lasso <- apply(upgraded_lasso$Weights, 1, function(x) sum(abs(x) < 1e-4) / length(x))
sparsity_enet <- apply(upgraded_enet$Weights, 1, function(x) sum(abs(x) < 1e-4) / length(x))

min_windows <- min(length(sparsity_l2), length(sparsity_lasso))
sparsity_l2 <- tail(sparsity_l2, min_windows)
sparsity_lasso <- tail(sparsity_lasso, min_windows)
sparsity_enet <- tail(sparsity_enet, min_windows)

plot_data_spars <- data.frame(
  Window = rep(1:min_windows, 3),
  Sparsity = c(sparsity_l2, sparsity_lasso, sparsity_enet),
  Strategy = rep(c("L2 (Ridge)", "L1 (Lasso)", "Elastic Net"), each = min_windows)
)

p_sparsity <- ggplot(plot_data_spars, aes(x = Window, y = Sparsity, fill = Strategy)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("L2 (Ridge)" = "#7B1113", "L1 (Lasso)" = "#1B3A61", "Elastic Net" = "#E3B505")) +
  labs(
    title = "Portfolio Sparsity: Percentage of Weights Set to Zero",
    x = "Rebalancing Window",
    y = "% of Assets Excluded (Zero Weight)"
  ) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(family = "serif")
  )

ggsave("Q4_Sparsity_Barplot.png", plot = p_sparsity, width = 10, height = 6, dpi = 300)

cat("\nQ4 Plots successfully generated!\n")


library(MASS)
library(tidyverse)
library(patchwork)
library(stcov)

portfolio_perf <- function(w, oos_returns, freq, Sigma) {
  oos_ret <- oos_returns %*% w
  data.frame(
    Return     = mean(oos_ret) * freq,
    Volatility = sd(oos_ret) * sqrt(freq),
    Sharpe     = if(sd(oos_ret) > 0) (mean(oos_ret) * freq) / (sd(oos_ret) * sqrt(freq)) else NA
  )
}

turnover <- function(w_mat) {
  if(nrow(w_mat) <= 1) return(0)
  diffs <- diff(w_mat)
  mean(rowSums(abs(diffs)))
}

iso_eig <- function(lambda, T_est) {
  diag_mat <- diag(lambda)
  shrunk_cov <- stcov::iso_cov(diag_mat, T_est)
  return(eigen(shrunk_cov)$values)
}

out_of_sample_pca_portfolio <- function(df, K, freq = 12, window_years = 10, step_period = NULL, Stein = FALSE){
  N <- ncol(df) - 1 
  results <- list() 
  
  if (is.null(step_period)) {
    step_period <- ifelse(freq == 12, 6, 126) 
  }
  
  weights_mv <- matrix(nrow = 0, ncol = N)
  weights_msr <- matrix(nrow = 0, ncol = N)
  
  for (start in seq(1, nrow(df) - window_years * freq - step_period + 1, by = step_period)) {
    
    est_returns <- df[start:(start + window_years * freq - 1), -1] %>% as.matrix() 
    oos_returns <- df[(start + window_years * freq):(start + window_years * freq + step_period - 1), -1] %>% as.matrix() 
    
    mu <- colMeans(est_returns) 
    Sigma <- cov(est_returns) 
    
    eig <- eigen(Sigma) 
    V <- eig$vectors 
    lambda <- eig$values 
    V_k <- V[, 1:K, drop = FALSE] 
    
    if (Stein) { 
      T_est <- nrow(est_returns)
      lambda_k <- iso_eig(lambda, T_est)[1:K] 
    } else {
      lambda_k <- lambda[1:K] 
    }
    
    inv_lambda_k <- if (length(lambda_k) == 1) {
      matrix(1 / lambda_k, nrow = 1, ncol = 1)
    } else {
      diag(1 / lambda_k)
    }
    
    Sigma_inv <- V_k %*% inv_lambda_k %*% t(V_k) 
    
    w_mv <- rowSums(Sigma_inv); w_mv <- w_mv / sum(w_mv) 
    w_msr <- Sigma_inv %*% mu; w_msr <- w_msr / sum(w_msr) 
    
    weights_mv <- rbind(weights_mv, t(w_mv))
    weights_msr <- rbind(weights_msr, t(w_msr))
    
    results[[length(results) + 1]] <- rbind( 
      "Minimum Variance" = portfolio_perf(w_mv, oos_returns, freq, Sigma), 
      "Mean-Variance"    = portfolio_perf(w_msr, oos_returns, freq, Sigma) 
    )
  }
  
  perf_oos <- do.call(rbind, results) 
  perf_oos_df <- as.data.frame(perf_oos) 
  perf_oos_df$Strategy <- rep(c("Minimum Variance", "Mean-Variance"), length(results)) 
  
  perf_summary <- perf_oos_df %>%
    group_by(Strategy) %>%
    summarise(across(c(Return, Volatility, Sharpe), mean), .groups = "drop") 
  
  turnovers <- tibble(
    Strategy = c("Minimum Variance", "Mean-Variance"),
    Turnover = c(
      turnover(weights_mv), 
      turnover(weights_msr) 
    )
  )
  
  perf_summary <- left_join(perf_summary, turnovers, by = "Strategy") 
  return(perf_summary) 
}

plot_perf_vs_K <- function(df, freq = 12, window_years = 10, step_period = NULL, Stein = FALSE, dataset_label = "") {
  
  max_K <- ncol(df) - 1
  
  perf_data <- map_dfr(1:(max_K - 1), function(K) {
    perf <- out_of_sample_pca_portfolio(df, K, freq, window_years, step_period, Stein)
    perf$K <- K
    perf
  })
  
  invalid_combos <- perf_data %>%
    filter(Return > 1 | Return < -2) %>%
    dplyr::select(K, Strategy)
  
  perf_data_filtered <- perf_data %>%
    anti_join(invalid_combos, by = c("K", "Strategy"))
  
  perf_long <- perf_data_filtered %>%
    pivot_longer(cols = c(Return, Volatility, Sharpe, Turnover),
                 names_to = "Metric", values_to = "Value")
  
  combined_plot <- ggplot(perf_long, aes(x = K, y = Value, color = Strategy)) +
    geom_line(linewidth = 1.2) + 
    geom_point(size = 2, alpha = 0.9) + 
    facet_wrap(~ Metric, scales = "free_y", ncol = 2) + 
    
    scale_color_manual(values = c("Minimum Variance" = "#003366", 
                                  "Mean-Variance"    = "#990000")) +
    
    labs(
      title = paste(dataset_label, if(Stein) "(Stein Shrunk)" else "(Pure PCA)"),
      x = "Number of Principal Components (K)",
      y = NULL 
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
      
      strip.background = element_rect(fill = "#f8f9fa", color = "grey80"),
      strip.text = element_text(face = "bold", size = 12, color = "#333333"),
      
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.border = element_rect(color = "grey80"),
      
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.margin = margin(t = -10),
      
      axis.text = element_text(color = "black")
    )
  
  return(combined_plot)
}

cat("\n--- AUTOMATED DATA-DRIVEN K SELECTION FOR ALL 8 PORTFOLIOS (90% THRESHOLD) ---\n")

extract_optimal_K <- function(df, threshold = 0.90) {
  returns <- df %>% select_at(vars(-matches("date"))) %>% as.matrix()
  Sigma <- cov(returns)
  eig <- eigen(Sigma)
  explained_var <- cumsum(eig$values) / sum(eig$values)
  K_opt <- which(explained_var >= threshold)[1]
  
  return(tibble(
    Optimal_K = K_opt,
    Actual_Explained_Var = explained_var[K_opt]
  ))
}

all_datasets <- list(
  "10 Industry Daily"    = X10_daily,
  "10 Industry Monthly"  = X10_monthly,
  "25 Size-BM Daily"     = X25_daily,
  "25 Size-BM Monthly"   = X25_monthly,
  "48 Industry Daily"    = X48_daily,
  "48 Industry Monthly"  = X48_monthly,
  "100 Size-BM Daily"    = X100_daily,
  "100 Size-BM Monthly"  = X100_monthly
)

optimal_K_df <- imap_dfr(all_datasets, function(df, name) {
  metrics <- extract_optimal_K(df, threshold = 0.90)
  metrics$Dataset <- name
  metrics
}) %>% dplyr::select(Dataset, Optimal_K, Actual_Explained_Var)

print(as.data.frame(optimal_K_df))

get_k <- function(dataset_name) {
  optimal_K_df$Optimal_K[optimal_K_df$Dataset == dataset_name]
}

p_10_m  <- plot_perf_vs_K(X10_monthly, freq = 12, dataset_label = "10 Industry Monthly")
p_10_d  <- plot_perf_vs_K(X10_daily, freq = 252, dataset_label = "10 Industry Daily")
p_25_m  <- plot_perf_vs_K(X25_monthly, freq = 12, dataset_label = "25 Size-BM Monthly")
p_25_d  <- plot_perf_vs_K(X25_daily, freq = 252, dataset_label = "25 Size-BM Daily")
p_48_m  <- plot_perf_vs_K(X48_monthly, freq = 12, dataset_label = "48 Industry Monthly")
p_48_d  <- plot_perf_vs_K(X48_daily, freq = 252, dataset_label = "48 Industry Daily")
p_100_m <- plot_perf_vs_K(X100_monthly, freq = 12, dataset_label = "100 Size-BM Monthly")
p_100_d <- plot_perf_vs_K(X100_daily, freq = 252, dataset_label = "100 Size-BM Daily")

window_pure_10  <- p_10_d  + p_10_m
window_pure_25  <- p_25_d  + p_25_m
window_pure_48  <- p_48_d  + p_48_m
window_pure_100 <- p_100_d + p_100_m

print(window_pure_10)
print(window_pure_25)
print(window_pure_48)
print(window_pure_100) 

p_10_m_stein  <- plot_perf_vs_K(X10_monthly, freq = 12, Stein = TRUE, dataset_label = "10 Industry Monthly")
p_10_d_stein  <- plot_perf_vs_K(X10_daily, freq = 252, Stein = TRUE, dataset_label = "10 Industry Daily")
p_25_m_stein  <- plot_perf_vs_K(X25_monthly, freq = 12, Stein = TRUE, dataset_label = "25 Size-BM Monthly")
p_25_d_stein  <- plot_perf_vs_K(X25_daily, freq = 252, Stein = TRUE, dataset_label = "25 Size-BM Daily")
p_48_m_stein  <- plot_perf_vs_K(X48_monthly, freq = 12, Stein = TRUE, dataset_label = "48 Industry Monthly")
p_48_d_stein  <- plot_perf_vs_K(X48_daily, freq = 252, Stein = TRUE, dataset_label = "48 Industry Daily")
p_100_m_stein <- plot_perf_vs_K(X100_monthly, freq = 12, Stein = TRUE, dataset_label = "100 Size-BM Monthly")
p_100_d_stein <- plot_perf_vs_K(X100_daily, freq = 252, Stein = TRUE, dataset_label = "100 Size-BM Daily")

window_stein_10  <- p_10_d_stein  + p_10_m_stein
window_stein_25  <- p_25_d_stein  + p_25_m_stein
window_stein_48  <- p_48_d_stein  + p_48_m_stein
window_stein_100 <- p_100_d_stein + p_100_m_stein

print(window_stein_10)
print(window_stein_25)
print(window_stein_48)
print(window_stein_100)

summary_10_m  <- out_of_sample_pca_portfolio(X10_monthly,  K = get_k("10 Industry Monthly"),   freq = 12)
summary_10_d  <- out_of_sample_pca_portfolio(X10_daily,    K = get_k("10 Industry Daily"),     freq = 252)
summary_25_m  <- out_of_sample_pca_portfolio(X25_monthly,  K = get_k("25 Size-BM Monthly"),    freq = 12)
summary_25_d  <- out_of_sample_pca_portfolio(X25_daily,    K = get_k("25 Size-BM Daily"),      freq = 252)
summary_48_m  <- out_of_sample_pca_portfolio(X48_monthly,  K = get_k("48 Industry Monthly"),   freq = 12)
summary_48_d  <- out_of_sample_pca_portfolio(X48_daily,    K = get_k("48 Industry Daily"),     freq = 252)
summary_100_m <- out_of_sample_pca_portfolio(X100_monthly, K = get_k("100 Size-BM Monthly"),   freq = 12)
summary_100_d <- out_of_sample_pca_portfolio(X100_daily,   K = get_k("100 Size-BM Daily"),     freq = 252)

summary_pure_pca_all <- bind_rows(
  summary_10_m  %>% mutate(Dataset = "10 Industry Monthly"),
  summary_10_d  %>% mutate(Dataset = "10 Industry Daily"),
  summary_25_m  %>% mutate(Dataset = "25 Size-BM Monthly"),
  summary_25_d  %>% mutate(Dataset = "25 Size-BM Daily"),
  summary_48_m  %>% mutate(Dataset = "48 Industry Monthly"),
  summary_48_d  %>% mutate(Dataset = "48 Industry Daily"),
  summary_100_m %>% mutate(Dataset = "100 Size-BM Monthly"),
  summary_100_d %>% mutate(Dataset = "100 Size-BM Daily")
) %>% 
  left_join(optimal_K_df, by = "Dataset")

stein_10_m  <- out_of_sample_pca_portfolio(X10_monthly,  K = get_k("10 Industry Monthly"),   freq = 12,  Stein = TRUE)
stein_10_d  <- out_of_sample_pca_portfolio(X10_daily,    K = get_k("10 Industry Daily"),     freq = 252, Stein = TRUE)
stein_25_m  <- out_of_sample_pca_portfolio(X25_monthly,  K = get_k("25 Size-BM Monthly"),    freq = 12,  Stein = TRUE)
stein_25_d  <- out_of_sample_pca_portfolio(X25_daily,    K = get_k("25 Size-BM Daily"),      freq = 252, Stein = TRUE)
stein_48_m  <- out_of_sample_pca_portfolio(X48_monthly,  K = get_k("48 Industry Monthly"),   freq = 12,  Stein = TRUE)
stein_48_d  <- out_of_sample_pca_portfolio(X48_daily,    K = get_k("48 Industry Daily"),     freq = 252, Stein = TRUE)
stein_100_m <- out_of_sample_pca_portfolio(X100_monthly, K = get_k("100 Size-BM Monthly"),   freq = 12,  Stein = TRUE)
stein_100_d <- out_of_sample_pca_portfolio(X100_daily,   K = get_k("100 Size-BM Daily"),     freq = 252, Stein = TRUE)

summary_stein_all <- bind_rows(
  stein_10_m  %>% mutate(Dataset = "10 Industry Monthly"),
  stein_10_d  %>% mutate(Dataset = "10 Industry Daily"),
  stein_25_m  %>% mutate(Dataset = "25 Size-BM Monthly"),
  stein_25_d  %>% mutate(Dataset = "25 Size-BM Daily"),
  stein_48_m  %>% mutate(Dataset = "48 Industry Monthly"),
  stein_48_d  %>% mutate(Dataset = "48 Industry Daily"),
  stein_100_m %>% mutate(Dataset = "100 Size-BM Monthly"),
  stein_100_d %>% mutate(Dataset = "100 Size-BM Daily")
) %>% 
  left_join(optimal_K_df, by = "Dataset")

print(as.data.frame(summary_pure_pca_all))
print(as.data.frame(summary_stein_all))

cat("\n--- EXPORTING PLOTS TO 'Q5 plots' DIRECTORY ---\n")

export_dir <- "Q5 plots"
if (!dir.exists(export_dir)) {
  dir.create(export_dir)
}

plot_width  <- 14
plot_height <- 7
plot_dpi    <- 300

ggsave(filename = file.path(export_dir, "Pure_PCA_10_Industry.png"),  plot = window_pure_10,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Pure_PCA_25_SizeBM.png"),    plot = window_pure_25,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Pure_PCA_48_Industry.png"),  plot = window_pure_48,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Pure_PCA_100_SizeBM.png"),   plot = window_pure_100, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

ggsave(filename = file.path(export_dir, "Stein_10_Industry.png"),  plot = window_stein_10,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Stein_25_SizeBM.png"),    plot = window_stein_25,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Stein_48_Industry.png"),  plot = window_stein_48,  width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
ggsave(filename = file.path(export_dir, "Stein_100_SizeBM.png"),   plot = window_stein_100, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

cat("Export complete! Check your working directory.\n")

library(fastICA)
library(PerformanceAnalytics)
library(tidyverse)
library(MASS)

calc_weights_icvp <- function(train_data, K) {
  N <- ncol(train_data)
  
  if (K == 1) {
    sigma <- stats::cov(train_data)
    eig <- eigen(sigma)
    A <- eig$vectors[, 1, drop=FALSE] * (eig$values[1]^(-0.5))
    sign_K <- sign(sum(A))
    w <- A * sign_K
    return(as.numeric(w / sum(w)))
  }
  
  set.seed(12345) 
  res <- tryCatch({
    fastICA(train_data, n.comp = K, alg.typ = "parallel", fun = "logcosh", 
            method = "C", row.norm = FALSE, maxit = 500, tol = 1e-04)
  }, error = function(e) NULL)
  
  if (is.null(res)) {
    sigma <- stats::cov(train_data)
    eig <- eigen(sigma)
    V_k <- eig$vectors[, 1:K, drop = FALSE]
    Lambda_k <- diag(eig$values[1:K]^(-0.5), nrow = K, ncol = K)
    A <- V_k %*% Lambda_k
  } else {
    A <- res$K %*% res$W 
  }
  
  sign_K <- sign(colSums(A))
  w_unnorm <- A %*% sign_K
  
  return(as.numeric(w_unnorm / sum(w_unnorm)))
}

calc_weights_minvar <- function(train_data) {
  sigma <- stats::cov(train_data)
  inv_cov <- ginv(sigma)
  ones <- rep(1, ncol(train_data))
  num <- inv_cov %*% ones
  den <- as.numeric(t(ones) %*% inv_cov %*% ones)
  return(as.numeric(num / den))
}

calc_higher_moments <- function(returns_vector, name) {
  ann_mean <- mean(returns_vector) * 12
  ann_vol  <- sd(returns_vector) * sqrt(12)
  
  skew     <- PerformanceAnalytics::skewness(returns_vector, method = "sample")
  kurt     <- PerformanceAnalytics::kurtosis(returns_vector, method = "excess")
  
  mVaR_95  <- as.numeric(PerformanceAnalytics::VaR(returns_vector, p = 0.95, method = "modified"))
  
  data.frame(
    Strategy = name,
    Ann_Mean = round(ann_mean, 4),
    Ann_Vol  = round(ann_vol, 4),
    Skewness = round(skew, 4),
    Exc_Kurtosis = round(kurt, 4),
    Mod_VaR_95 = round(mVaR_95, 4)
  )
}

evaluate_q6_strategies <- function(raw_returns_df, K_value = NULL, is_baseline = FALSE) {
  dates <- as.Date(paste0(raw_returns_df$date, "01"), format = "%Y%m%d")
  ret_data <- as.matrix(raw_returns_df[, -1])
  
  outer_windows <- make_calendar_rolling_windows(dates, train_years = 10, test_months = 6, step_months = 6)
  oos_ret <- c()
  
  for (i in seq_len(nrow(outer_windows))) {
    train_idx <- which(dates >= outer_windows$train_start[i] & dates <= outer_windows$train_end[i])
    test_idx  <- which(dates >= outer_windows$test_start[i]  & dates <= outer_windows$test_end[i])
    
    if(length(train_idx) == 0 || length(test_idx) == 0) next
    train_data <- ret_data[train_idx, , drop = FALSE]
    test_data  <- ret_data[test_idx, , drop = FALSE]
    
    if (is_baseline) {
      w_oos <- calc_weights_minvar(train_data)
    } else {
      w_oos <- calc_weights_icvp(train_data, K_value)
    }
    
    oos_ret <- c(oos_ret, test_data %*% w_oos)
  }
  
  strat_name <- if(is_baseline) "Sample MinVar" else paste("IC-Variance-Parity K =", K_value)
  return(calc_higher_moments(oos_ret, strat_name))
}

k_grid <- c(1, 5, 10, 25, 50, 75, 100)
results_list <- list()

cat("\nEvaluating IC-Variance-Parity across K values...\n")
for (k in k_grid) {
  results_list[[length(results_list) + 1]] <- evaluate_q6_strategies(X100_monthly, K_value = k)
}

cat("Evaluating Sample Min-Var Baseline...\n")
results_list[[length(results_list) + 1]] <- evaluate_q6_strategies(X100_monthly, is_baseline = TRUE)

q6_master_table <- do.call(rbind, results_list)

cat("\n--- QUESTION 6: HIGHER MOMENTS & mVaR MASTER TABLE ---\n")
print(q6_master_table)

evaluate_icmv_components <- function(raw_returns_df, K_value) {
  dates <- as.Date(paste0(raw_returns_df$date, "01"), format = "%Y%m%d")
  ret_data <- as.matrix(raw_returns_df[, -1])
  asset_names <- colnames(ret_data)
  
  outer_windows <- make_calendar_rolling_windows(dates, train_years = 10, test_months = 6, step_months = 6)

  mv_weights_history <- list()
  ic_weights_history <- list()
  
  for (i in seq_len(nrow(outer_windows))) {
    train_idx <- which(dates >= outer_windows$train_start[i] & dates <= outer_windows$train_end[i])
    if(length(train_idx) == 0) next
    
    train_data <- ret_data[train_idx, , drop = FALSE]
    
    w_mv <- calc_weights_minvar(train_data)
    
    w_ic <- calc_weights_icvp(train_data, K_value)

    mv_weights_history[[i]] <- data.frame(
      Window = i, Date = outer_windows$test_start[i], Strategy = "w_MV", 
      Asset = asset_names, Weight = w_mv
    )
    ic_weights_history[[i]] <- data.frame(
      Window = i, Date = outer_windows$test_start[i], Strategy = "w_IC", 
      Asset = asset_names, Weight = w_ic
    )
  }
  
  df_w_mv <- bind_rows(mv_weights_history)
  df_w_ic <- bind_rows(ic_weights_history)
  
  return(list(w_MV_history = df_w_mv, w_IC_history = df_w_ic))
}


# ==============================================================================
# STEP 6: ICMV SHRINKAGE OPTIMIZATION ENGINE & BACKTEST
# ==============================================================================
cat("\nRunning In-Sample delta calibration and Out-of-Sample ICMV evaluation...\n")

# --- Helper Function: Calculate Portfolio Higher Moments & MVaR ---
calc_portfolio_mvar <- function(p_returns, alpha = 0.01) {
  mu_p    <- mean(p_returns)
  sigma_p <- sd(p_returns)
  
  # Extract sample skewness and excess kurtosis safely
  skew_p <- PerformanceAnalytics::skewness(p_returns, method = "sample")
  kurt_p <- PerformanceAnalytics::kurtosis(p_returns, method = "excess")
  
  # Gaussian quantile (e.g., -2.326348 for alpha = 1%)
  z_alpha <- qnorm(alpha)
  
  # Cornish-Fisher parameter adjustment terms
  term1 <- z_alpha
  term2 <- (1/6) * (z_alpha^2 - 1) * skew_p
  term3 <- (1/24) * (z_alpha^3 - 3 * z_alpha) * kurt_p
  term4 <- - (1/36) * (2 * z_alpha^3 - 5 * z_alpha) * (skew_p^2)
  
  # Compute Modified VaR
  mvar <- - (mu_p + sigma_p * (term1 + term2 + term3 + term4))
  
  return(list(mu = mu_p, mvar = mvar))
}

# --- Main Out-of-Sample Evaluation Engine ---
evaluate_icmv_portfolio <- function(raw_returns_df, K_value) {
  dates <- as.Date(paste0(raw_returns_df$date, "01"), format = "%Y%m%d")
  ret_data <- as.matrix(raw_returns_df[, -1])
  
  # Setup rolling calendar windows (matching your prior configurations)
  outer_windows <- make_calendar_rolling_windows(dates, train_years = 10, test_months = 6, step_months = 6)
  
  oos_icmv_returns <- c()
  calibrated_deltas <- c()
  
  # Define optimization grid search resolution for delta
  delta_grid <- seq(0, 1, by = 0.05)
  
  for (i in seq_len(nrow(outer_windows))) {
    train_idx <- which(dates >= outer_windows$train_start[i] & dates <= outer_windows$train_end[i])
    test_idx  <- which(dates >= outer_windows$test_start[i]  & dates <= outer_windows$test_end[i])
    
    if(length(train_idx) == 0 || length(test_idx) == 0) next
    
    train_data <- ret_data[train_idx, , drop = FALSE]
    test_data  <- ret_data[test_idx, , drop = FALSE]
    
    # 1. Compute baseline weight components for this window
    w_mv <- calc_weights_minvar(train_data)
    w_ic <- calc_weights_icvp(train_data, K_value)
    
    # 2. IN-SAMPLE CALIBRATION: Find delta that maximizes Modified Sharpe Ratio
    best_delta <- 0
    max_mod_sharpe <- -Inf
    
    for (d in delta_grid) {
      # Combine weights for in-sample testing
      w_blend_is <- (1 - d) * w_mv + d * w_ic
      is_portfolio_returns <- train_data %*% w_blend_is
      
      # Evaluate Cornish-Fisher metric parameters
      metrics <- calc_portfolio_mvar(is_portfolio_returns, alpha = 0.01)
      
      # Guard against negative or zero VaR to prevent division errors
      if (metrics$mvar > 0) {
        mod_sharpe <- metrics$mu / metrics$mvar
      } else {
        mod_sharpe <- -Inf
      }
      
      if (mod_sharpe > max_mod_sharpe) {
        max_mod_sharpe <- mod_sharpe
        best_delta <- d
      }
    }
    
    calibrated_deltas <- c(calibrated_deltas, best_delta)
    
    # 3. OUT-OF-SAMPLE IMPLEMENTATION: Apply optimized delta to unseen forward test data
    w_icmv_oos <- (1 - best_delta) * w_mv + best_delta * w_ic
    
    # Capitalize returns over the forward 6-month testing block
    window_oos_returns <- test_data %*% w_icmv_oos
    oos_icmv_returns <- c(oos_icmv_returns, window_oos_returns)
  }
  
  return(list(returns = oos_icmv_returns, deltas = calibrated_deltas))
}

# ==============================================================================
# STEP 7: EXECUTION AND BENCHMARK SUMMARY COMPARISON
# ==============================================================================

# Run the full optimization using K = 25 components for example
icmv_results_10_daily <- evaluate_icmv_portfolio(X10_daily, K_value = 5)
icmv_results_25_daily <- evaluate_icmv_portfolio(X25_daily, K_value = 10)
icmv_results_48_daily <- evaluate_icmv_portfolio(X48_daily, K_value = 20)
icmv_results_100_daily <- evaluate_icmv_portfolio(X100_daily, K_value = 25)
#monthly
icmv_results_10_monthly <- evaluate_icmv_portfolio(X10_monthly, K_value = 5)
icmv_results_25_monthly <- evaluate_icmv_portfolio(X25_monthly, K_value = 10)
icmv_results_48_monthly<- evaluate_icmv_portfolio(X48_monthly, K_value = 20)
icmv_results_100_monthly <- evaluate_icmv_portfolio(X100_monthly, K_value = 25)
# Summarize metrics using your historical moment engine
cat("\n--- FINAL EVALUATION METRICS FOR ICMV PORTFOLIO ---\n")


icmv_results <- list(
  "10_daily"    = icmv_results_10_daily$deltas,
  "25_daily"    = icmv_results_25_daily$deltas,
  "48_daily"    = icmv_results_48_daily$deltas,
  "100_daily"   = icmv_results_100_daily$deltas,
  "10_monthly"  = icmv_results_10_monthly$deltas, 
  "25_monthly"  = icmv_results_25_monthly$deltas,
  "48_monthly"  = icmv_results_48_monthly$deltas,
  "100_monthly" = icmv_results_100_monthly$deltas
)

for (name in names(icmv_results)) {
  cat("\n--- Results for:", name, "---\n")
  print(summary(icmv_results[[name]]))
}


icmv_metrics_10_daily <- calc_higher_moments(icmv_results_10_daily$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_25_daily <- calc_higher_moments(icmv_results_25_daily$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_48_daily <- calc_higher_moments(icmv_results_48_daily$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_100_daily <- calc_higher_moments(icmv_results_100_daily$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_10_mothly <- calc_higher_moments(icmv_results_10_monthly$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_25_monthly <- calc_higher_moments(icmv_results_25_monthly$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_48_monthly <- calc_higher_moments(icmv_results_48_monthly$returns, "ICMV Shrinkage Portfolio")
icmv_metrics_100_monthly <- calc_higher_moments(icmv_results_100_monthly$returns, "ICMV Shrinkage Portfolio")

icmv_metrics_all <- list(
  "10_daily"    = icmv_results_10_daily$returns,
  "25_daily"    = icmv_results_25_daily$returns,
  "48_daily"    = icmv_results_48_daily$returns,
  "100_daily"   = icmv_results_100_daily$returns,
  "10_monthly"  = icmv_results_10_monthly$returns, 
  "25_monthly"  = icmv_results_25_monthly$returns,
  "48_monthly"  = icmv_results_48_monthly$returns,
  "100_monthly" = icmv_results_100_monthly$returns
)

for (name in names(icmv_metrics_all)) {
  cat("\n--- Metrics for:", name, "---\n")
  
  metrics <- calc_higher_moments(icmv_metrics_all[[name]], "ICMV Shrinkage Portfolio")
  print(metrics)
}


# ==============================================================================
# STEP 8: VISUALIZING THE 8 EFFICIENT FRONTIERS WITH ICMV LOCATIONS
# ==============================================================================
cat("\nGenerating all 8 Efficient Frontier plots comparing GMV and ICMV...\n")

library(ggplot2)
library(quadprog)
library(dplyr)

# --- Configuration Mapping List ---
# This matches your exact variable names and sets correct parameters dynamically
frontier_config <- list(
  list(df = X10_daily,    freq = "daily",   name = "10_daily",    label = "10 Industry (Daily)",   K = 5,   res_obj = icmv_results_10_daily),
  list(df = X25_daily,    freq = "daily",   name = "25_daily",    label = "25 Portfolios (Daily)", K = 10,  res_obj = icmv_results_25_daily),
  list(df = X48_daily,    freq = "daily",   name = "48_daily",    label = "48 Industry (Daily)",   K = 20,  res_obj = icmv_results_48_daily),
  list(df = X100_daily,   freq = "daily",   name = "100_daily",   label = "100 Portfolios (Daily)",K = 25,  res_obj = icmv_results_100_daily),
  
  list(df = X10_monthly,  freq = "monthly", name = "10_monthly",  label = "10 Industry (Monthly)", K = 5,   res_obj = icmv_results_10_monthly),
  list(df = X25_monthly,  freq = "monthly", name = "25_monthly",  label = "25 Portfolios (Monthly)",K = 10, res_obj = icmv_results_25_monthly),
  list(df = X48_monthly,  freq = "monthly", name = "48_monthly",  label = "48 Industry (Monthly)", K = 20,  res_obj = icmv_results_48_monthly),
  list(df = X100_monthly, freq = "monthly", name = "100_monthly", label = "100 Portfolios (Monthly)",K = 25,res_obj = icmv_results_100_monthly)
)

# Global holder list in the environment to access plots later
all_frontier_plots <- list()

for (track in frontier_config) {
  cat(paste("Computing frontier coordinates for:", track$label, "...\n"))
  
  # 1. Establish Annualization Scaling Factors
  freq_mult <- if (track$freq == "daily") 252 else 12
  
  ret_data <- as.matrix(track$df[, -1])
  N <- ncol(ret_data)
  
  mu    <- colMeans(ret_data) * freq_mult
  sigma <- cov(ret_data) * freq_mult
  
  # Regularize if near-singular fallback safeguard to prevent quadprog crashes
  if (rcond(sigma) < 1e-12) {
    sigma <- sigma + 1e-6 * diag(diag(sigma))
  }
  
  # 2. Map the Hyperbola Frontier Line
  ones <- rep(1, N)
  inv_sigma <- solve(sigma)
  
  target_returns <- seq(min(mu), max(mu) * 1.1, length.out = 100)
  frontier_risks <- numeric(100)
  Amat <- cbind(ones, mu)
  
  for (i in seq_along(target_returns)) {
    bvec <- c(1, target_returns[i])
    qp_sol <- tryCatch({
      solve.QP(Dmat = sigma, dvec = rep(0, N), Amat = Amat, bvec = bvec, meq = 2)
    }, error = function(e) NULL)
    
    frontier_risks[i] <- if (!is.null(qp_sol)) sqrt(2 * qp_sol$value) else NA
  }
  
  df_frontier <- data.frame(Risk = frontier_risks, Return = target_returns) %>% 
    dplyr::filter(!is.na(Risk))
  
  # 3. Calculate Global Min-Variance (GMV) Coordinates
  w_gmv <- as.numeric((inv_sigma %*% ones) / as.numeric(t(ones) %*% inv_sigma %*% ones))
  
  # 4. Calculate Steady-State ICMV Coordinates using your generated list object deltas
  avg_delta <- mean(track$res_obj$deltas)
  w_ic_full <- calc_weights_icvp(ret_data, K = track$K)
  w_icmv_fixed <- (1 - avg_delta) * w_gmv + avg_delta * w_ic_full
  
  df_points <- data.frame(
    Risk = c(
      sqrt(as.numeric(t(w_gmv) %*% sigma %*% w_gmv)), 
      sqrt(as.numeric(t(w_icmv_fixed) %*% sigma %*% w_icmv_fixed))
    ),
    Return = c(
      as.numeric(t(w_gmv) %*% mu), 
      as.numeric(t(w_icmv_fixed) %*% mu)
    ),
    Label = c("Global Min-Variance (GMV)", "ICMV Shrinkage Portfolio")
  )
  
  df_assets <- data.frame(Risk = sqrt(diag(sigma)), Return = mu)
  
  # 5. Build Individual Chart Frame
  p <- ggplot() +
    geom_path(data = df_frontier, aes(x = Risk, y = Return), color = "#34495E", linewidth = 1) +
    geom_point(data = df_assets, aes(x = Risk, y = Return), color = "#BDC3C7", alpha = 0.5, size = 1.2) +
    geom_point(data = df_points, aes(x = Risk, y = Return, fill = Label, shape = Label), color = "black", size = 4) +
    scale_shape_manual(values = c("Global Min-Variance (GMV)" = 23, "ICMV Shrinkage Portfolio" = 24)) +
    scale_fill_manual(values = c("Global Min-Variance (GMV)" = "#E74C3C", "ICMV Shrinkage Portfolio" = "#F1C40F")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title = paste("Efficient Frontier:", track$label),
      subtitle = paste("Avg. Calibrated Intensity \u03b4 =", round(avg_delta, 3)),
      x = "Annualized Volatility (Standard Deviation)",
      y = "Annualized Expected Return",
      fill = "Portfolios", shape = "Portfolios"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "#7F8C8D", hjust = 0.5),
      legend.position = "bottom"
    )
  
  # Print directly to RStudio screen plots layout view as it runs
  print(p)
  
  # Save to list framework inside your environment memory
  all_frontier_plots[[track$name]] <- p
}

cat("\nAll 8 frontiers plotted successfully and stored in environment list object 'all_frontier_plots'!\n")



# ==============================================================================
# STEP 9: COMPILE ALL 8 ICMV COORDINATES & DELTAS INTO A UNIFIED TABLE
# ==============================================================================
cat("\nCompiling all 8 ICMV allocations, risk profiles, and delta factors...\n")

library(dplyr)
library(purrr)

# --- Define Mapping Matrix for Full Extraction ---
icmv_table_config <- list(
  list(df = X10_daily,    freq = "Daily",   name = "10 Industry",    K = 5,   res_obj = icmv_results_10_daily),
  list(df = X25_daily,    freq = "Daily",   name = "25 Portfolios",  K = 10,  res_obj = icmv_results_25_daily),
  list(df = X48_daily,    freq = "Daily",   name = "48 Industry",    K = 20,  res_obj = icmv_results_48_daily),
  list(df = X100_daily,   freq = "Daily",   name = "100 Portfolios", K = 25,  res_obj = icmv_results_100_daily),
  
  list(df = X10_monthly,  freq = "Monthly", name = "10 Industry",    K = 5,   res_obj = icmv_results_10_monthly),
  list(df = X25_monthly,  freq = "Monthly", name = "25 Portfolios",  K = 10,  res_obj = icmv_results_25_monthly),
  list(df = X48_monthly,  freq = "Monthly", name = "48 Industry",    K = 20,  res_obj = icmv_results_48_monthly),
  list(df = X100_monthly, freq = "Monthly", name = "100 Portfolios", K = 25,  res_obj = icmv_results_100_monthly)
)

# --- Process and Map Elements directly into an Active Tibble ---
icmv_summary_coordinates_table <- map_df(icmv_table_config, function(track) {
  freq_mult <- if (track$freq == "Daily") 252 else 12
  
  ret_data <- as.matrix(track$df[, -1])
  N        <- ncol(ret_data)
  
  mu    <- colMeans(ret_data) * freq_mult
  sigma <- cov(ret_data) * freq_mult
  
  # Regularize if near-singular fallback safeguard to match plotting engines
  if (rcond(sigma) < 1e-12) {
    sigma <- sigma + 1e-6 * diag(diag(sigma))
  }
  
  ones <- rep(1, N)
  inv_sigma <- solve(sigma)
  
  # Calculate Baseline GMV and IC parameters
  w_gmv     <- as.numeric((inv_sigma %*% ones) / as.numeric(t(ones) %*% inv_sigma %*% ones))
  w_ic_full <- calc_weights_icvp(ret_data, K = track$K)
  
  # Extract long-term average shrinkage factor
  avg_delta    <- mean(track$res_obj$deltas)
  w_icmv_fixed <- (1 - avg_delta) * w_gmv + avg_delta * w_ic_full
  
  # Annualized coordinates derivations
  gmv_risk   <- sqrt(as.numeric(t(w_gmv) %*% sigma %*% w_gmv))
  gmv_return <- as.numeric(t(w_gmv) %*% mu)
  
  icmv_risk   <- sqrt(as.numeric(t(w_icmv_fixed) %*% sigma %*% w_icmv_fixed))
  icmv_return <- as.numeric(t(w_icmv_fixed) %*% mu)
  
  # Return row elements
  data.frame(
    Universe        = track$name,
    Frequency       = track$freq,
    K_Components    = track$K,
    Avg_Delta       = round(avg_delta, 4),
    GMV_Ann_Risk    = round(gmv_risk, 4),
    GMV_Ann_Return  = round(gmv_return, 4),
    ICMV_Ann_Risk   = round(icmv_risk, 4),
    ICMV_Ann_Return = round(icmv_return, 4)
  )
})

# --- Display Content directly to R Workspace Session ---
print(icmv_summary_coordinates_table)
cat("\nSummary frame stored in environment variable 'icmv_summary_coordinates_table'!\n")


# ==============================================================================
# STEP 10: OOS PERFORMANCE COMPARISON MATRIX & CUMULATIVE WEALTH PLOTS (FIXED)
# ==============================================================================
cat("\nExecuting multi-strategy OOS backtests (GMV vs IC vs ICMV)...\n")

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(PerformanceAnalytics)

calc_oos_table_metrics <- function(returns_vector, strategy_label, portfolio_label, freq_mult) {
  ann_mean <- mean(returns_vector) * freq_mult
  ann_vol  <- sd(returns_vector) * sqrt(freq_mult)
  sharpe   <- if(ann_vol > 0) ann_mean / ann_vol else NA
  
  data.frame(
    Portfolio = portfolio_label,
    Strategy  = strategy_label,
    Mean      = round(ann_mean, 4),
    Vol       = round(ann_vol, 4),
    Sharpe    = round(sharpe, 4)
  )
}

evaluate_triple_oos_strategies <- function(data_df, frequency, universe_name, K_value) {
  cat(paste("\nSimulating forward OOS loops for:", universe_name, "(", tools::toTitleCase(frequency), ") ...\n"))
  
  if (frequency == "monthly") {
    freq_mult <- 12; window_size <- 120; roll_step <- 6
    dates <- as.Date(paste0(data_df$date, "01"), format = "%Y%m%d")
  } else {
    freq_mult <- 252; window_size <- 2520; roll_step <- 126
    dates <- as.Date(as.character(data_df$date), format = "%Y%m%d")
  }
  
  ret_data <- as.matrix(data_df[, -1])
  T_total  <- nrow(ret_data)
  
  oos_dates   <- c()
  oos_ret_mv  <- c()
  oos_ret_ic  <- c()
  oos_ret_icmv <- c()
  
  delta_grid <- seq(0, 1, by = 0.05)
  
  for (t_end in seq(window_size, T_total - 1, by = roll_step)) {
    t_start    <- t_end - window_size + 1
    train_data <- ret_data[t_start:t_end, , drop = FALSE]
    
    test_start <- t_end + 1
    test_end   <- min(t_end + roll_step, T_total)
    test_data  <- ret_data[test_start:test_end, , drop = FALSE]
    
    if(length(train_data) == 0 || length(test_data) == 0) next
    
    w_mv <- calc_weights_minvar(train_data)
    w_ic <- calc_weights_icvp(train_data, K_value)
    
    best_delta <- 0
    max_mod_sharpe <- -Inf
    
    for (d in delta_grid) {
      w_blend_is <- (1 - d) * w_mv + d * w_ic
      is_returns <- train_data %*% w_blend_is
      
      mu_p   <- mean(is_returns)
      skew_p <- PerformanceAnalytics::skewness(is_returns, method = "sample")
      kurt_p <- PerformanceAnalytics::kurtosis(is_returns, method = "excess")
      
      z_alpha <- qnorm(0.01)
      mvar_is <- - (mu_p + sd(is_returns) * (z_alpha + (1/6)*(z_alpha^2-1)*skew_p + 
                                               (1/24)*(z_alpha^3-3*z_alpha)*kurt_p - (1/36)*(2*z_alpha^3-5*z_alpha)*(skew_p^2)))
      
      mod_sharpe <- if (mvar_is > 0) mu_p / mvar_is else -Inf
      
      if (mod_sharpe > max_mod_sharpe) {
        max_mod_sharpe <- mod_sharpe
        best_delta <- d
      }
    }
    
    w_icmv <- (1 - best_delta) * w_mv + best_delta * w_ic
    
    oos_ret_mv   <- c(oos_ret_mv,   test_data %*% w_mv)
    oos_ret_ic   <- c(oos_ret_ic,   test_data %*% w_ic)
    oos_ret_icmv <- c(oos_ret_icmv, test_data %*% w_icmv)
    
    oos_dates <- c(oos_dates, dates[test_start:test_end])
  }
  
  class(oos_dates) <- "Date"
  portfolio_label <- paste(universe_name, tools::toTitleCase(frequency), sep = " ")
  
  matrix_table <- bind_rows(
    calc_oos_table_metrics(oos_ret_mv,   "OOS: Min-Variance (GMV)", portfolio_label, freq_mult),
    calc_oos_table_metrics(oos_ret_ic,   "OOS: IC-Variance-Parity", portfolio_label, freq_mult),
    calc_oos_table_metrics(oos_ret_icmv, "OOS: ICMV Shrinkage",      portfolio_label, freq_mult)
  )
  
  print(matrix_table)
  
  # FIXED HERE: Clean column names to ensure stable pivoting behavior
  df_wealth <- data.frame(
    Date          = oos_dates,
    Min_Variance  = cumprod(1 + oos_ret_mv),
    IC_Parity     = cumprod(1 + oos_ret_ic),
    ICMV_Shrink   = cumprod(1 + oos_ret_icmv)
  ) %>%
    pivot_longer(-Date, names_to = "Strategy", values_to = "Wealth")
  
  # Aesthetic dictionary configurations
  strategy_colors <- c("Min_Variance" = "#619CFF", "IC_Parity" = "#F8766D", "ICMV_Shrink" = "#F1C40F")
  pretty_labels   <- c("Min_Variance" = "Min-Variance (GMV)", "IC_Parity" = "IC-Variance-Parity", "ICMV_Shrink" = "ICMV Shrinkage")
  
  p <- ggplot(df_wealth, aes(x = Date, y = Wealth, color = Strategy)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = strategy_colors, labels = pretty_labels) + # Clean values matching with labels mapping
    scale_y_continuous(labels = scales::dollar_format()) +
    labs(
      title = paste("OOS Cumulative Growth comparison: ", universe_name, " (", tools::toTitleCase(frequency), ")", sep=""),
      x = "Year", y = "Growth of $1 Investment", color = "Strategy Framework"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "bottom")
  
  print(p)
  
  return(list(Table = matrix_table, Plot = p))
}

# ==============================================================================
# ENVIRONMENT MULTI-TRACK EXECUTION CONTROL
# ==============================================================================
master_backtest_configs <- list(
  list(df = X10_monthly,  freq = "monthly", name = "10 Ind",   K = 5),
  list(df = X10_daily,    freq = "daily",   name = "10 Ind",   K = 5),
  list(df = X25_monthly,  freq = "monthly", name = "25 Port",  K = 10),
  list(df = X25_daily,    freq = "daily",   name = "25 Port",  K = 10),
  list(df = X48_monthly,  freq = "monthly", name = "48 Ind",   K = 20),
  list(df = X48_daily,    freq = "daily",   name = "48 Ind",   K = 20),
  list(df = X100_monthly, freq = "monthly", name = "100 Port", K = 25),
  list(df = X100_daily,   freq = "daily",   name = "100 Port", K = 25)
)

# Active workspace cache list object
oos_triple_results <- list()

for (run in master_backtest_configs) {
  storage_key <- paste(gsub(" ", "_", run$name), run$freq, sep = "_")
  
  oos_triple_results[[storage_key]] <- evaluate_triple_oos_strategies(
    data_df = run$df, 
    frequency = run$freq, 
    universe_name = run$name, 
    K_value = run$K
  )
}

cat("\nAll OOS matrices compiled and plots rendered into active memory framework objects!\n")







