# Author: Jean de Biolley #


# The following code backtests an industry momentum strategy by ranking sectors on 12-month rolling returns (Q3.2.1),
# evaluates the monthly and annualized Sharpe ratios of Winner and Loser portfolios against the market benchmark (Q3.2.2–3),
# plots cumulative wealth curves on a log scale (Q3.2.4), and analyzes severe historical drawdowns (Q3.2.5)


source("code/05_data_cleaning.R")

library(zoo)

# Sort target underlying datasets sequentially by date records
portfolio_data <- X30_Industry_Portfolios %>% arrange(date)
factors_data <- FF_Research_Data_Factors %>% arrange(date)

# Isolate risk-free asset return rates from the macro factor table
rf_data <- factors_data %>%
  select(date, matches("(?i)^rf$")) %>%
  rename(rf = 2)

# Calculate asset net excess returns over risk-free base rates
industry_returns <- portfolio_data %>%
  left_join(rf_data, by = "date") %>%
  mutate(across(-c(date, rf), ~ . - rf)) %>%
  select(-rf) %>%
  filter(date >= 192607)

# Average historical returns across a rolling twelve month backward window
rolling_mom <- industry_returns %>%
  mutate(across(-date, ~ rollmean(lag(.), k = 12, fill = NA, align = "right"))) %>%
  filter(date >= 192707)

# Reshape the output data matrix to long form for ranking operations
mom_long <- rolling_mom %>%
  pivot_longer(cols = -date, names_to = "Industry", values_to = "Past_12M_Return") %>%
  filter(!is.na(Past_12M_Return))

# Assign sorting ranks within each independent monthly snapshot
mom_ranked <- mom_long %>%
  group_by(date) %>%
  mutate(Rank = rank(Past_12M_Return, ties.method = "first")) %>%
  ungroup()

# Summarize the lifetime average rank locations across industry categories
industry_rank_summary <- mom_ranked %>%
  group_by(Industry) %>%
  summarize(Avg_Rank = mean(Rank)) %>%
  arrange(Avg_Rank)

lowest_industry <- slice_head(industry_rank_summary, n = 1)
highest_industry <- slice_tail(industry_rank_summary, n = 1)

cat("Industry rank standings\n")
print(paste("Lowest average rank industry:", lowest_industry$Industry, "with an average rank of:", round(lowest_industry$Avg_Rank, 2)))
print(paste("Highest average rank industry:", highest_industry$Industry, "with an average rank of:", round(highest_industry$Avg_Rank, 2)))

# Measure the absolute chronological differences in industry ranks between months
turnover_analysis <- mom_ranked %>%
  group_by(Industry) %>%
  arrange(date) %>%
  mutate(Rank_Change = abs(Rank - lag(Rank))) %>%
  ungroup() %>%
  filter(!is.na(Rank_Change))

avg_monthly_rank_shift <- mean(turnover_analysis$Rank_Change)
print(paste("Average monthly rank position shift per industry:", round(avg_monthly_rank_shift, 2)))

# Map date fields to standard datetime formatting profiles for plotting
mom_ranked_plot_df <- mom_ranked %>%
  mutate(Date_Parsed = as.Date(paste0(date, "01"), format = "%Y%m%d"))

target_industry <- "Games"
ggplot(filter(mom_ranked_plot_df, Industry == target_industry), aes(x = Date_Parsed, y = Rank)) +
  geom_line(color = "blue", alpha = 0.6) +
  scale_y_continuous(breaks = seq(1, 30, by = 5)) +
  labs(
    title = paste("Rank evolution of", target_industry, "industry"),
    subtitle = "Monthly rank positions (1 = Worst 12m return, 30 = Best 12m return)",
    x = "Year",
    y = "Rank position"
  ) +
  theme_minimal(base_size = 12)

ggplot(industry_rank_summary, aes(x = reorder(Industry, Avg_Rank), y = Avg_Rank)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Average rank across all 30 industries, based on past average returns",
    x = "Industry",
    y = "Average historical rank"
  ) +
  theme_minimal()


# Reshape asset historical returns to facilitate data join operations
industry_returns_long <- industry_returns %>%
  pivot_longer(cols = -date, names_to = "Industry", values_to = "Excess_Return")

# Merge dynamic factor ranks with contemporaneous performance data
portfolio_mapping <- mom_ranked %>%
  left_join(industry_returns_long, by = c("date", "Industry"))

# Classify industries into discrete winner and loser equal-weighted buckets
strategy_returns <- portfolio_mapping %>%
  mutate(Portfolio_Type = if_else(Rank > 15, "Winner", "Loser")) %>%
  group_by(date, Portfolio_Type) %>%
  summarize(Port_Excess_Return = mean(Excess_Return, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Portfolio_Type, values_from = Port_Excess_Return)

# Pull down historical benchmark equity premium returns
market_benchmark <- FF_Research_Data_Factors %>%
  arrange(date) %>%
  filter(date >= 192707) %>%
  select(date, `Mkt-RF`)

# Construct a single reference table holding all structural returns
final_returns_df <- strategy_returns %>%
  left_join(market_benchmark, by = "date")

# Compute mean, deviation, and Sharpe statistics across strategy vectors
performance_summary <- final_returns_df %>%
  summarize(
    Winner_Avg_Excess    = mean(Winner),
    Winner_Std_Dev       = sd(Winner),
    Winner_Monthly_SR    = Winner_Avg_Excess / Winner_Std_Dev,
    Winner_Annualized_SR = Winner_Monthly_SR * sqrt(12),
    Loser_Avg_Excess     = mean(Loser),
    Loser_Std_Dev        = sd(Loser),
    Loser_Monthly_SR     = Loser_Avg_Excess / Loser_Std_Dev,
    Loser_Annualized_SR  = Loser_Monthly_SR * sqrt(12),
    Market_Avg_Excess    = mean(`Mkt-RF`),
    Market_Std_Dev       = sd(`Mkt-RF`),
    Market_Monthly_SR    = Market_Avg_Excess / Market_Std_Dev,
    Market_Annualized_SR = Market_Monthly_SR * sqrt(12)
  )

cat("\nWinner portfolio performance metrics\n")
cat(paste("Avg monthly excess return:  ", round(performance_summary$Winner_Avg_Excess, 4), "\n"))
cat(paste("Monthly standard deviation: ", round(performance_summary$Winner_Std_Dev, 4), "\n"))
cat(paste("Monthly sharpe ratio:        ", round(performance_summary$Winner_Monthly_SR, 4), "\n"))
cat(paste("Annualized sharpe ratio:     ", round(performance_summary$Winner_Annualized_SR, 4), "\n\n"))

cat("Loser portfolio performance metrics\n")
cat(paste("Avg monthly excess return:  ", round(performance_summary$Loser_Avg_Excess, 4), "\n"))
cat(paste("Monthly standard deviation: ", round(performance_summary$Loser_Std_Dev, 4), "\n"))
cat(paste("Monthly sharpe ratio:        ", round(performance_summary$Loser_Monthly_SR, 4), "\n"))
cat(paste("Annualized sharpe ratio:     ", round(performance_summary$Loser_Annualized_SR, 4), "\n\n"))

cat("Market index benchmark\n")
cat(paste("Annualized sharpe ratio:     ", round(performance_summary$Market_Annualized_SR, 4), "\n"))


# Reintegrate risk-free rates to build total nominal returns
rf_reconstruct <- FF_Research_Data_Factors %>%
  arrange(date) %>%
  filter(date >= 192707) %>%
  select(date, RF)

total_returns_df <- final_returns_df %>%
  left_join(rf_reconstruct, by = "date") %>%
  mutate(
    Winner_Total = Winner + RF,
    Loser_Total = Loser + RF,
    Market_Total = `Mkt-RF` + RF,
    Ind_Mom_Total = Winner_Total - Loser_Total
  )

# Calculate compound product curves tracking growth from base investment
cumulative_df <- total_returns_df %>%
  arrange(date) %>%
  mutate(
    Cum_Winner = cumprod(1 + Winner_Total),
    Cum_Loser = cumprod(1 + Loser_Total),
    Cum_Market = cumprod(1 + Market_Total),
    Cum_Ind_Mom = cumprod(1 + Ind_Mom_Total),
    Date_Parsed = as.Date(paste0(date, "01"), format = "%Y%m%d")
  )

# Transform structural accumulation vectors into tidy plotting shapes
cumulative_long <- cumulative_df %>%
  select(Date_Parsed, Cum_Winner, Cum_Loser, Cum_Market, Cum_Ind_Mom) %>%
  pivot_longer(cols = -Date_Parsed, names_to = "Portfolio", values_to = "Cumulative_Value") %>%
  mutate(Portfolio = case_when(
    Portfolio == "Cum_Winner" ~ "Winner Portfolio",
    Portfolio == "Cum_Loser" ~ "Loser Portfolio",
    Portfolio == "Cum_Market" ~ "Market Index",
    Portfolio == "Cum_Ind_Mom" ~ "Long/Short Ind-Mom"
  ))

ggplot(cumulative_long, aes(x = Date_Parsed, y = Cumulative_Value, color = Portfolio)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(labels = scales::dollar_format()) +
  labs(
    title = "Cumulative returns of $1 invested",
    x = "Year",
    y = "Log cumulative wealth ($)",
    color = "Strategy"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Filter timelines to examine strategy drawdowns during the 2009 shift
crash_2009 <- total_returns_df %>%
  filter(str_detect(as.character(date), "^2009")) %>%
  select(date, Ind_Mom_Total, Market_Total, Winner_Total, Loser_Total)

cat("\nIndustry momentum performance in 2009\n")
print(crash_2009)

# Isolate the worst standalone performance periods across history
worst_months <- total_returns_df %>%
  select(date, Ind_Mom_Total, Market_Total, Winner_Total, Loser_Total) %>%
  arrange(Ind_Mom_Total) %>%
  slice_head(n = 10)

cat("\nTop 10 worst months for l/s ind-mom\n")
print(worst_months)
