# Author: Jean de Biolley


# The following code downloads normalized monthly price charts for five technology stocks (GOOG, CSCO, LOGI, AMZN, AAPL), 
# constructs the 5-asset efficient frontier combining the Global Minimum Variance Portfolio (GMVP) and the Tangency Portfolio, 
# computes the tangency portfolio weights, and plots the efficient frontier together with the Capital Allocation Line (CAL). 


import yfinance as yf
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tickers = ['GOOG', 'CSCO', 'LOGI','AMZN','AAPL']
start_date = '2015-01-01'

data = yf.download(tickers, start=start_date, auto_adjust=True)['Close']

monthly_prices = data.resample('BM').last()
monthly_returns = monthly_prices.pct_change().dropna()

print("First 5 rows of monthly returns:\n", monthly_returns.head())

# Compute the vector of sample means
sample_means = monthly_returns.mean()
print("\nSample Means of Monthly Returns:\n", sample_means)

# Compute the matrix of sample variances and covariances
sample_cov_matrix = monthly_returns.cov()
print("\nSample Variance-Covariance Matrix of Monthly Returns:\n", sample_cov_matrix)

# Normalize prices (base = 100 at the start date)
normalized_prices = (monthly_prices / monthly_prices.iloc[0]) * 100

# Plot the normalized price chart
plt.figure(figsize=(10, 6))
for ticker in tickers:
    plt.plot(normalized_prices.index, normalized_prices[ticker], label=ticker)

plt.title('Normalized Monthly Price Chart (Base 100)')
plt.xlabel('Date')
plt.ylabel('Normalized Price')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Mean returns and variances-covariances of the assets
num_assets = len(tickers)
mean_returns = sample_means.values
cov_matrix = sample_cov_matrix.values
inv_cov_matrix = np.linalg.inv(cov_matrix)
ones = np.ones(num_assets)

# Global Minimum Variance (GMV) Portfolio
w_gmv = (inv_cov_matrix @ ones) / (ones.T @ inv_cov_matrix @ ones)
gmv_mean = w_gmv.T @ mean_returns
gmv_var = w_gmv.T @ cov_matrix @ w_gmv
gmv_stdev = np.sqrt(gmv_var)

print("\n--- Global Minimum Variance Portfolio ---")
for i, ticker in enumerate(tickers):
    print(f"Weight {ticker}: {w_gmv[i]:.4f}")
print(f"Expected Return: {gmv_mean:.6f}")
print(f"Standard Deviation: {gmv_stdev:.6f}")

# Tangency Portfolio
rf = 0.00025
excess_returns = mean_returns - rf
w_tan = (inv_cov_matrix @ excess_returns) / (ones.T @ inv_cov_matrix @ excess_returns)
tan_mean = w_tan.T @ mean_returns
tan_var = w_tan.T @ cov_matrix @ w_tan
tan_stdev = np.sqrt(tan_var)

print("\n--- Tangency Portfolio ---")
for i, ticker in enumerate(tickers):
    print(f"Weight {ticker}: {w_tan[i]:.4f}")
print(f"Expected Return: {tan_mean:.6f}")
print(f"Standard Deviation: {tan_stdev:.6f}")

# Efficient frontier (combinations of GMV and Tangency)
weights = np.linspace(-1.0, 2.5, 200)
port_returns = []
port_stdevs = []
for w in weights:
    w_port = w * w_tan + (1 - w) * w_gmv
    port_returns.append(w_port.T @ mean_returns)
    port_stdevs.append(np.sqrt(w_port.T @ cov_matrix @ w_port))

port_returns = np.array(port_returns)
port_stdevs = np.array(port_stdevs)

# Capital Allocation Line (CAL)
cal_stdevs = np.linspace(0, max(port_stdevs)*1.1, 100)
sharpe_ratio = (tan_mean - rf) / tan_stdev
cal_returns = rf + sharpe_ratio * cal_stdevs

# Plot Efficient Frontier
plt.figure(figsize=(10, 6))
plt.plot(port_stdevs, port_returns, label='Efficient Frontier', color='blue')
plt.scatter(gmv_stdev, gmv_mean, color='red', marker='o', s=100, label='GMV Portfolio', zorder=5)
plt.scatter(tan_stdev, tan_mean, color='green', marker='*', s=150, label='Tangency Portfolio', zorder=5)
plt.plot(cal_stdevs, cal_returns, label='Capital Allocation Line (CAL)', color='orange', linestyle='--')

# Plot the individual assets
for i, ticker in enumerate(tickers):
    std = np.sqrt(cov_matrix[i, i])
    mean = mean_returns[i]
    plt.scatter(std, mean, color='black', marker='x', s=100, zorder=5)
    plt.annotate(ticker, (std, mean), textcoords="offset points", xytext=(0,10), ha='center')

plt.title('Efficient Frontier and Portfolio Analysis')
plt.xlabel('Standard Deviation (Risk)')
plt.ylabel('Expected Return')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()