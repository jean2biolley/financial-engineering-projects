# Author: Jean de Biolley

# The following code downloads European equity (VEUSX) and long-term bond (VBLIX) price data via yfinance, 
# generates normalized monthly performance charts, computes sample return and covariance matrices, 
# and constructs the two-asset Markowitz efficient frontier, Global Minimum Variance Portfolio (GMVP), 
# and Capital Allocation Line (CAL). 


import yfinance as yf
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tickers = ['VEUSX', 'VBLIX']
start_date = '2015-01-01'

# Download with auto_adjust to avoid the KeyError
data = yf.download(tickers, start=start_date, auto_adjust=True)['Close']

# The rest of your logic remains the same
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
asset1 = 'VEUSX'
asset2 = 'VBLIX'

mean1 = sample_means[asset1]
mean2 = sample_means[asset2]

var1 = sample_cov_matrix.loc[asset1, asset1]
var2 = sample_cov_matrix.loc[asset2, asset2]
cov12 = sample_cov_matrix.loc[asset1, asset2]

# Efficient frontier
weights1 = np.linspace(-0.5, 1.5, 200)
port_returns = weights1 * mean1 + (1 - weights1) * mean2
port_vars = (weights1**2) * var1 + ((1 - weights1)**2) * var2 + 2 * weights1 * (1 - weights1) * cov12
port_stdevs = np.sqrt(port_vars)

# Global Minimum Variance (GMV) Portfolio
w1_gmv = (var2 - cov12) / (var1 + var2 - 2 * cov12)
w2_gmv = 1 - w1_gmv

gmv_mean = w1_gmv * mean1 + w2_gmv * mean2
gmv_var = (w1_gmv**2) * var1 + (w2_gmv**2) * var2 + 2 * w1_gmv * w2_gmv * cov12
gmv_stdev = np.sqrt(gmv_var)

print("\n--- Global Minimum Variance Portfolio ---")
print(f"Weight {asset1}: {w1_gmv:.4f}")
print(f"Weight {asset2}: {w2_gmv:.4f}")
print(f"Expected Return: {gmv_mean:.6f}")
print(f"Standard Deviation: {gmv_stdev:.6f}")

# Tangency Portfolio
rf = 0.00025
num1 = (mean1 - rf) * var2 - (mean2 - rf) * cov12
den = (mean1 - rf) * var2 + (mean2 - rf) * var1 - (mean1 - rf + mean2 - rf) * cov12
w1_tan = num1 / den
w2_tan = 1 - w1_tan

tan_mean = w1_tan * mean1 + w2_tan * mean2
tan_var = (w1_tan**2) * var1 + (w2_tan**2) * var2 + 2 * w1_tan * w2_tan * cov12
tan_stdev = np.sqrt(tan_var)

print("\n--- Tangency Portfolio ---")
print(f"Weight {asset1}: {w1_tan:.4f}")
print(f"Weight {asset2}: {w2_tan:.4f}")
print(f"Expected Return: {tan_mean:.6f}")
print(f"Standard Deviation: {tan_stdev:.6f}")

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
plt.scatter([np.sqrt(var1), np.sqrt(var2)], [mean1, mean2], color='black', marker='x', s=100, zorder=5)
plt.annotate(asset1, (np.sqrt(var1), mean1), textcoords="offset points", xytext=(0,10), ha='center')
plt.annotate(asset2, (np.sqrt(var2), mean2), textcoords="offset points", xytext=(0,10), ha='center')

plt.title('Efficient Frontier and Portfolio Analysis')
plt.xlabel('Standard Deviation (Risk)')
plt.ylabel('Expected Return')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()