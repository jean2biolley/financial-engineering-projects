# Author: Jean de Biolley

# The following code simulates and plots the impulse response function of returns and 
# prices for an MA(1) model.

import matplotlib.pyplot as plt
import numpy as np

deltas = [0.5, 0, -0.5]
t_max = 5
t = np.arange(0, t_max + 1)


epsilon = np.zeros(t_max + 2) 
epsilon[1] = 1 

def get_ma1_irf(delta):
    r = np.zeros(len(t))
    for i in range(len(t)):
        idx = t[i]
        val = epsilon[idx] + delta * epsilon[idx-1] if idx > 0 else 0
        r[i] = val
    return r

# Plot 1: Return IRF
plt.figure(figsize=(8, 5))
for d in deltas:
    r_irf = get_ma1_irf(d)
    plt.plot(t, r_irf, marker='o', label=rf'$\delta = {d}$')

plt.title('Impulse Response Function of Returns ($r_t$)', fontsize=14)
plt.xlabel('Time (t)', fontsize=12)
plt.ylabel('Return Response', fontsize=12)
plt.xticks(t)
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('return_irf.png')
plt.close()

# Plot 2: Price IRF (Cumulative Log Returns)
plt.figure(figsize=(8, 5))
for d in deltas:
    r_irf = get_ma1_irf(d)
    p_irf = np.cumsum(r_irf)
    plt.plot(t, p_irf, marker='s', label=rf'$\delta = {d}$')

plt.title('Impulse Response Function of Price ($p_t$)', fontsize=14)
plt.xlabel('Time (t)', fontsize=12)
plt.ylabel('Cumulative Log Return (Price)', fontsize=12)
plt.xticks(t)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('price_irf.png')
plt.close()