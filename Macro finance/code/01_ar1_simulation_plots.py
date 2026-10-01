# Author : Jean de Biolley


# The following code simulates and plots the conditional mean (Q1.c.i), 1 standard deviation bands (Q1.c.ii), 
# Chebyshev versus normal prediction cones (Q1.c.iii), and 20 realization paths for 
# an AR(1) stock price process (Q1.c.iv).   


############# Q1.c.i##############

import numpy as np
import matplotlib.pyplot as plt

rho = 0.95
x1 = 1
t_max = 20

t_values = np.arange(1, t_max + 1)

conditional_means = [x1 * (rho**(t - 1)) for t in t_values]

plt.figure(figsize=(10, 6))
plt.plot(t_values, conditional_means, marker='o', linestyle='-', color='b', label=f'Conditional Mean (rho={rho})')

plt.title('Conditional Mean $E_1(x_t)$ for AR(1) Process')
plt.xlabel('Time (t)')
plt.ylabel('Expected Value')
plt.xticks(t_values)
plt.grid(True, linestyle='--', alpha=0.7)
plt.axhline(0, color='black', linewidth=0.8) # Long-run mean
plt.legend()

plt.show()


###########3 Q1.c.ii ##############

sigma_eps = 0.25

cond_mean = x1 * (rho**(t_values - 1))

cond_var = (sigma_eps**2) * (1 - rho**(2 * (t_values - 1))) / (1 - rho**2)
cond_std = np.sqrt(cond_var)

plt.figure(figsize=(10, 6))
plt.plot(t_values, cond_mean, marker='o', label='Conditional Mean $E_1(x_t)$', color='blue')

plt.fill_between(t_values, cond_mean - cond_std, cond_mean + cond_std, 
                 color='blue', alpha=0.15, label=r'$\pm 1$ Cond. Std. Dev. Band')

plt.axhline(0, color='black', linestyle='--', linewidth=0.8)
plt.title(r'AR(1) Conditional Mean and $\pm 1$ Standard Deviation Bands')
plt.xlabel('Time (t)')
plt.ylabel('$x_t$')
plt.xticks(t_values)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

plt.show()

####### Q1.c.iii ##############


k95 = 1 / np.sqrt(0.05)  
k90 = 1 / np.sqrt(0.10)  

plt.figure(figsize=(10, 6))

plt.plot(t_values, cond_mean, marker='o', color='blue', label='Conditional Mean')
plt.fill_between(t_values, cond_mean - cond_std, cond_mean + cond_std, 
                 color='blue', alpha=0.15, label=r'$\pm 1$ Cond. Std. Dev. Band')

plt.fill_between(t_values, cond_mean - k95*cond_std, cond_mean + k95*cond_std, 
                 color='red', alpha=0.1, label='Chebyshev 95% Band')

plt.fill_between(t_values, cond_mean - k90*cond_std, cond_mean + k90*cond_std, 
                 color='red', alpha=0.2, label='Chebyshev 90% Band')

plt.title('AR(1) Prediction Cones (Chebyshev Inequality)')
plt.xlabel('Time (t)')
plt.ylabel('$x_t$')
plt.axhline(0, color='black', linestyle='--', alpha=0.5)
plt.xticks(t_values)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()

plt.show()

####### Q1.c.iv ##############

num_simulations = 20

z95 = 1.96  # 95% confidence
z90 = 1.645 # 90% confidence

all_sims = np.zeros((num_simulations, t_max))
all_sims[:, 0] = x1 

for i in range(num_simulations):
    for t in range(1, t_max):
        epsilon = np.random.normal(0, sigma_eps)
        all_sims[i, t] = rho * all_sims[i, t-1] + epsilon

plt.figure(figsize=(12, 7))

plt.fill_between(t_values, cond_mean - z95*cond_std, cond_mean + z95*cond_std, 
                 color='blue', alpha=0.1, label='Normal 95% Prediction Cone')
plt.fill_between(t_values, cond_mean - z90*cond_std, cond_mean + z90*cond_std, 
                 color='blue', alpha=0.2, label='Normal 90% Prediction Cone')

for i in range(num_simulations):
    plt.plot(t_values, all_sims[i, :], color='black', alpha=0.2, linewidth=1)

plt.plot(t_values, all_sims[0, :], color='black', alpha=0.6, linewidth=1, label='Simulated Paths ($x_t$)')
plt.plot(t_values, cond_mean, marker='o', color='blue', linewidth=2, label='Conditional Mean $E_1(x_t)$')

plt.title(f'AR(1) Simulations vs. Theoretical Mean and Normal Prediction Cones\n($\\rho={rho}, \\sigma_\\epsilon={sigma_eps}$)')
plt.xlabel('Time (t)')
plt.ylabel('$x_t$')
plt.axhline(0, color='black', linestyle='--', alpha=0.5)
plt.xticks(t_values)
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(loc='upper right')

plt.show()