import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import json

# Load mocked data from JSON file
with open('diffusion_data.json', 'r') as f:
    data = json.load(f)

# Prepare the data
data_sorted = sorted(data, key=lambda item: np.mean(item['experiment_coef']))

# Calculating error in a pythonic way
errors = []
exp_means = []
for item in data_sorted:
    calculated_coefs = np.array([coef for coef in item['calculated_coef'] if coef is not None])
    if calculated_coefs.size > 0:
        exp_mean = np.mean(item['experiment_coef'])
        exp_means.append(exp_mean)
        errors.append(np.mean((calculated_coefs - exp_mean) ** 2))
mean_squared_error = np.mean(errors)
rmse = np.sqrt(mean_squared_error)
print("Error: ", rmse)

# Plotting
plt.figure(figsize=(10, 6))
plt.grid(True, zorder=0)

for item in data_sorted:
    exp_mean = np.mean(item['experiment_coef'])
    for y_sample in [coef for coef in item['calculated_coef'] if coef is not None]:
        single_label = 'Calculated Coefficient' if item == data_sorted[0] and y_sample == item['calculated_coef'][0] else None
        plt.scatter(exp_mean, y_sample, alpha=0.7, edgecolor='red', facecolor='none', zorder=3, label=single_label)

x_values = np.array([20, 170])
plt.plot(x_values, x_values, 'k-', zorder=2, label='y = x', linewidth=1)
plt.plot(x_values, x_values+rmse, 'k--', zorder=2, label='+Error', linewidth=1)
plt.plot(x_values, x_values-rmse, 'k--', zorder=2, label='-Error', linewidth=1)

plt.ylabel('Calculated Coefficient')
plt.xlabel('Experimental Coefficient')
plt.gca().set_aspect('equal', adjustable='box')
plt.ylim(20, 180)
plt.xlim(20, 180)
plt.tight_layout()
plt.legend()
plt.savefig('plot6.png', bbox_inches='tight')

# Finally plot errors against the mean experimental coefficient
plt.figure(figsize=(6, 4))
plt.grid(True, zorder=0)
plt.xlim(20, 180)
x_values = np.array([20, 170])
plt.plot(x_values, [rmse]*len(x_values), 'k--', zorder=2, label='Error', linewidth=1)
plt.scatter(exp_means, np.sqrt(errors), edgecolor='none', facecolor='red', zorder=2)
plt.ylabel('Mean Squared Error')
plt.xlabel('Experimental Coefficient')
plt.tight_layout()
plt.savefig('plot7.png', bbox_inches='tight')
