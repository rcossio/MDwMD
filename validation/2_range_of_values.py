import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import json

# Load mocked data from JSON file
with open('diffusion_data.json', 'r') as f:
    data = json.load(f)

# Prepare the data
data_sorted = sorted(data, key=lambda item: np.mean(item['experiment_coef']))
names = [item['name'][0:20] for item in data_sorted] 

# Prepare the data for plotting
xtics = np.arange(len(data_sorted))  

# Plotting
plt.figure(figsize=(10, 6))
plt.grid(True, zorder=0)

# Plot individual y_samples for each condition
for i, item in enumerate(data_sorted):
    y_samples = item['experiment_coef']
    plt.scatter([i] * len(y_samples), y_samples, alpha=0.7, edgecolor='blue', facecolor='none', label='Exp Samples' if i == 0 else "", zorder=2)

# Customizing the plot
plt.ylabel('Diffusion Coefficient')
plt.xticks(xtics, names, rotation='vertical') 

plt.legend()

# Adjust layout to prevent clipping of tick-labels
plt.tight_layout()

# Save plot
plt.savefig('plot1.png', bbox_inches='tight')

# Plotting the calculated coefficients
for i, item in enumerate(data_sorted):
    x_samples = item['calculated_coef']
    plt.scatter([i] * len(x_samples), x_samples, alpha=0.7, edgecolor='red', facecolor='none', label='Calc Samples' if i == 0 else "", zorder=3, marker='s')
plt.legend()
plt.savefig('plot2.png', bbox_inches='tight')


y_means = [np.mean(item['experiment_coef']) for item in data_sorted]
y_means_weights = [len(item['experiment_coef']) for item in data_sorted]

plt.figure(figsize=(10, 6))
plt.grid(True, zorder=0)
for i,item in enumerate(y_means):
    plt.scatter(i, item, edgecolor='blue', facecolor='blue', label='Exp Average' if i == 0 else "", zorder=2, marker='D', s=y_means_weights[i]*12)
for i, item in enumerate(data_sorted):
    x_samples = item['calculated_coef']
    plt.scatter([i] * len(x_samples), x_samples, alpha=0.7, edgecolor='red', facecolor='none', label='Calc Samples' if i == 0 else "", zorder=3, marker='s')
    plt.ylabel('Diffusion Coefficient')
plt.xticks(xtics, names, rotation='vertical') 
plt.legend()
plt.tight_layout()
plt.savefig('plot3.png', bbox_inches='tight')