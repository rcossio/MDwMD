import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import json

# Load mocked data from JSON file
with open('diffusion_data.json', 'r') as f:
    data = json.load(f)

# Prepare the data
data_sorted = sorted(data, key=lambda item: np.mean(item['experimentalDiffusionCoefficient']))
y_means = np.array([np.mean(item['experimentalDiffusionCoefficient']) for item in data_sorted])
accession_numbers = [item['accessionNumber'] for item in data_sorted]  # Collect accession numbers

# Prepare the data for plotting
conditions = np.arange(len(data_sorted))  # Use numpy to create an array of conditions

# Plotting
plt.figure(figsize=(10, 6))

# Plot individual y_samples for each condition
for i, item in enumerate(data_sorted):
    y_samples = item['experimentalDiffusionCoefficient']
    # Scatter plot for y_samples with a lower zorder to ensure they are behind the main markers
    plt.scatter([i] * len(y_samples), y_samples, alpha=0.7, edgecolor='blue', facecolor='none', label='Exp Samples' if i == 0 else "", zorder=2)

# Customizing the plot
plt.xlabel('Accession Number')
plt.ylabel('Experimental Diffusion Coefficient')
plt.title('Plot of Experimental Diffusion Coefficients')
plt.xticks(conditions, accession_numbers, rotation='vertical')  # Apply accession numbers as x-tick labels

plt.legend()
plt.grid(True)

# Adjust layout to prevent clipping of tick-labels
plt.tight_layout()

# Save plot
plt.savefig('plot5.png', bbox_inches='tight')
