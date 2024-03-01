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
x_predictions = np.array([item['calculatedDiffusionCoefficient'] for item in data_sorted])

# Since we are dropping the weights (or making them equal to 1), we use OLS instead of WLS
X = sm.add_constant(x_predictions)  # Add a constant to the predictor variable array for the intercept

# Perform Ordinary Least Squares (OLS) regression
model = sm.OLS(y_means, X).fit()

print(model.summary())

# Optionally, you can print the coefficients directly
print(f"Coefficient for calculatedDiffusionCoefficient: {model.params[1]:.2f}, Intercept: {model.params[0]:.2f}")

# Generate a range of x values for plotting the model line
x_range = np.linspace(min(x_predictions), max(x_predictions), 100)
X_range = sm.add_constant(x_range)  # Add constant for intercept

# Calculate predicted y values across the x range
y_pred_range = model.predict(X_range)



# Prepare the data for plotting
conditions = range(1, len(data_sorted) + 1)  # Simple numeric labels for conditions

# Plotting
plt.figure(figsize=(10, 6))

# Plot individual y_samples for each condition
for i, item in enumerate(data_sorted, start=1):
    y_samples = item['experimentalDiffusionCoefficient']
    # Scatter plot for y_samples with a lower zorder to ensure they are behind the main markers
    plt.scatter([i] * len(y_samples), y_samples, alpha=0.7, edgecolor='blue', facecolor='none', label='Exp Samples' if i == 1 else "", zorder=2)

# Plot x_prediction with higher zorder to ensure it is on top
plt.scatter(conditions, x_predictions, marker='D', color='red', edgecolor='none', label='X Prediction', zorder=3)

# Customizing the plot
plt.xlabel('Condition')
plt.ylabel('Value')
plt.title('Plot of Y Means with Std Dev, X Predictions, and Y Samples')
plt.xticks(conditions)  # Set x-ticks to condition labels
plt.legend()
plt.grid(True)

# Save plot
plt.savefig('/home/radossio/MDwMD/validation_set/plot3.png', bbox_inches='tight')




# Plotting the model
plt.figure(figsize=(8, 8))  # Adjust figure size to make it square
plt.plot(x_range, y_pred_range, 'g-', label='Regression Line', zorder=4)

# Scatter plot for y_means vs. x_predictions to show the relationship
plt.scatter(x_predictions, y_means, marker='o', color='blue', edgecolor='black', label='Data Points', zorder=5)

# Determine the maximum range to set equal limits for both axes
max_range = max(max(x_predictions) - min(x_predictions), max(y_means) - min(y_means))
x_lim = (min(x_predictions) - 0.1 * max_range, max(x_predictions) + 0.1 * max_range)
y_lim = (min(y_means) - 0.1 * max_range, max(y_means) + 0.1 * max_range)

# Set the same limits for both axes to ensure delta x equals delta y
plt.xlim(x_lim)
plt.ylim(y_lim)

# Customizing the plot
plt.xlabel('Calculated Diffusion Coefficient')
plt.ylabel('Experimental Diffusion Coefficient (Mean)')
plt.title('Relationship between Calculated and Experimental Diffusion Coefficients')
plt.legend()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')  # Ensure the aspect ratio is equal

# Show or save the plot
plt.savefig('plot4.png', bbox_inches='tight')  # Adjust the path as necessary

