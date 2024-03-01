
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


mocked_data = [
    {'y_samples': [2.55, 2.84, 3.18, 3.52], 'x_prediction': 2.70},
    {'y_samples': [4.15, 4.26, 3.92, 3.88, 4.64], 'x_prediction': 4.19},
    {'y_samples': [6.50, 5.62, 6.28], 'x_prediction': 5.57},
    {'y_samples': [7.25, 6.56, 7.25, 8.06, 6.98], 'x_prediction': 7.21},
    {'y_samples': [9.90, 5.71, 9.16, 8.42, 8.03], 'x_prediction': 8.40},
    {'y_samples': [9.74, 8.14], 'x_prediction': 9.31},
    {'y_samples': [10.80, 11.33, 12.36], 'x_prediction': 11.06},
    {'y_samples': [12.18, 12.09, 12.18], 'x_prediction': 12.26},
    {'y_samples': [14.08, 13.82, 13.43], 'x_prediction': 13.53},
    {'y_samples': [15.43, 15.08, 15.82, 14.41, 14.72], 'x_prediction': 15.12}
]

# Prepare the data
mocked_data_sorted = sorted(mocked_data, key=lambda item: np.mean(item['y_samples']))
y_means = np.array([np.mean(item['y_samples']) for item in mocked_data_sorted])
x_predictions = np.array([item['x_prediction'] for item in mocked_data_sorted])
y_stdevs = np.array([np.std(item['y_samples'], ddof=1) for item in mocked_data_sorted])

# The weights for WLS are the inverse of the variance (1 / standard deviation^2)
weights = 1 / y_stdevs**2

# Add a constant to the predictor variable array for the intercept
X = sm.add_constant(x_predictions)

# Perform Weighted Least Squares (WLS) regression
model = sm.WLS(y_means, X, weights=weights).fit()

print(model.summary())

# Optionally, you can print the coefficients directly
print(f"Coefficient for x_prediction: {model.params[1]:.2f}, Intercept: {model.params[0]:.2f}")


# Prepare the data for plotting
conditions = range(1, len(mocked_data_sorted) + 1)  # Simple numeric labels for conditions

# Plotting
plt.figure(figsize=(10, 6))

# Plot individual y_samples for each condition
for i, item in enumerate(mocked_data_sorted, start=1):
    y_samples = item['y_samples']
    # Scatter plot for y_samples with a lower zorder to ensure they are behind the main markers
    plt.scatter([i] * len(y_samples), y_samples, alpha=0.7, edgecolor='blue', facecolor='none', label='Y Samples' if i == 1 else "", zorder=2)

# Plot mean ± stdev as error bars
#plt.errorbar(conditions, y_means, yerr=y_stdevs, fmt='o', capsize=5, color='blue', label='Mean ± Std Dev', zorder=1)

# Replace errorbar with a box indicating the standard deviation and mean
for i, (mean, stdev) in enumerate(zip(y_means, y_stdevs), start=1):
    # Draw a rectangle for the standard deviation range
    plt.gca().add_patch(patches.Rectangle((i - 0.2, mean - stdev), 0.4, 2 * stdev, color='blue', alpha=0.2, zorder=1))
    # Mark the mean
    plt.plot(i, mean, 'o', color='blue', label='Mean' if i == 1 else "")

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
plt.savefig('/home/radossio/MDwMD/validation_set/plot.png', bbox_inches='tight')

# Generate a range of x values for plotting the model line
x_range = np.linspace(min(x_predictions), max(x_predictions), 100)
X_range = sm.add_constant(x_range)  # Add constant for intercept

# Calculate predicted y values across the x range
y_pred_range = model.predict(X_range)

# Plotting the model
plt.figure(figsize=(10, 6))
plt.plot(x_range, y_pred_range, 'g-', label='WLS Model', zorder=4)

# Scatter plot for y_means vs. x_predictions to show the relationship
plt.scatter(x_predictions, y_means, marker='o', color='blue', edgecolor='black', label='Y Means', zorder=5)

# Customizing the plot
plt.xlabel('X Prediction')
plt.ylabel('Y Mean')
plt.title('Relationship between X Predictions and Y Means with Model')
plt.legend()
plt.grid(True)

# Show or save the plot
plt.savefig('/home/radossio/MDwMD/validation_set/plot2.png', bbox_inches='tight')