import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Load data
data = np.loadtxt('proj.xvg')
em_data = np.loadtxt('em_proj.xvg')
nvt_data = np.loadtxt('nvt_proj.xvg')
npt_data = np.loadtxt('npt_proj.xvg')
rst200_data = np.loadtxt('rst_200_proj.xvg')
rst50_data   = np.loadtxt('rst_50_proj.xvg')
rst10_data   = np.loadtxt('rst_10_proj.xvg')
rst0_data   = np.loadtxt('rst_0_proj.xvg')
md_data  = np.loadtxt('md_proj.xvg')


# Data columns for the main data
x = data[:, 0]
y = data[:, 1]

# Data columns for the additional data
# Ensure it's treated as a two-dimensional array with a single point
em_x, em_y = em_data  # Direct unpacking since it's just one pair
nvt_x, nvt_y = nvt_data 
npt_x, npt_y = npt_data
rst200_x, rst200_y = rst200_data
rst50_x, rst50_y = rst50_data
rst10_x, rst10_y = rst10_data
rst0_x, rst0_y = rst0_data
md_x, md_y = md_data


# Create a subset of data points, taking every second point to plot
sub_x = x[::1]  # Every other x value
sub_y = y[::1]  # Every other y value

# Determine ranges for the plot
max_val = max(np.max(np.abs(x)), np.max(np.abs(y)), abs(em_x), abs(em_y))
range_val = max_val

# Define the plot range
plot_range = [-range_val * 1.05, range_val * 1.05]

# Create a custom color map
cmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'blue', 'violet'])

# Normalize for colormap: create a continuous normalized value array for line
norm = plt.Normalize(0, len(sub_x))

# Plotting
plt.figure(figsize=(10, 8))  # Adjusted size for space for color bar

# Line plot with color mapping for each segment
for i in range(len(sub_x) - 1):
    plt.plot(sub_x[i:i+2], sub_y[i:i+2], color=cmap(norm(i)), linewidth=0.5, zorder=1)

# Scatter plot every 100 points with custom size and color palette
sample_indices = np.arange(0, len(x), 100)
sc = plt.scatter(x[sample_indices], y[sample_indices], c=norm(np.linspace(0, 1, len(sample_indices))), cmap=cmap, s=40, edgecolor='black', zorder=2)

# Add dots at the beginning and end of the plot
start_color = cmap(norm(0))
end_color = cmap(norm(len(sub_x)-1))
plt.scatter(sub_x[0], sub_y[0], color=start_color, s=100, zorder=3, edgecolor='black')  # Start point
plt.scatter(sub_x[-1], sub_y[-1], color=end_color, s=100, zorder=3, edgecolor='black')  # End point

# Add the special point from em_proj.xvg
plt.text(em_x, em_y, 'M', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)
plt.text(npt_x, npt_y, 'P', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)
plt.text(rst200_x, rst200_y, '1', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)
plt.text(rst50_x, rst50_y, '2', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)
plt.text(rst10_x, rst10_y, '3', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)
plt.text(rst0_x, rst0_y, '4', fontsize=10, ha='center', va='center', color='black', fontweight='bold', zorder=4)

# Set plot limits
plt.xlim(plot_range)
plt.ylim(plot_range)

# Set grid
plt.grid(True)

# Set aspect ratio
plt.gca().set_aspect('equal', adjustable='box')

# Set labels
plt.xlabel('PC1')
plt.ylabel('PC2')

# Set legend
plt.legend(['Trajectory'], loc='upper right')

# Add color bar
cbar = plt.colorbar(sc, ax=plt.gca(), orientation='vertical')
cbar.set_label('Normalized Color Index')

# Save the plot
plt.savefig('pca.png', dpi=300)
