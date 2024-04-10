import pandas as pd
import matplotlib.pyplot as plt
import argparse
from io import StringIO

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process and plot data.')
parser.add_argument('-i', '--input', required=True, help='Input file name')
parser.add_argument('--cols', type=int, default=2, help='Number of columns to process, default is 2')
parser.add_argument('--moving_average', action='store_true', help='Plot moving average')

# Parse arguments
args = parser.parse_args()

filename = args.input
num_cols = args.cols

# Preprocess to skip lines starting with @ or #
with open(filename, 'r') as file:
    lines = file.readlines()
lines = [line for line in lines if not line.startswith(('@', '#'))]

# Convert lines back to a string to use with pandas
data_string = ''.join(lines)
data = pd.read_csv(StringIO(data_string), sep='\s+', header=None, usecols=range(num_cols))

# Assign column names
column_names = ['Index'] + [f'Value{i}' for i in range(1, num_cols)]
data.columns = column_names

# Define colors
mean_colors = ['darkgreen', 'darkblue', 'darkred', 'darkviolet', 'goldenrod']
colors =['seagreen', 'cornflowerblue', 'lightcoral', 'mediumpurple', 'khaki']

# Plotting
plt.figure(figsize=(10, 6))
for i in range(1, num_cols):
    # Calculate mean and standard deviation
    mean_value = data[column_names[i]].mean()
    std_deviation = data[column_names[i]].std()
    label_str = f'{column_names[i]} Mean(std) = {mean_value:.2f} ({std_deviation:.2f})'

    # Plot original values and label with stats
    plt.plot(data['Index'], data[column_names[i]], label=label_str, color=colors[i % len(colors)])

    # Check if moving average is requested
    if args.moving_average:
        data['Moving_Average_' + str(i)] = data[column_names[i]].rolling(window=10).mean()
        plt.plot(data['Index'], data['Moving_Average_' + str(i)], label=f'Moving Average {column_names[i]}', color=mean_colors[i % len(mean_colors)])

plt.xlabel('Index')
plt.ylabel('Values')
plt.legend()
plt.grid(True)
plt.savefig(filename.replace('.xvg', '.png'))
