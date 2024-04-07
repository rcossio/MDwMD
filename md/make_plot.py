import pandas as pd
import matplotlib.pyplot as plt
import sys

# Check for the input file
if len(sys.argv) < 2:
    print("Usage: python script.py <filename> [moving_average]")
    sys.exit(1)

filename = sys.argv[1]

# Preprocess to skip lines starting with @ or #
with open(filename, 'r') as file:
    lines = file.readlines()
lines = [line for line in lines if not line.startswith(('@', '#'))]

# Convert lines back to a string to use with pandas
from io import StringIO
data_string = ''.join(lines)
data = pd.read_csv(StringIO(data_string), sep='\\s+', header=None)

data.columns = ['Index', 'Value']

# Calculate the mean and standard deviation and print them
mean_value = data['Value'].mean()
std_deviation = data['Value'].std()

print(f"({filename}) Mean of Values: {mean_value}")
print(f"({filename}) Standard Deviation of Values: {std_deviation}")

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(data['Index'], data['Value'], label='Original Values')

# Check if moving average is requested
if len(sys.argv) > 2 and sys.argv[2] == "moving_average":
    data['Moving_Average'] = data['Value'].rolling(window=10).mean()
    plt.plot(data['Index'], data['Moving_Average'], label='Moving Average', color='red')

plt.xlabel('Index')
plt.ylabel('Values')
plt.legend()
plt.grid(True, zorder=5)
plt.savefig(filename.replace('.xvg', '.png'))
