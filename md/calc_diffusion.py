import numpy as np
from scipy.stats import linregress
import sys

def load_data(file_path):
    """ Load data from a file into a numpy array. """
    return np.loadtxt(file_path, comments=["@", "#"])

def perform_regressions(data):
    """ Perform linear regressions, excluding the first N points progressively. """
    n_points = len(data)
    results = []
    for i in range(n_points - 3):
        x = data[i:, 0]
        y = data[i:, 1]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        diff_coeff = (slope *1e3)/ 6.

        # Split the data in halves, calculate the diffusion coefficient for each half, and take the difference
        half = len(x) // 2
        half1 = x[:half], y[:half]
        half2 = x[half:], y[half:]
        slope1, _, _, _, std_err1 = linregress(half1[0], half1[1])
        slope2, _, _, _, std_err2 = linregress(half2[0], half2[1])
        error = np.abs(slope1 - slope2)*1e3/6.0

        results.append((diff_coeff, error))
    return results

def main():
    file_path = sys.argv[1]
    data = load_data(file_path)
    results = perform_regressions(data)
    
    for index, (slope, std_dev) in enumerate(results):
        print(f"{index} {slope:.5f} {std_dev:.5f}")

if __name__ == "__main__":
    main()