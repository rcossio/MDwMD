import numpy as np

def load_data(filepath):
    """ Load the COM data from an xvg file. """
    return np.loadtxt(filepath, comments=["@", "#"])

def calculate_squared_distances(data):
    """ Calculate the Euclidean distance from the initial COM position for each point. """
    initial_pos = data[0, 1:4]  # The initial x, y, z position is on the first row
    squared_distances = np.sum((data[:, 1:4] - initial_pos) ** 2, axis=1)
    return squared_distances

def main():
    filepath = 'com.xvg'  # Change this to the path of your xvg file
    data = load_data(filepath)
    sq_distances = calculate_squared_distances(data)

    # Print or save the results
    for time, sq_distance in zip(data[:, 0], sq_distances):
        print(f"{time} {sq_distance:.6f}")

if __name__ == "__main__":
    main()
