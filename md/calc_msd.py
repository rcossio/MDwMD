import numpy as np

def load_data(filepath):
    """ Load the COM data from an xvg file. """
    return np.loadtxt(filepath, comments=["@", "#"])

def estimate_msd(data,target_time):
    accum = 0
    n_intervals = len(data) - target_time
    for start_index in range(n_intervals):
        x_f = data[start_index+target_time, 1:4]
        x_i = data[start_index, 1:4]
        accum += np.linalg.norm(x_f-x_i)**2
    msd = accum / (n_intervals+1)
    return msd

def main():
    filepath = 'prod_com.xvg'
    data = load_data(filepath)
    times = data[:, 0]

    n_frames = len(data)
    msds = []
    for target_time in range(0, n_frames):
        msds.append(estimate_msd(data, target_time))
    
    for (time, msd) in zip(times, msds):
        print(f"{time} {msd}")

if __name__ == "__main__":
    main()
