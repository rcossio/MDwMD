import numpy as np
from scipy.stats import linregress
import sys
import matplotlib.pyplot as plt

def load_data(filepath):
    """Load the COM data from an xvg file."""
    try:
        return np.loadtxt(filepath, comments=["@", "#"])
    except Exception as e:
        print(f"Failed to load data: {e}")
        sys.exit(1)

def calculate_msd(data):
    """Calculate the Euclidean distance from the initial COM position for each point."""
    initial_pos = data[0, 1:4]  # The initial x, y, z position is on the first row
    msd = np.sum((data[:, 1:4] - initial_pos) ** 2, axis=1)
    return msd


def plotMSD(name, time, msd, slope, intercept):
    plt.figure(figsize=(10, 6))
    plt.plot(time, msd, label='MSD')
    plt.plot(time, slope*time + intercept, label='Fit')
    plt.title('Mean Squared Displacement')
    plt.grid(True)
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD (nm$^{2}$)')
    plt.legend()
    plt.savefig(name)
    plt.close()
    return None

def calculate_diffusion(com_data):
    n_points = len(com_data)
    diffusion_data = []
    for i in range(n_points - 3):
        time = com_data[i:, 0]
        print(com_data[i:,:].shape)
        msd = calculate_msd(com_data[i:,:])
        slope, intercept, r_value, p_value, std_err = linregress(time, msd)
        diff_coeff = (slope *1e6)/ 6.

        if i % 100 == 0:
            plotMSD(f"tmp/plot_{i}.png",time,msd, slope, intercept)

        # Split the data in halves, calculate the diffusion coefficient for each half, and take the difference
        half = len(time) // 2
        time_half1 = com_data[i:i+half,0]
        msd_half1  = calculate_msd(com_data[i:i+half,:])

        time_half2 = com_data[i+half:,0]
        msd_half2  = calculate_msd(com_data[i+half:,:])

        slope1, _, _, _, std_err1 = linregress(time_half1, msd_half1)
        slope2, _, _, _, std_err2 = linregress(time_half2, msd_half2)
        error = np.abs(slope1 - slope2)*1e6/6.0

        diffusion_data.append((com_data[i, 0], diff_coeff, error))

    return diffusion_data

def main():
    filepath = 'prod_com.xvg'
    com_data = load_data(filepath)
    
    diffusion_data = calculate_diffusion(com_data)
    time = [item[0] for item in diffusion_data]
    diffusion_coeffs = [item[1] for item in diffusion_data]
    diffusion_errors = [item[2] for item in diffusion_data]

    plt.figure(figsize=(10, 6))
    plt.plot(time, diffusion_coeffs, label='Diffusion Coefficient')
    plt.plot(time, diffusion_errors, label='Error')
    plt.title('Diffusion Coefficient from last frames')
    plt.grid(True)
    plt.xlabel('Time (ps)')
    plt.ylabel('D(um$^{2}$/s)')
    plt.legend()
    plt.savefig('diff_from_last.png')

if __name__ == "__main__":
    main()
