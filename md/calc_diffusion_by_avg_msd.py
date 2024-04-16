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

def calculate_diffusion(time, msd):
    """Calculate the diffusion coefficient from the given MSD data."""
    if len(time) != len(msd):
        print("Error: Time and MSD arrays must have the same length.")
        sys.exit(1)
    slope, intercept, r_value, p_value, std_err = linregress(time, msd)
    diff_coeff = (slope * 1e6) / 6.0  # Convert to appropriate units and calculate

    half = len(time) // 2
    half1 = time[:half], msd[:half]
    half2 = time[half:], msd[half:]
    slope1, _, _, _, std_err1 = linregress(half1[0], half1[1])
    slope2, _, _, _, std_err2 = linregress(half2[0], half2[1])
    error = np.abs(slope1 - slope2)*1e6/6.0

    return diff_coeff, error

def main():
    filepath = 'prod_com.xvg'
    com_data = load_data(filepath)
    
    num_chunks = 10  # Set the number of chunks here
    total_frames = len(com_data)
    frame_step = (total_frames - 1) // num_chunks
    frame_size = frame_step + 1
    
    chunked_com_data = [com_data[i:i + frame_size] for i in range(0, total_frames - frame_size + 1, frame_step)]
    chunked_msd = [calculate_msd(chunk) for chunk in chunked_com_data if len(chunk) == frame_size]

    if len(chunked_msd) != len(chunked_com_data):
        print("Error: Some chunks do not have the full number of frames.")
        sys.exit(1)

    times = chunked_com_data[0][:, 0]
    final_msd = np.mean(chunked_msd, axis=0)
    diffusion_coeff, diffusion_error = calculate_diffusion(times, final_msd)

    plt.figure(figsize=(10, 6))
    # Plot each chunk's MSD
    for msd in chunked_msd:
        plt.plot(times, msd, color='blue', linewidth=0.5, alpha=0.5)
    
    # Plot the average MSD
    plt.plot(times, final_msd, color='blue', linewidth=2, marker='o', label=f'Average MSD (D={diffusion_coeff:.1f}$\pm${diffusion_error:.1f} um$^{2}$/s)')
    plt.title('Mean Squared Displacements')
    plt.grid(True)
    plt.xlabel('Time (ps)')
    plt.ylabel('MSD (nm$^{2}$)')
    plt.legend()
    plt.savefig('final_msd.png')

if __name__ == "__main__":
    main()
