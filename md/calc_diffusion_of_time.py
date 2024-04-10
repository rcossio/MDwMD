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

def main():
    filepath = 'com.xvg'
    com_data = load_data(filepath)
    
    num_chunks = 1  # Set the number of chunks here
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
    diffusion_coeff_of_t = 1e6*final_msd/6/times
    
    print(f"Average diffusion coefficient: {np.mean(diffusion_coeff_of_t[1:]):.2f} um^2/s")
    print(f"Standard deviation: {np.std(diffusion_coeff_of_t[1:]):.2f} um^2/s")
    print(f"Number of chunks: {len(chunked_msd)}")
    print(f"Number of frames per chunk: {frame_size}")

    plt.figure(figsize=(10, 6))
    plt.plot(times, diffusion_coeff_of_t, label='Diffusion Coefficient')
    plt.title('D(t)')
    plt.grid(True)
    plt.xlabel('Time (ps)')
    plt.ylabel('DiffusionMSD (um$^{2}$/s)')
    plt.legend()
    plt.savefig('diffusion_of_t.png')

if __name__ == "__main__":
    main()
